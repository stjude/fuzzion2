//------------------------------------------------------------------------------------
//
// match.cpp - module supporting pattern matching
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "match.h"
#include "window.h"
#include <algorithm>
#include <cmath>

//------------------------------------------------------------------------------------

struct LocationCompare // for sorting a LocationVector
{
   bool operator()(const Location& a, const Location& b) const
   {
      // sort by ascending index, then by ascending offset
      return (a.index == b.index ? a.offset < b.offset : a.index < b.index);
   }
};

//------------------------------------------------------------------------------------

struct MatchCompare // for sorting a MatchVector
{
   bool operator()(const Match& a, const Match& b) const
   {
      // sort by descending number of matching bases, then by ascending insert size,
      // then by ascending pattern index

      int key1 = a.matchingBases() - b.matchingBases();

      if (key1 != 0)
         return (key1 > 0);

      int key2 = a.insertSize() - b.insertSize();

      if (key2 != 0)
         return (key2 < 0);

      return (a.c1.index < b.c1.index);
   }
};

//------------------------------------------------------------------------------------
// Match::possible() returns the maximum possible number of matching bases

int Match::possible() const
{
   if (c1.matchingBases == 0)      // c1 is an unmatched mate
      return c2.length;
   else if (c2.matchingBases == 0) // c2 is an unmatched mate
      return c1.length;
   else // matching read-pair
      return c1.length + c2.length;
}

//------------------------------------------------------------------------------------
// Match::numSpanning() returns the number of junction-spanning reads (0, 1 or 2)

int Match::numSpanning() const
{
   if (c1.junctionSpanning && c2.junctionSpanning)
      return 2;
   else if (c1.junctionSpanning || c2.junctionSpanning)
      return 1;
   else
      return 0;
}

//------------------------------------------------------------------------------------
// Match::insertSize() returns the insert size of this match

int Match::insertSize() const
{
   if (c1.matchingBases == 0)       // c1 is an unmatched mate
      return c2.length;
   else if (c2.matchingBases == 0)  // c2 is an unmatched mate
      return c1.length;
   else if (c1.offset <= c2.offset) // c1 is aligned ahead of c2
      return std::max(c1.length, c2.offset - c1.offset + c2.length);
   else                             // c2 is aligned ahead of c1
      return std::max(c2.length, c1.offset - c2.offset + c1.length);
}

//------------------------------------------------------------------------------------
// computeMinMatches() returns the minimum number of matching bases for the given
// sequence length

static int computeMinMatches(int seqlen, double minBases)
{
   return std::ceil((minBases / 100) * seqlen);
}

//------------------------------------------------------------------------------------
// lengthOfLCS() returns the length of a longest common subsequence of substrings of
// two strings; it is a measure of similarity of these substrings

static int lengthOfLCS(const std::string& strA, int offsetA, int lenA,
		       const std::string& strB, int offsetB, int lenB)
{
   if (lenA <= 0 || lenB <= 0)
      return 0;

   const char *cstrA = strA.c_str();
   const char *cstrB = strB.c_str();

   const char *a = &cstrA[offsetA];
   const char *b = &cstrB[offsetB];

   int *previous = new int[lenB + 1]();
   int *current  = new int[lenB + 1]();

   for (int i = 0; i < lenA; i++)
   {
      for (int j = 0; j < lenB; j++)
         if (a[i] == b[j])
            current[j + 1] = previous[j] + 1;
         else
            current[j + 1] = std::max(current[j], previous[j + 1]);

      int *temp = previous;
      previous  = current;
      current   = temp;
   }

   int matches = previous[lenB];

   delete[] previous;
   delete[] current;

   return matches;
}

//------------------------------------------------------------------------------------
// getLocations() extracts minimizers from a sequence, looks them up in a pattern map,
// and obtains locations of the minimizers within patterns; only eligible patterns are
// considered; passing NULL in the last parameter means all patterns are eligible

static void getLocations(const std::string& sequence, PatternMap *patternMap,
                         MinimizerWindowLength w, const KmerRankTable *rankTable,
			 Minimizer maxMinimizer, LocationVector& locationVector,
			 BoolVector *eligiblePattern=NULL)
{
   bool allPatternsEligible = (eligiblePattern == NULL);

   WindowVector windowVector;

   getWindows(sequence, w, rankTable, windowVector);

   int numWindows = windowVector.size();

   for (int i = 0; i < numWindows; i++)
   {
      Minimizer minimizer = windowVector[i].minimizer;
      if (minimizer > maxMinimizer)
         continue; // ignore common minimizer

      PatternMap::iterator ppos = patternMap->find(minimizer);
      if (ppos == patternMap->end())
         continue; // minimizer is not in map

      const LocationVector& location = ppos->second;
      int numLocations = location.size();

      for (int j = 0; j < numLocations; j++)
      {
         int index = location[j].index;

	 if (allPatternsEligible || (*eligiblePattern)[index])
	 {
            // determine starting offset of pattern's matching substring
	    int offset = std::max(0, location[j].offset - windowVector[i].offset);

	    locationVector.push_back(Location(index, offset));
	 }
      }
   }

   if (locationVector.size() > 1)
      std::sort(locationVector.begin(), locationVector.end(), LocationCompare());
}

//------------------------------------------------------------------------------------
// getCandidates() identifies candidate matches of a sequence and stores them in a
// map; only eligible patterns are considered; passing NULL in the last parameter
// means all patterns are eligible

static void getCandidates(const std::string& sequence,
                          const PatternVector *patternVector, PatternMap *patternMap,
			  MinimizerWindowLength w, const KmerRankTable *rankTable,
			  Minimizer maxMinimizer, double minBases, int minMins,
			  CandidateMap& cmap, BoolVector *eligiblePattern=NULL)
{
   LocationVector locationVector;

   getLocations(sequence, patternMap, w, rankTable, maxMinimizer, locationVector,
                eligiblePattern);

   int seqlen     = sequence.length();
   int minMatches = computeMinMatches(seqlen, minBases);

   int i = 0, numLocations = locationVector.size();

   while (i < numLocations)
   {
      const Location& location = locationVector[i];

      int count = 1;

      while (++i < numLocations &&
             locationVector[i].index  == location.index &&
	     locationVector[i].offset == location.offset)
         count++;

      if (count < minMins)
         continue; // not enough matching minimizers

      const std::string& psequence = (*patternVector)[location.index].sequence;

      int pseqlen = psequence.length();
      int pcmplen = std::min(seqlen, pseqlen - location.offset);

      int matchingBases = lengthOfLCS(sequence, 0, seqlen,
                                      psequence, location.offset, pcmplen);

      if (matchingBases < minMatches)
         continue; // not enough matching bases

      // found a candidate; put it in the map
      CandidateMap::iterator cpos = cmap.find(location.index);

      if (cpos == cmap.end())
      {
         cmap.insert(std::make_pair(location.index, CandidateVector()));
	 cpos = cmap.find(location.index);
      }

      CandidateVector& candidateVector = cpos->second;

      candidateVector.push_back(Candidate(location, seqlen, matchingBases));
   }
}

//------------------------------------------------------------------------------------
// getBestPair() finds the best read-pair match for each pattern if bestOverall is
// false, or the best overall read-pair match, across all patterns, if bestOverall is
// true

static void getBestPair(CandidateMap& cmap1, CandidateMap& cmap2, int maxInsert,
                        int maxTrim, bool bestOverall, MatchVector& matchVector)
{
   int best = 0;

   for (CandidateMap::iterator cpos1 = cmap1.begin(); cpos1 != cmap1.end(); ++cpos1)
   {
      CandidateMap::iterator cpos2 = cmap2.find(cpos1->first);
      if (cpos2 == cmap2.end())
         continue;

      const CandidateVector& cv1 = cpos1->second; int n1 = cv1.size();
      const CandidateVector& cv2 = cpos2->second; int n2 = cv2.size();

      for (int i = 0; i < n1; i++)
         for (int j = 0; j < n2; j++)
	 {
            Match match(cv1[i], cv2[j]);

	    if (match.insertSize() > maxInsert ||
                match.c1.offset - match.c2.offset > maxTrim)
               continue; // insert size too large or second read aligned too far
                         // ahead of first read

	    if (best == matchVector.size())
               matchVector.push_back(match); // save first match, best so far
	    else
	    {
               int matchingBases     = match.matchingBases();
	       int mostMatchingBases = matchVector[best].matchingBases();

	       if (matchingBases >  mostMatchingBases ||
                   matchingBases == mostMatchingBases &&
		   match.insertSize() < matchVector[best].insertSize())
                  matchVector[best] = match; // overwrite best match with new best
	    }
	 }

      if (!bestOverall && best < matchVector.size())
         best++; // advance for next pattern
   }
}

//------------------------------------------------------------------------------------
// getBestSingle() finds the best single-read match for each pattern if bestOverall is
// false, or the best overall single-read match, across all patterns, if bestOverall
// is true

static void getBestSingle(CandidateMap& cmap, bool bestOverall, bool firstRead,
                          int mateLength, MatchVector& matchVector)
{
   int best = (bestOverall ? 0 : matchVector.size());

   for (CandidateMap::iterator cpos = cmap.begin(); cpos != cmap.end(); ++cpos)
   {
      const CandidateVector& cv = cpos->second; int n = cv.size();

      for (int i = 0; i < n; i++)
         if (best == matchVector.size() ||
             cv[i].matchingBases > matchVector[best].matchingBases())
	 {
            Location location(cv[i].index, cv[i].offset);
	    Candidate mate(location, mateLength, 0); // unmatched mate

	    const Candidate& c1 = (firstRead ? cv[i] : mate);
	    const Candidate& c2 = (firstRead ? mate  : cv[i]);

	    Match match(c1, c2);

	    if (best == matchVector.size())
               matchVector.push_back(match); // save first match, best so far
	    else
               matchVector[best] = match;    // overwrite best match with new best
	 }

      if (!bestOverall && best < matchVector.size())
         best++; // advance for next pattern
   }
}

//------------------------------------------------------------------------------------
// getMatches() finds pattern matches for the given read pair and sorts them by
// descending number of matching bases; if bestOverall is true, only the best overall
// match is identified; if findSingle is true, single-read matches are sought when
// there are no read-pair matches

void getMatches(const std::string& sequence1, const std::string& sequence2,
                const PatternVector *patternVector, PatternMap *patternMap,
		MinimizerWindowLength w, const KmerRankTable *rankTable,
		Minimizer maxMinimizer, double minBases, int minMins, int maxInsert,
		int maxTrim, bool bestOverall, bool findSingle,
		MatchVector& matchVector)
{
   CandidateMap cmap1, cmap2;

   getCandidates(sequence1, patternVector, patternMap, w, rankTable, maxMinimizer,
                 minBases, minMins, cmap1);

   if (cmap1.size() == 0 && !findSingle)
      return;

   std::string revcomp = stringReverseComplement(sequence2);

   if (findSingle)
      getCandidates(revcomp, patternVector, patternMap, w, rankTable, maxMinimizer,
                    minBases, minMins, cmap2);
   else // we only need to find candidates in certain patterns
   {
      BoolVector eligiblePattern(patternVector->size(), false);

      for (CandidateMap::iterator cpos1 = cmap1.begin();
           cpos1 != cmap1.end(); ++cpos1)
         eligiblePattern[cpos1->first] = true;

      getCandidates(revcomp, patternVector, patternMap, w, rankTable, maxMinimizer,
                    minBases, minMins, cmap2, &eligiblePattern);
   }

   if (cmap1.size() > 0 && cmap2.size() > 0)
      getBestPair(cmap1, cmap2, maxInsert, maxTrim, bestOverall, matchVector);

   if (matchVector.size() == 0 && findSingle)
   {
      // didn't find a matching read-pair; look for single-read matches
      if (cmap1.size() > 0)
         getBestSingle(cmap1, bestOverall, true,  sequence2.length(), matchVector);

      if (cmap2.size() > 0)
         getBestSingle(cmap2, bestOverall, false, sequence1.length(), matchVector);
   }

   if (matchVector.size() > 1)
      std::sort(matchVector.begin(), matchVector.end(), MatchCompare());
}

//------------------------------------------------------------------------------------
// overlapLeft() returns true if the given sequence aligned at the specified offset
// overlaps the left side of the pattern

static bool overlapLeft(const std::string& sequence, const Pattern& pattern,
                        int offset, double minBases, int minOverlap)
{
   int maxPatternLen = pattern.sequence.length() - offset;
   int patternLen    = std::min((int)sequence.length(), maxPatternLen);
   int overlapBases  = std::min(patternLen, pattern.leftBases - offset);

   if (overlapBases < minOverlap)
      return false; // insufficient number of overlapping bases

   if (overlapBases == patternLen)
      return true;  // full overlap; we have already verified agreement of bases

   // verify agreement of bases for partial overlap
   int lcslen =
      lengthOfLCS(sequence, 0, overlapBases, pattern.sequence, offset, overlapBases);

   return (lcslen >= computeMinMatches(overlapBases, minBases));
}

//------------------------------------------------------------------------------------
// overlapRight() returns true if the given sequence aligned at the specified offset
// overlaps the right side of the pattern

static bool overlapRight(const std::string& sequence, const Pattern& pattern,
                         int offset, double minBases, int minOverlap)
{
   int maxPatternLen = pattern.sequence.length() - offset;
   int patternLen    = std::min((int)sequence.length(), maxPatternLen);
   int overlapBases  = patternLen + std::min(0, pattern.rightBases - maxPatternLen);

   if (overlapBases < minOverlap)
      return false; // insufficient number of overlapping bases

   if (overlapBases == patternLen)
      return true;  // full overlap; we have already verified agreement of bases

   // verify agreement of bases for partial overlap
   int lcslen =
      lengthOfLCS(sequence, sequence.length() - overlapBases, overlapBases,
                  pattern.sequence, offset + patternLen - overlapBases, overlapBases);

   return (lcslen >= computeMinMatches(overlapBases, minBases));
}

//------------------------------------------------------------------------------------
// Match::validOverlaps() returns true if the given match satisfies the overlap
// requirements

bool Match::validOverlaps(const std::string& sequence1, const std::string& sequence2,
                          const PatternVector *patternVector, double minBases,
                          int minOverlap)
{
   const Pattern& pattern = (*patternVector)[c1.index];

   bool left1 = false, right1 = false;
   bool left2 = false, right2 = false;

   if (c1.matchingBases > 0)
   {
      left1  = overlapLeft (sequence1, pattern, c1.offset, minBases, minOverlap);
      right1 = overlapRight(sequence1, pattern, c1.offset, minBases, minOverlap);
   }

   if (c2.matchingBases > 0)
   {
      left2  = overlapLeft (sequence2, pattern, c2.offset, minBases, minOverlap);
      right2 = overlapRight(sequence2, pattern, c2.offset, minBases, minOverlap);
   }

   c1.junctionSpanning = (left1 && right1);
   c2.junctionSpanning = (left2 && right2);

   if (c1.junctionSpanning || c2.junctionSpanning)
      return true;

   if (pattern.hasBraces)
      return false; // a junction-spanning read is required in this case

   return (left1 && right2 || left2 && right1);
}
