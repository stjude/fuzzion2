//------------------------------------------------------------------------------------
//
// match.cpp - module supporting pattern matching
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "match.h"
#include "window.h"
#include <algorithm>
#include <cmath>

typedef std::vector<bool> BoolVector;

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

      int key2 = a.offsetDiff() - b.offsetDiff();

      if (key2 != 0)
         return (key2 < 0);

      return (a.c1.index < b.c1.index);
   }
};

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
// considered; passing NULL in place of a BoolVector means all patterns are eligible

static void getLocations(const std::string& sequence, PatternMap *patternMap,
                         MinimizerWindowLength w, const KmerRankTable *rankTable,
			 Minimizer maxMinimizer, BoolVector *eligiblePattern,
			 LocationVector& locationVector)
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
// getCandidates() identifies candidate matches of a sequence and stores them in a map

static void getCandidates(const std::string& sequence,
                          const PatternVector *patternVector, double minBases,
			  int minMins, const LocationVector& locationVector,
			  CandidateMap& candidateMap)
{
   int minMatches = computeMinMatches(sequence.length(), minBases);

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
         continue; // not enough matching uncommon minimizers

      const std::string& patternSequence = (*patternVector)[location.index].sequence;

      int patternLen =
         std::min(sequence.length(), patternSequence.length() - location.offset);

      int matchingBases =
         lengthOfLCS(sequence, 0, sequence.length(), patternSequence, location.offset,
                     patternLen);

      if (matchingBases < minMatches)
         continue; // not enough matching bases

      // found a candidate; put it in the map
      CandidateMap::iterator cpos = candidateMap.find(location.index);

      if (cpos == candidateMap.end())
      {
         candidateMap.insert(std::make_pair(location.index, CandidateVector()));
	 cpos = candidateMap.find(location.index);
      }

      CandidateVector& candidateVector = cpos->second;

      candidateVector.push_back(Candidate(location, matchingBases));
   }
}

//------------------------------------------------------------------------------------
// getAllCandidates() identifies candidate matches for each read of a read pair and
// puts them in a separate map for each read; true is returned if there is at least
// one candidate match for each read

static bool getAllCandidates(const std::string& sequence1,
                             const std::string& sequence2,
			     const PatternVector *patternVector,
			     PatternMap *patternMap,
			     MinimizerWindowLength w, const KmerRankTable *rankTable,
			     Minimizer maxMinimizer, double minBases, int minMins,
			     CandidateMap& cmap1, CandidateMap& cmap2)
{
   LocationVector locationVector;

   // get candidate matches for the first read

   getLocations(sequence1, patternMap, w, rankTable, maxMinimizer, NULL,
                locationVector);

   getCandidates(sequence1, patternVector, minBases, minMins, locationVector, cmap1);

   if (cmap1.size() == 0)
      return false; // didn't find any

   locationVector.clear();

   // get candidate matches for the reverse complement of the second read; only
   // patterns with candidate matches of the first read are eligible

   std::string revcomp = stringReverseComplement(sequence2);

   BoolVector eligiblePattern(patternVector->size(), false);

   for (CandidateMap::iterator cpos = cmap1.begin(); cpos != cmap1.end(); ++cpos)
      eligiblePattern[cpos->first] = true;

   getLocations(revcomp, patternMap, w, rankTable, maxMinimizer, &eligiblePattern,
                locationVector);

   getCandidates(revcomp, patternVector, minBases, minMins, locationVector, cmap2);

   return (cmap2.size() > 0);
}

//------------------------------------------------------------------------------------
// pairupCandidates() finds the best match for each pattern if bestOverall is false;
// otherwise, it finds only the best overall match across all patterns

static void pairupCandidates(CandidateMap& cmap1, CandidateMap& cmap2,
                             int maxOffsetDiff, bool bestOverall,
			     MatchVector& matchVector)
{
   Location  location(0, 0);
   Candidate candidate(location, 0);
   Match bestMatch(candidate, candidate); // will hold best match found so far
   bool foundMatch = false;               // nothing found yet

   for (CandidateMap::iterator cpos1 = cmap1.begin(); cpos1 != cmap1.end(); ++cpos1)
   {
      CandidateMap::iterator cpos2 = cmap2.find(cpos1->first);
      if (cpos2 == cmap2.end())
         continue; // nothing for this pattern

      const CandidateVector& cv1 = cpos1->second; int n1 = cv1.size();
      const CandidateVector& cv2 = cpos2->second; int n2 = cv2.size();

      for (int i = 0; i < n1; i++)
         for (int j = 0; j < n2; j++)
	 {
            Match match(cv1[i], cv2[j]);

	    int offsetDiff    = match.offsetDiff();
	    int matchingBases = match.matchingBases();

	    if (offsetDiff >= 0 && offsetDiff <= maxOffsetDiff && (!foundMatch ||
                matchingBases >  bestMatch.matchingBases() ||
		matchingBases == bestMatch.matchingBases() &&
		offsetDiff < bestMatch.offsetDiff()))
	    {
               foundMatch = true;
	       bestMatch  = match;
	    }
	 }

      if (!bestOverall && foundMatch)
      {
         matchVector.push_back(bestMatch); // best match for this pattern
	 foundMatch = false;               // reset for next pattern
      }
   }

   if (foundMatch)
      matchVector.push_back(bestMatch);    // best overall match

   if (matchVector.size() > 1)
      std::sort(matchVector.begin(), matchVector.end(), MatchCompare());
}

//------------------------------------------------------------------------------------
// getMatches() finds pattern matches for the given read pair and sorts them by
// descending number of matching bases; if bestOverall is true, only the best overall
// match is identified

void getMatches(const std::string& sequence1, const std::string& sequence2,
                const PatternVector *patternVector, PatternMap *patternMap,
		MinimizerWindowLength w, const KmerRankTable *rankTable,
		Minimizer maxMinimizer, double minBases, int minMins, int maxInsert,
		bool bestOverall, MatchVector& matchVector)
{
   CandidateMap cmap1, cmap2;

   if (getAllCandidates(sequence1, sequence2, patternVector, patternMap, w, rankTable,
                        maxMinimizer, minBases, minMins, cmap1, cmap2))
      pairupCandidates(cmap1, cmap2, maxInsert - sequence2.length(), bestOverall,
                       matchVector);
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
// validOverlaps() returns true if the given match satisfies the overlap requirements

bool validOverlaps(const std::string& sequence1, const std::string& sequence2,
                   const PatternVector *patternVector, double minBases,
                   int minOverlap, const Match& match)
{
   const Pattern& pattern = (*patternVector)[match.c1.index];

   int offset1 = match.c1.offset;
   int offset2 = match.c2.offset;

   if (pattern.hasBraces)
      return (overlapLeft (sequence1, pattern, offset1, minBases, minOverlap) &&
              overlapRight(sequence1, pattern, offset1, minBases, minOverlap) ||
	      overlapLeft (sequence2, pattern, offset2, minBases, minOverlap) &&
	      overlapRight(sequence2, pattern, offset2, minBases, minOverlap));
   else // pattern has brackets
      return (overlapLeft (sequence1, pattern, offset1, minBases, minOverlap) &&
              overlapRight(sequence2, pattern, offset2, minBases, minOverlap));
}
