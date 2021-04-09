//------------------------------------------------------------------------------------
//
// pattern.cpp - module supporting pattern sequences
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "pattern.h"
#include "util.h"
#include <algorithm>
#include <fstream>
#include <stdexcept>

//------------------------------------------------------------------------------------

class Window // represents a window of a sequence
{
public:
   Window(Minimizer inMinimizer, int inOffset)
      : minimizer(inMinimizer), offset(inOffset) { }

   virtual ~Window() { }

   Minimizer minimizer; // window minimizer
   int       offset;    // offset of first base of minimizing k-mer
};

typedef std::vector<Window> WindowVector;

//------------------------------------------------------------------------------------

class Finder : public RankMinimizerFinder
{
public:
   Finder(const std::string& sequence, MinimizerWindowLength w,
          const KmerRankTable *rankTable, WindowVector *inWindowVector)
      : RankMinimizerFinder(sequence.c_str(), sequence.length(), w, rankTable),
	windowVector(inWindowVector) { }

   virtual ~Finder() { }

   virtual bool reportMinimizer(Minimizer minimizer, int startIndex, int windowID,
                                bool finalMinimizer)
   {
      windowVector->push_back(Window(minimizer, startIndex));
      return true;
   }

   WindowVector *windowVector;
};

//------------------------------------------------------------------------------------
// getWindows() returns a vector of minimizers extracted from the given sequence; it
// is the caller's obligation to de-allocate it

static WindowVector *getWindows(const std::string& sequence, MinimizerWindowLength w,
                                const KmerRankTable *rankTable)
{
   WindowVector *windowVector = new WindowVector();

   Finder finder(sequence, w, rankTable, windowVector);
   finder.find();

   return windowVector;
}

//------------------------------------------------------------------------------------
// Pattern::Pattern() constructs an object representing the named pattern

Pattern::Pattern(const std::string& inName, const std::string& inSequence)
   : name(inName), displaySequence(inSequence)
{
   bool foundDelimiters = false;
   int  delim1, delim2;

   if ((delim1 = inSequence.find(']')) != std::string::npos &&
       (delim2 = inSequence.find('[')) != std::string::npos)
   {
      foundDelimiters = true;
      hasBraces       = false;
   }
   else if ((delim1 = inSequence.find('}')) != std::string::npos &&
            (delim2 = inSequence.find('{')) != std::string::npos)
   {
      foundDelimiters = true;
      hasBraces       = true;
   }

   if (!foundDelimiters || delim1 == 0 || delim2 == inSequence.length() - 1 ||
       delim1 > delim2)
      throw std::runtime_error("invalid pattern " + inSequence);

   int middleBases = delim2 - delim1 - 1;

   leftBases  = delim1;
   rightBases = inSequence.length() - 1 - delim2;

   sequence   = inSequence.substr(0, leftBases) +
                inSequence.substr(delim1 + 1, middleBases) +
		inSequence.substr(delim2 + 1, rightBases);
}

//------------------------------------------------------------------------------------
// readPatterns() returns a vector of patterns read from the named file; it is the
// caller's obligation to de-allocate it

PatternVector *readPatterns(const std::string& filename)
{
   PatternVector *patternVector = new PatternVector();

   std::ifstream infile(filename.c_str());
   if (!infile.is_open())
      throw std::runtime_error("unable to open " + filename);

   std::string line;

   while (std::getline(infile, line))
   {
      StringVector column;

      if (splitString(line, column) != 2)
         throw std::runtime_error("unexpected #columns in " + filename);

      patternVector->push_back(new Pattern(column[0], column[1]));
   }

   infile.close();

   return patternVector;
}

//------------------------------------------------------------------------------------
// createPatternMap() returns a map constructed from the given pattern vector; it is
// the caller's obligation to de-allocate it

PatternMap *createPatternMap(const PatternVector *patternVector,
                             MinimizerWindowLength w,
			     const KmerRankTable *rankTable,
			     Minimizer maxMinimizer)
{
   PatternMap *patternMap = new PatternMap();

   int numPatterns = patternVector->size();

   for (int i = 0; i < numPatterns; i++)
   {
      const Pattern *pattern = (*patternVector)[i];

      WindowVector *windowVector = getWindows(pattern->sequence, w, rankTable);

      int numWindows = windowVector->size();

      for (int j = 0; j < numWindows; j++)
      {
         Minimizer minimizer = (*windowVector)[j].minimizer;
	 if (minimizer > maxMinimizer)
            continue; // don't put common minimizer in map

	 PatternMap::iterator mpos = patternMap->find(minimizer);

	 if (mpos == patternMap->end())
	 {
            patternMap->insert(std::make_pair(minimizer, LocationVector()));
	    mpos = patternMap->find(minimizer);
	 }

	 LocationVector& locationVector = mpos->second;

	 locationVector.push_back(Location(i, (*windowVector)[j].offset));
      }

      delete windowVector;
   }

   return patternMap;
}

//------------------------------------------------------------------------------------
// lengthOfLCS() returns the length of a longest common subsequence of substrings of
// two strings; it is a measure of similarity of these substrings

int lengthOfLCS(const std::string& strA, int offsetA, int lenA,
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

   int matches = current[lenB];

   delete[] previous;
   delete[] current;

   return matches;
}

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
// getHits() looks for minimizer hits in the pattern map and stores them in the given
// LocationVector

static void getHits(const std::string& sequence, PatternMap *patternMap,
                    const PatternVector *patternVector, MinimizerWindowLength w,
		    const KmerRankTable *rankTable, Minimizer maxMinimizer,
		    LocationVector& hit)
{
   WindowVector *windowVector = getWindows(sequence, w, rankTable);

   int numWindows = windowVector->size();

   for (int i = 0; i < numWindows; i++)
   {
      Minimizer minimizer = (*windowVector)[i].minimizer;
      if (minimizer > maxMinimizer)
         continue; // ignore common minimizer

      PatternMap::iterator mpos = patternMap->find(minimizer);
      if (mpos == patternMap->end())
         continue; // minimizer is not in map

      const LocationVector& location = mpos->second;
      int numLocations = location.size();

      for (int j = 0; j < numLocations; j++)
         if ((*patternVector)[location[j].index]) // this is an eligible pattern
	 {
            // determine the starting offset of the pattern's matching substring
            int offset = std::max(0, location[j].offset - (*windowVector)[i].offset);

	    hit.push_back(Location(location[j].index, offset));
	 }
   }

   delete windowVector;

   if (hit.size() > 1)
      std::sort(hit.begin(), hit.end(), LocationCompare());
}

//------------------------------------------------------------------------------------
// findCandidates() finds candidates for best pattern match of the given sequence;
// only eligible patterns, i.e., those with non-NULL entries in patternVector, are
// considered; a vector of candidates is returned; it is the caller's obligation to
// de-allocate it

CandidateVector *findCandidates(const std::string& sequence, PatternMap *patternMap,
                                const PatternVector *patternVector,
			        MinimizerWindowLength w,
				const KmerRankTable *rankTable,
			        Minimizer maxMinimizer, int minMinimizers,
				int minMatches)
{
   CandidateVector *candidate = new CandidateVector();
   LocationVector   hit;

   getHits(sequence, patternMap, patternVector, w, rankTable, maxMinimizer, hit);

   int i = 0, numHits = hit.size();

   while (i < numHits)
   {
      Location currentLocation = hit[i];
      int      currentCount    = 1;

      while (++i < numHits && hit[i].index == currentLocation.index &&
             hit[i].offset == currentLocation.offset)
         currentCount++;

      if (currentCount < minMinimizers)
         continue; // not enough matching uncommon minimizers

      const std::string& patternSequence =
         (*patternVector)[currentLocation.index]->sequence;

      int patternLen = std::min(sequence.length(),
                                patternSequence.length() - currentLocation.offset);

      int matches = lengthOfLCS(sequence, 0, sequence.length(),
                                patternSequence, currentLocation.offset, patternLen);

      if (matches >= minMatches) // found enough base matches
         candidate->push_back(Candidate(currentLocation, matches));
   }

   return candidate;
}
