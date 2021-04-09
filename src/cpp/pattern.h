//------------------------------------------------------------------------------------
//
// pattern.h - module supporting pattern sequences
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef PATTERN_H
#define PATTERN_H

#include "minimizer.h"
#include <unordered_map>
#include <vector>

//------------------------------------------------------------------------------------

class Pattern
{
public:
   Pattern(const std::string& inName, const std::string& inSequence);

   virtual ~Pattern() { }

   std::string name;            // pattern name
   std::string sequence;        // sequence w/o brackets or braces
   std::string displaySequence; // sequence w/  brackets or braces
   int  leftBases;              // #bases in left  side of sequence
   int  rightBases;             // #bases in right side of sequence
   bool hasBraces;              // true = braces, false = brackets
};

typedef std::vector<Pattern *> PatternVector;

//------------------------------------------------------------------------------------

class Location // indicates a location within a pattern sequence
{
public:
   Location(int inIndex, int inOffset)
      : index(inIndex), offset(inOffset) { }

   virtual ~Location() { }

   int index;  // index of pattern in a PatternVector
   int offset; // offset within the pattern sequence
};

typedef std::vector<Location> LocationVector;

// this maps a minimizer to locations of the minimizer in patterns
typedef std::unordered_map<Minimizer, LocationVector> PatternMap;

//------------------------------------------------------------------------------------

class Candidate // candidate for best alignment of a read to a pattern
{
public:
   Candidate(const Location& inLocation, int inMatches)
      : location(inLocation), matches(inMatches) { }

   virtual ~Candidate() { }

   Location location; // location within a pattern sequence
   int matches;       // number of bases matching the read
};

typedef std::vector<Candidate> CandidateVector;

//------------------------------------------------------------------------------------

PatternVector *readPatterns(const std::string& filename);

PatternMap *createPatternMap(const PatternVector *patternVector,
                             MinimizerWindowLength w,
			     const KmerRankTable *rankTable,
			     Minimizer maxMinimizer);

int lengthOfLCS(const std::string& strA, int offsetA, int lenA,
		const std::string& strB, int offsetB, int lenB);

CandidateVector *findCandidates(const std::string& sequence,
	 	                PatternMap *patternMap,
                                const PatternVector *patternVector,
			        MinimizerWindowLength w,
			        const KmerRankTable *rankTable,
			        Minimizer maxMinimizer,
			        int minMinimizers, int minMatches);

//------------------------------------------------------------------------------------
#endif
