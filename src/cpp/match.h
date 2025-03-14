//------------------------------------------------------------------------------------
//
// match.h - module supporting pattern matching
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef MATCH_H
#define MATCH_H

#include "pattern.h"
#include <map>

//------------------------------------------------------------------------------------

class Candidate : public Location // represents a match of a read to a pattern
{
public:
   Candidate(const Location& location, int inLength, int inMatchingBases)
      : Location(location), length(inLength), matchingBases(inMatchingBases),
	leftOverlap(0), leftMatching(0), rightOverlap(0), rightMatching(0),
	junctionSpanning(false) { }

   virtual ~Candidate() { }

   int  length;                      // read length
   int  matchingBases;               // #matching bases (zero for an unmatched mate)

   void setLeftRight(const std::string& sequence, const PatternVector *patternVector);
   int  leftOverlap,  leftMatching;  // #overlapping and matching bases on the left
   int  rightOverlap, rightMatching; // #overlapping and matching bases on the right

   void setJunctionSpanning(double minBases, int minOverlap);
   bool junctionSpanning;            // true if the read spans the junction
};

typedef std::vector<Candidate> CandidateVector;

typedef std::map<int, CandidateVector> CandidateMap; // key is pattern index

//------------------------------------------------------------------------------------

class Match // represents a match of a read pair to a pattern
{
public:
   Match(const Candidate& candidate1, const Candidate& candidate2)
      : c1(candidate1), c2(candidate2) { }

   virtual ~Match() { }

   int matchingBases() const { return c1.matchingBases + c2.matchingBases; }
   int possible()      const;
   int numSpanning()   const;
   int insertSize()    const;

   bool validOverlaps(const std::string& sequence1, const std::string& sequence2,
                      const PatternVector *patternVector, double minBases,
		      int minOverlap);

   Candidate c1, c2;
};

typedef std::vector<Match> MatchVector;

//------------------------------------------------------------------------------------

void getMatches(const std::string& sequence1, const std::string& sequence2,
                const PatternVector *patternVector, PatternMap *patternMap,
		MinimizerWindowLength w, const KmerRankTable *rankTable,
		Minimizer maxMinimizer, double minBases, int minMins, int maxInsert,
		int maxTrim, bool bestOverall, bool findSingle,
		MatchVector& matchVector);

#endif
