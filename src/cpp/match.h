//------------------------------------------------------------------------------------
//
// match.h - module supporting pattern matching
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
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
      : Location(location), length(inLength), matchingBases(inMatchingBases) { }

   virtual ~Candidate() { }

   int length;        // read length
   int matchingBases; // number of matching bases; it is zero for an unmatched mate
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

   int  matchingBases() const { return c1.matchingBases + c2.matchingBases; }
   int  possible()   const;
   int  insertSize() const;

   bool leftmost1()  const;
   bool rightmost2() const;

   Candidate c1, c2; // if a single-read match, c2 represents the unmatched mate
};

typedef std::vector<Match> MatchVector;

//------------------------------------------------------------------------------------

void getMatches(const std::string& sequence1, const std::string& sequence2,
                const PatternVector *patternVector, PatternMap *patternMap,
		MinimizerWindowLength w, const KmerRankTable *rankTable,
		Minimizer maxMinimizer, double minBases, int minMins, int maxInsert,
		int maxTrim, bool bestOverall, bool findSingle,
		MatchVector& matchVector);

bool validOverlaps(const std::string& sequence1, const std::string& sequence2,
                   const PatternVector *patternVector, double minBases,
		   int minOverlap, const Match& match);

//------------------------------------------------------------------------------------
#endif
