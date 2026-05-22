//------------------------------------------------------------------------------------
//
// match.h - module supporting pattern matching
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef MATCH_H
#define MATCH_H

#include "align.h"

const char ELLIPSIS_MARK = '.'; // indicates an ellipsis in a visualization

enum class MatchQual { QUALIFIED, DISQUALIFIED, UNDETERMINED };

//------------------------------------------------------------------------------------

class SingleMatch // match of a single read to a pattern
{
public:
   SingleMatch(const Align *inLeft, const Align *inMiddle, const Align *inRight,
               bool includeMiddleStats);
   SingleMatch(const SingleMatch& other);

   virtual ~SingleMatch() { delete left; delete middle; delete right; }

   bool sufficientAgreement(double minPercentAgreement) const;
   bool spanning(double minPercentAgreement, int minOverlap) const;
   bool betterThan(const SingleMatch& other) const;

   void setQual(const Pattern& pattern, double minPercentAgreement, int minOverlap);

   void getMiddleVis(const Pattern& pattern, int maxmidlen,
                     String& mvis1, String& mvis2) const;
   void getVis(const Pattern& pattern, int maxmidlen,
               String& vis1, String& vis2) const;

   Align *left;    // left or unsided alignment (nullptr if N/A)
   Align *middle;  // middle alignment          (nullptr if N/A)
   Align *right;   // right  alignment          (nullptr if N/A)

   int possible;   // possible #matching bases
   int matches;    // actual   #matching bases
   int score;      // score used to rank single matches

   int loverlap;   // length of left/unsided alignment
   int roverlap;   // length of right alignment

   MatchQual qual; // indicates whether the match qualifies as a hit

private:
   MatchQual leftRightQual(const Pattern& pattern);
};

typedef std::vector<SingleMatch *> SingleMatchVector;
typedef std::map<int, SingleMatchVector> SingleMatchMap; // key is pattern index

void selectBestSingleMatches(bool bestOverall, SingleMatchMap& mmap,
                             SingleMatchMap& bestMap);
void freeSingleMatchMap(SingleMatchMap& mmap);

//------------------------------------------------------------------------------------

class PairMatch // match of a read pair to a pattern
{
public:
   PairMatch(const Pattern& pattern, const SingleMatch *inMatch1,
             const SingleMatch *inMatch2);
   PairMatch(const PairMatch& other);

   virtual ~PairMatch() { delete match1; delete match2; }

   bool isPlausible(int maxInsert, int maxTrim) const
   { return (insertSize <= maxInsert && misalignment <= maxTrim); }

   bool betterThan(const PairMatch& other) const;

   void setQual(const Pattern& pattern, double minPercentAgreement, int minOverlap);

   void getVis(const Pattern& pattern, int maxmidlen, StringVector& vis) const;

   SingleMatch *match1; // match of read1 to pattern
   SingleMatch *match2; // match of read2 to pattern

   int possible;        // possible #matching bases
   int matches;         // actual   #matching bases
   int score;           // score used to rank pair matches

   int insertSize;      // insert size of read pair aligned to pattern
   int misalignment;    // #bases read2 is aligned ahead of read1 (normally is zero)

   MatchQual qual;      // indicates whether the match qualifies as a hit
};

typedef std::map<int, PairMatch *> PairMatchMap; // maps pattern index to best match

void selectBestPairMatch(PairMatchMap& pairMap, PairMatchMap& bestMap);
void freePairMatchMap(PairMatchMap& pairMap);

//------------------------------------------------------------------------------------

int getSingleMatches(const String& readStr, MinimizerWindowLength w,
                     const KmerRankTable *rankTable, Minimizer maxMinimizer,
                     const PatternMap& pmap, const BoolVector& inPmap,
                     const PatternVector& pvector, double minPercentAgreement,
                     int minOverlap, Seq *& readSeq, SingleMatchMap& mmap,
                     const SingleMatchMap *mateMap=nullptr);

int getPairMatches(const String& readStr1, const String& readStr2,
                   MinimizerWindowLength w, const KmerRankTable *rankTable,
                   Minimizer maxMinimizer, const PatternMap& pmap,
                   const BoolVector& inPmap, const PatternVector& pvector,
                   double minPercentAgreement, int minOverlap, int maxInsert,
                   int maxTrim, Seq *& readSeq1, Seq *& readSeq2,
                   PairMatchMap& pairMap);

#endif
