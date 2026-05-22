//------------------------------------------------------------------------------------
//
// window.cpp - module supporting minimizer windows
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "window.h"

//------------------------------------------------------------------------------------

class Finder : public RankMinimizerFinder
{
public:
   Finder(const char *sequence, const int sequenceLen, const MinimizerWindowLength w,
          const KmerRankTable *rankTable, WindowVector& inWindowVector)
      : RankMinimizerFinder(sequence, sequenceLen, w, rankTable),
        windowVector(inWindowVector) { }

   virtual ~Finder() { }

   virtual bool reportMinimizer(const Minimizer minimizer,
                                const int startIndex, const int windowID,
                                const bool finalMinimizer) override
   {
      windowVector.emplace_back(minimizer, startIndex);
      return true;
   }

   WindowVector& windowVector;
};

//------------------------------------------------------------------------------------
// getWindows() partitions a sequence into consecutive windows and stores them in a
// vector

void getWindows(const char *sequence, const int sequenceLen,
                const MinimizerWindowLength w, const KmerRankTable *rankTable,
                WindowVector& windowVector)
{
   Finder finder(sequence, sequenceLen, w, rankTable, windowVector);
   finder.find();
}
