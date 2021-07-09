//------------------------------------------------------------------------------------
//
// window.cpp - module supporting minimizer windows
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "window.h"

//------------------------------------------------------------------------------------

class Finder : public RankMinimizerFinder
{
public:
   Finder(const std::string& sequence, MinimizerWindowLength w,
          const KmerRankTable *rankTable, WindowVector& inWindowVector)
      : RankMinimizerFinder(sequence.c_str(), sequence.length(), w, rankTable),
	windowVector(inWindowVector) { }

   virtual ~Finder() { }

   virtual bool reportMinimizer(Minimizer minimizer, int startIndex, int windowID,
                                bool finalMinimizer)
   {
      windowVector.push_back(Window(minimizer, startIndex));
      return true;
   }

   WindowVector& windowVector;
};

//------------------------------------------------------------------------------------
// getWindows() partitions a sequence into consecutive windows and stores them in a
// vector

void getWindows(const std::string& sequence, MinimizerWindowLength w,
                const KmerRankTable *rankTable, WindowVector& windowVector)
{
   Finder finder(sequence, w, rankTable, windowVector);
   finder.find();
}
