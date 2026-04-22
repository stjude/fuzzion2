//------------------------------------------------------------------------------------
//
// window.h - module supporting minimizer windows
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef WINDOW_H
#define WINDOW_H

#include "minimizer.h"

//------------------------------------------------------------------------------------

class Window // represents a window of a sequence
{
public:
   Window(const Minimizer inMinimizer, const int inOffset)
      : minimizer(inMinimizer), offset(inOffset) { }

   virtual ~Window() { }

   const Minimizer minimizer; // window minimizer
   const int       offset;    // offset of first base of minimizing k-mer
};

typedef std::vector<Window> WindowVector;

//------------------------------------------------------------------------------------

void getWindows(const char *sequence, int sequenceLen, MinimizerWindowLength w,
                const KmerRankTable *rankTable, WindowVector& windowVector);

#endif
