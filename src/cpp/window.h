//------------------------------------------------------------------------------------
//
// window.h - module supporting minimizer windows
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef WINDOW_H
#define WINDOW_H

#include "minimizer.h"
#include <vector>

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

void getWindows(const std::string& sequence, MinimizerWindowLength w,
                const KmerRankTable *rankTable, WindowVector& windowVector);

//------------------------------------------------------------------------------------
#endif
