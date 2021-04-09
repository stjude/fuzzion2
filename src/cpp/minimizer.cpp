//------------------------------------------------------------------------------------
//
// minimizer.cpp - module defining classes for finding minimizers
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "minimizer.h"
#include <stdexcept>

//------------------------------------------------------------------------------------
// MinimizerFinder::MinimizerFinder() checks for a valid window length and initializes
// data members

MinimizerFinder::MinimizerFinder(const char *sequence, int sequenceLen,
		                 KmerLength kmerLen, MinimizerWindowLength windowLen)
   : KmerFinder(sequence, sequenceLen, kmerLen), w(windowLen), currentWindowID(-1),
     currentMinimizer(0), currentStartIndex(-1)
{
   if (w < 1)
      throw std::runtime_error("invalid minimizer window length");
}

//------------------------------------------------------------------------------------
// MinimizerFinder::find() searches for minimizers in the sequence

void MinimizerFinder::find()
{
   KmerFinder::find();

   if (currentWindowID >= 0) // report final minimizer
      reportMinimizer(currentMinimizer, currentStartIndex, currentWindowID, true);
}

//------------------------------------------------------------------------------------
// MinimizerFinder::reportKmer() examines each k-mer found in the sequence;
// minimizers are constructed and reported

bool MinimizerFinder::reportKmer(Kmer kmer, int startIndex)
{
   KmerHash kmerHash = hash(kmer);

   int windowID = minimizerWindowID(startIndex, w);

   if (windowID == currentWindowID)
   {
      if (kmerHash < currentMinimizer)
      {
         currentMinimizer  = kmerHash;
	 currentStartIndex = startIndex;
      }
   }
   else // this is the first k-mer of a new window
   {
      // report minimizer for the previous window
      if (currentWindowID >= 0 &&
          !reportMinimizer(currentMinimizer, currentStartIndex, currentWindowID,
                           false))
      {
         // discontinue the search
	 currentWindowID = -1;
	 return false;
      }

      currentWindowID   = windowID;
      currentMinimizer  = kmerHash;
      currentStartIndex = startIndex;
   }

   return true;
}
