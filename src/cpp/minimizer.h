//------------------------------------------------------------------------------------
//
// minimizer.h - module defining classes for finding minimizers
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "rank.h"

typedef Kmer KmerHash;
typedef KmerHash Minimizer;

typedef uint8_t MinimizerWindowLength;

inline int minimizerWindowID(int startIndex, MinimizerWindowLength windowLen)
{ return startIndex / windowLen; }

//------------------------------------------------------------------------------------

class MinimizerFinder : public KmerFinder // abstract class for finding minimizers
			                  // in a sequence
{
public:
   MinimizerFinder(const char *sequence, int sequenceLen, KmerLength kmerLen,
		   MinimizerWindowLength windowLen);

   virtual ~MinimizerFinder() { }

   virtual void find(); // finds and reports minimizers in the given sequence

   virtual bool reportKmer(Kmer kmer, int startIndex);

   // override this function to compute a hash value for a k-mer; a minimizer is the
   // smallest hash value for all k-mers in a window
   virtual KmerHash hash(Kmer kmer) = 0;

   // override this function to do whatever you want with each minimizer found;
   // return true to continue, false to discontinue the search
   virtual bool reportMinimizer(Minimizer minimizer, int startIndex,
		                int windowID, bool finalMinimizer) = 0;

   MinimizerWindowLength w;       // length of each window in bases
   int currentWindowID;           // ID of current window or (-1) if none
   Minimizer currentMinimizer;    // minimum hash found so far
   int currentStartIndex;         // start index of the k-mer having the minimum hash
};

//------------------------------------------------------------------------------------

class RankMinimizerFinder : public MinimizerFinder // abstract class for finding
			                           // rank minimizers in a sequence
{
public:
   RankMinimizerFinder(const char *sequence, int sequenceLen,
		       MinimizerWindowLength windowLen,
                       const KmerRankTable *rankTable)
      : MinimizerFinder(sequence, sequenceLen, rankTable->k, windowLen),
	table(rankTable) { }

   virtual ~RankMinimizerFinder() { }

   virtual KmerHash hash(Kmer kmer) { return table->rank[kmer]; }

   // override this function to do whatever you want with each rank minimizer found;
   // return true to continue, false to discontinue the search
   virtual bool reportMinimizer(Minimizer minimizer, int startIndex, int windowID,
		                bool finalMinimizer) = 0;

   const KmerRankTable *table; // lookup table holding a rank for each k-mer
};

//------------------------------------------------------------------------------------
#endif
