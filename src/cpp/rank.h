//------------------------------------------------------------------------------------
//
// rank.h - module with logic for k-mer rank tables
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef RANK_H
#define RANK_H

#include "kmer.h"

typedef Kmer KmerRank;

//------------------------------------------------------------------------------------

class KmerRankTable // a lookup table holding a rank for each k-mer
{
public:
   KmerRankTable(KmerLength kmerLen);

   virtual ~KmerRankTable() { delete[] rank; }

   void writeText(const String& textFilename) const;
   void writeBinary(const String& binaryFilename) const;

   const KmerLength k; // length of each k-mer
   KmerRank *rank;     // lookup table indexed by k-mer
};

KmerRankTable *readRankTable(const String& binaryFilename);

KmerRankTable *createRankTable(KmerLength k, const String& refGenFilename);

#endif
