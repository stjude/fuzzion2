//------------------------------------------------------------------------------------
//
// rank.h - module with logic for k-mer rank tables
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef RANK_H
#define RANK_H

#include "kmer.h"

typedef Kmer KmerRank;

inline double kmerRankPercentile(KmerLength k, KmerRank rank)
{
   return (k <= MAX_KMER_LENGTH ?
           (static_cast<double>(rank) / numKmers(k)) * 100.0 : -1.0);
}

//------------------------------------------------------------------------------------

class KmerRankTable // a lookup table holding a rank for each k-mer
{
public:
   KmerRankTable(KmerLength kmerLen);

   virtual ~KmerRankTable() { delete[] rank; }

   void writeText(const std::string& textFilename) const;

   void writeBinary(const std::string& binaryFilename) const;

   KmerLength  k;    // length of each k-mer
   KmerRank   *rank; // lookup table indexed by k-mer
};

KmerRankTable *readRankTable(const std::string& binaryFilename);

KmerRankTable *createRankTable(KmerLength k, const std::string& refGenFilename);

//------------------------------------------------------------------------------------

class KmerRankInverter // a lookup table holding a k-mer for each rank
{
public:
   KmerRankInverter(const KmerRankTable *rankTable);

   virtual ~KmerRankInverter() { delete[] kmer; }

   void getKmers(const std::string& rank,
                 std::string& kmer1, std::string& kmer2) const;

   KmerLength  k;    // length of each k-mer
   Kmer       *kmer; // lookup table indexed by rank
};

//------------------------------------------------------------------------------------
#endif
