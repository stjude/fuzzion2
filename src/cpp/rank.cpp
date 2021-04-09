//------------------------------------------------------------------------------------
//
// rank.cpp - module with logic for k-mer rank tables
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "rank.h"
#include "refgen.h"
#include "util.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

const uint32_t RANK_FILE_SIGNATURE_NOSWAP = 0x17D26E39;
const uint32_t RANK_FILE_SIGNATURE_SWAP   = 0x396ED217;

//------------------------------------------------------------------------------------

class KmerCount // holds #occurrences of one k-mer
{
public:
   KmerCount()
      : kmer(0), count(0) { }

   virtual ~KmerCount() { }

   inline void increment() { if (count < UINT32_MAX) count++; } // avoid overflow

   Kmer     kmer;  // 2-bit representation of k-mer
   uint32_t count; // #occurrences of the k-mer
};

typedef std::vector<KmerCount> KmerCountVector;

struct CompareKmerCounts // used to sort counts in ascending order
{
   bool operator()(const KmerCount& a, const KmerCount& b) const
   {
      return (a.count == b.count ? a.kmer < b.kmer : a.count < b.count);
   }
};

//------------------------------------------------------------------------------------

class KmerFinderCounter : public KmerFinder // finds and counts k-mers
{
public:
   KmerFinderCounter(const char *sequence, int sequenceLen, KmerLength kmerLen,
                     KmerCountVector *countVector)
      : KmerFinder(sequence, sequenceLen, kmerLen), cv(countVector) { }

   virtual ~KmerFinderCounter() { }

   virtual bool reportKmer(Kmer kmer, int startIndex)
   {
      Kmer revcomp = kmerReverseComplement(k, kmer);

      (*cv)[kmer].increment();
      (*cv)[revcomp].increment();

      return true;
   }

   KmerCountVector *cv;
};

//------------------------------------------------------------------------------------
// KmerRankTable::KmerRankTable() allocates but does not initialize the lookup table

KmerRankTable::KmerRankTable(KmerLength kmerLen)
   : k(kmerLen)
{
   if (k < 4 || k > MAX_KMER_LENGTH)
      throw std::runtime_error("unsupported k-mer length");

   rank = new KmerRank[numKmers(k)];
}

//------------------------------------------------------------------------------------
// KmerRankTable::writeText() writes a text file containing the k-mers and their ranks

void KmerRankTable::writeText(const std::string& textFilename) const
{
   std::ofstream outfile(textFilename.c_str());

   if (!outfile.is_open())
      throw std::runtime_error("unable to create " + textFilename);

   Kmer n = numKmers(k);

   for (Kmer kmer = 0; kmer < n; kmer++)
      outfile << kmerToString(k, kmer) << TAB << rank[kmer] << NEWLINE;

   outfile.close();
}

//------------------------------------------------------------------------------------
// KmerRankTable::writeBinary() writes a binary rank file

void KmerRankTable::writeBinary(const std::string& binaryFilename) const
{
   BinWriter writer;

   writer.open(binaryFilename);

   writer.writeUint32(RANK_FILE_SIGNATURE_NOSWAP);
   writer.writeUint8(k);

   // write the lookup table in four equal-sized chunks

   int quarter  = numKmers(k) >> 2;
   int numBytes = quarter * sizeof(KmerRank);

   for (int i = 0; i < 4; i++)
      writer.writeBuffer(&rank[i * quarter], numBytes);

   writer.close();
}

//------------------------------------------------------------------------------------
// readRankTable() returns a lookup table containing the ranks read from a binary rank
// file; it is the caller's obligation to de-allocate it

KmerRankTable *readRankTable(const std::string& binaryFilename)
{
   BinReader reader;

   reader.open(binaryFilename);

   uint32_t signature;
   KmerLength k;

   if (!reader.readUint32(signature) ||
       signature != RANK_FILE_SIGNATURE_NOSWAP &&
       signature != RANK_FILE_SIGNATURE_SWAP ||
       !reader.readUint8(k))
      throw std::runtime_error(binaryFilename + " is not a k-mer rank file");

   KmerRankTable *table = new KmerRankTable(k);

   // initialize the lookup table in four equal-sized chunks

   Kmer n = numKmers(k);

   int quarter  = n >> 2;
   int numBytes = quarter * sizeof(KmerRank);

   for (int i = 0; i < 4; i++)
      if (!reader.readBuffer(&table->rank[i * quarter], numBytes))
         throw std::runtime_error("truncated k-mer rank file " + binaryFilename);

   // make sure there are no additional bytes in the file
   uint8_t extraByte;
   if (reader.readUint8(extraByte))
      throw std::runtime_error("invalid k-mer rank file " + binaryFilename);

   reader.close();

   if (signature == RANK_FILE_SIGNATURE_SWAP) // fix the byte ordering
      for (Kmer kmer = 0; kmer < n; kmer++)
         swapBytes(&table->rank[kmer], sizeof(KmerRank));

   return table;
}

//------------------------------------------------------------------------------------
// createRankTable() returns a lookup table containing the ranks derived from a
// reference genome; it is the caller's obligation to de-allocate it

KmerRankTable *createRankTable(KmerLength k, const std::string& refGenFilename)
{
   Kmer n = numKmers(k);

   KmerCountVector *cv = new KmerCountVector(n);

   for (Kmer kmer = 0; kmer < n; kmer++)
      (*cv)[kmer].kmer = kmer;

   RefGenReader reader;

   reader.open(refGenFilename);

   for (int i = 0; i < reader.numref; i++) // for each reference
   {
      RefGenSeq *rgs = reader.getRefGenSeq(reader.refName[i], 1, 1000000000);

      KmerFinderCounter fc(rgs->seq, rgs->end, k, cv);
      fc.find();

      delete rgs;
   }

   reader.close();

   // sort the k-mers by ascending count
   std::sort(cv->begin(), cv->end(), CompareKmerCounts());

   KmerRankTable *table = new KmerRankTable(k);

   for (KmerRank rank = 0; rank < n; rank++)
      table->rank[(*cv)[rank].kmer] = rank;

   delete cv;

   return table;
}

//------------------------------------------------------------------------------------
// KmerRankInverter::KmerRankInverter() allocates and initializes the lookup table

KmerRankInverter::KmerRankInverter(const KmerRankTable *rankTable)
   : k(rankTable->k)
{
   Kmer n = numKmers(k);

   kmer = new Kmer[n];

   for (Kmer i = 0; i < n; i++)
      kmer[rankTable->rank[i]] = i;
}

//------------------------------------------------------------------------------------
// KmerRankInverter::getKmers() gets the k-mer sequence and its reverse complement for
// the given rank

void KmerRankInverter::getKmers(const std::string& rank,
                                std::string& kmer1, std::string& kmer2) const
{
   int ranklen = rank.length();

   if (ranklen == 0)
      throw std::runtime_error("empty rank string");

   for (int i = 0; i < ranklen; i++)
      if (!std::isdigit(rank[i]))
         throw std::runtime_error("invalid rank string " + rank);

   // convert rank to integer
   std::istringstream stream(rank);
   KmerRank  myRank;
   stream >> myRank;

   if (myRank >= numKmers(k))
      throw std::runtime_error("invalid rank string " + rank);

   Kmer myKmer1 = kmer[myRank];
   Kmer myKmer2 = kmerReverseComplement(k, myKmer1);

   kmer1 = kmerToString(k, myKmer1);
   kmer2 = kmerToString(k, myKmer2);
}
