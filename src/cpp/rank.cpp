//------------------------------------------------------------------------------------
//
// rank.cpp - module for k-mer rank tables
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "rank.h"
#include "refgen.h"

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
   KmerFinderCounter(const char *sequence, const int sequenceLen,
                     const KmerLength kmerLen, KmerCountVector *countVector)
      : KmerFinder(sequence, sequenceLen, kmerLen), cv(countVector) { }

   virtual ~KmerFinderCounter() { }

   virtual bool reportKmer(const Kmer kmer, const int startIndex) override
   {
      const Kmer revcomp = kmerReverseComplement(k, kmer);

      (*cv)[kmer].increment();
      (*cv)[revcomp].increment();

      return true;
   }

   KmerCountVector *cv;
};

//------------------------------------------------------------------------------------
// KmerRankTable::KmerRankTable() allocates but does not initialize the lookup table

KmerRankTable::KmerRankTable(const KmerLength kmerLen)
   : k(kmerLen)
{
   if (k < 4 || k > MAX_KMER_LENGTH)
      throw Error("unsupported k-mer length");

   rank = new KmerRank[numKmers(k)];
}

//------------------------------------------------------------------------------------
// KmerRankTable::writeText() writes a text file containing the k-mers and their ranks

void KmerRankTable::writeText(const String& textFilename) const
{
   std::ofstream outfile(textFilename.c_str());

   if (!outfile.is_open())
      throw Error("unable to create " + textFilename);

   const Kmer n = numKmers(k);

   for (Kmer kmer = 0; kmer < n; kmer++)
      outfile << kmerToString(k, kmer) << TAB << rank[kmer] << NEWLINE;

   outfile.close();
}

//------------------------------------------------------------------------------------
// KmerRankTable::writeBinary() writes a binary rank file

void KmerRankTable::writeBinary(const String& binaryFilename) const
{
   BinWriter writer;

   writer.open(binaryFilename);

   writer.writeUint32(RANK_FILE_SIGNATURE_NOSWAP);
   writer.writeUint8(k);

   // write the lookup table in four equal-sized chunks

   const int quarter  = numKmers(k) >> 2;
   const int numBytes = quarter * sizeof(KmerRank);

   for (int i = 0; i < 4; i++)
      writer.writeBuffer(&rank[i * quarter], numBytes);

   writer.close();
}

//------------------------------------------------------------------------------------
// readRankTable() returns a lookup table containing the ranks read from a binary rank
// file; it is the caller's obligation to de-allocate it

KmerRankTable *readRankTable(const String& binaryFilename)
{
   BinReader reader;

   reader.open(binaryFilename);

   uint32_t signature;
   KmerLength k;

   if (!reader.readUint32(signature) ||
       signature != RANK_FILE_SIGNATURE_NOSWAP &&
       signature != RANK_FILE_SIGNATURE_SWAP ||
       !reader.readUint8(k))
      throw Error(binaryFilename + " is not a k-mer rank file");

   KmerRankTable *table = new KmerRankTable(k);

   // initialize the lookup table in four equal-sized chunks

   const Kmer n = numKmers(k);

   const int quarter  = n >> 2;
   const int numBytes = quarter * sizeof(KmerRank);

   for (int i = 0; i < 4; i++)
      if (!reader.readBuffer(&table->rank[i * quarter], numBytes))
         throw Error("truncated k-mer rank file " + binaryFilename);

   // make sure there are no additional bytes in the file
   uint8_t extraByte;
   if (reader.readUint8(extraByte))
      throw Error("invalid k-mer rank file " + binaryFilename);

   reader.close();

   if (signature == RANK_FILE_SIGNATURE_SWAP) // fix the byte ordering
      for (Kmer kmer = 0; kmer < n; kmer++)
         swapBytes(&table->rank[kmer], sizeof(KmerRank));

   return table;
}

//------------------------------------------------------------------------------------
// createRankTable() returns a lookup table containing the ranks derived from a
// reference genome; it is the caller's obligation to de-allocate it

KmerRankTable *createRankTable(const KmerLength k, const String& refGenFilename)
{
   const Kmer n = numKmers(k);

   KmerCountVector *cv = new KmerCountVector(n);

   for (Kmer kmer = 0; kmer < n; kmer++)
      (*cv)[kmer].kmer = kmer;

   RefGenReader reader;

   reader.open(refGenFilename);

   for (int i = 0; i < reader.numref; i++) // for each reference
   {
      const RefGenSeq *rgs = reader.getRefGenSeq(reader.refName[i], 1, 1000000000);

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
