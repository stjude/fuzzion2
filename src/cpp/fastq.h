//------------------------------------------------------------------------------------
//
// fastq.h - module for reading FASTQ files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef FASTQ_H
#define FASTQ_H

#include "pairread.h"
#include <cstdio>

//------------------------------------------------------------------------------------

class FastqReader // for reading a single FASTQ file
{
public:
   FastqReader()
      : f(NULL), isPipe(false) { }

   virtual ~FastqReader() { }

   void open(const std::string& filename);

   bool getNext(std::string& name, std::string& sequence);

   void close();

   std::FILE *f;        // is non-NULL when the file is open
   bool       isPipe;   // is true when reading from a pipe
};

//------------------------------------------------------------------------------------

class FastqPairReader : public PairReader
{
public:
   FastqPairReader()
      : reader1(), reader2() { }

   virtual ~FastqPairReader() { }

   void open(const std::string& filename1, const std::string& filename2)
   { reader1.open(filename1); reader2.open(filename2); }

   bool getNextPair(std::string& name1, std::string& sequence1,
                    std::string& name2, std::string& sequence2);

   void close() { reader1.close(); reader2.close(); }

   FastqReader reader1, reader2;
};

//------------------------------------------------------------------------------------

class InterleavedFastqPairReader : public PairReader
{
public:
   InterleavedFastqPairReader()
      : reader() { }

   virtual ~InterleavedFastqPairReader() { }

   void open(const std::string& filename) { reader.open(filename); }

   bool getNextPair(std::string& name1, std::string& sequence1,
                    std::string& name2, std::string& sequence2);

   void close() { reader.close(); }

   FastqReader reader;
};

//------------------------------------------------------------------------------------
#endif
