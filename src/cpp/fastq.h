//------------------------------------------------------------------------------------
//
// fastq.h - module for reading FASTQ files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
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
   FastqReader(const std::string& inFilename)
      : filename(inFilename), f(NULL), isPipe(false) { }

   virtual ~FastqReader() { }

   void open();

   bool getNext(std::string& name, std::string& sequence);

   void close();

   std::string  filename; // name of file to read
   std::FILE   *f;        // is non-NULL when the file is open
   bool         isPipe;   // is true when reading from a pipe
};

//------------------------------------------------------------------------------------

class FastqPairReader : public PairReader
{
public:
   FastqPairReader(const std::string& filename1, const std::string& filename2)
      : reader1(filename1), reader2(filename2) { }

   virtual ~FastqPairReader() { }

   void open()  { reader1.open(); reader2.open(); }

   bool getNextPair(std::string& name1, std::string& sequence1,
                    std::string& name2, std::string& sequence2);

   void close() { reader1.close(); reader2.close(); }

   FastqReader reader1, reader2;
};

//------------------------------------------------------------------------------------

class InterleavedFastqPairReader : public PairReader
{
public:
   InterleavedFastqPairReader(const std::string& filename)
      : reader(filename) { }

   virtual ~InterleavedFastqPairReader() { }

   void open()  { reader.open(); }

   bool getNextPair(std::string& name1, std::string& sequence1,
                    std::string& name2, std::string& sequence2);

   void close() { reader.close(); }

   FastqReader reader;
};

//------------------------------------------------------------------------------------

bool isFastqFile(const std::string& filename, std::string& name1, std::string& name2,
                 bool& interleaved);

#endif
