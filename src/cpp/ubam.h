//------------------------------------------------------------------------------------
//
// ubam.h - module for reading unaligned Bam files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef UBAM_H
#define UBAM_H

#include "bamread.h"
#include "pairread.h"
#include "util.h"

//------------------------------------------------------------------------------------

class UbamPairReader : public PairReader
{
public:
   UbamPairReader()
      : reader(), filename(), currentFile(-1) { }

   virtual ~UbamPairReader() { }

   void open(const StringVector& inFilename);

   bool getNextPair(std::string& name1, std::string& sequence1,
                    std::string& name2, std::string& sequence2);

   void close() { reader.close(); }

   BamReader    reader;
   StringVector filename;
   int          currentFile; // index into filename of file currently being read
};

//------------------------------------------------------------------------------------
#endif
