//------------------------------------------------------------------------------------
//
// ubam.h - module for reading unaligned Bam files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef UBAM_H
#define UBAM_H

#include "bamread.h"
#include "pairread.h"

//------------------------------------------------------------------------------------

class UbamPairReader : public PairReader
{
public:
   UbamPairReader(const std::string& inFilename)
      : filename(inFilename), reader() { }

   virtual ~UbamPairReader() { }

   void open()  { reader.open(filename); }

   bool getNextPair(std::string& name1, std::string& sequence1,
                    std::string& name2, std::string& sequence2);

   void close() { reader.close(); }

   std::string filename;
   BamReader  reader;
};

//------------------------------------------------------------------------------------

bool isUbamFile(const std::string& filename);

#endif
