//------------------------------------------------------------------------------------
//
// ubam.cpp - module for reading unaligned Bam files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "ubam.h"
#include <stdexcept>

//------------------------------------------------------------------------------------
// UbamPairReader::open() opens the first of the unaligned Bam files

void UbamPairReader::open(const StringVector& inFilename)
{
   filename = inFilename;

   int numFiles = filename.size();
   if (numFiles == 0)
      throw std::runtime_error("no file names specified");

   reader.open(filename[0]);

   currentFile = 0;
}

//------------------------------------------------------------------------------------
// UbamPairReader::getNextPair() gets the next read pair and returns true, or returns
// false when end-of-file on the last file is reached

bool UbamPairReader::getNextPair(std::string& name1, std::string& sequence1,
                                 std::string& name2, std::string& sequence2)
{
   BamRead read1, read2;

   int numFiles = filename.size();

   while (!reader.getNext(read1) && ++currentFile < numFiles)
   {
      // reached end-of-file on the current file
      reader.close(); // close the current file
      reader.open(filename[currentFile]); // open the next file
   }

   if (currentFile == numFiles)
      return false; // reached EOF on the last file

   if (!reader.getNext(read2))
      throw std::runtime_error("mismatched number of reads");

   name1 = read1.constName();
   name2 = read2.constName();

   if (!namesMatch(name1, name2))
      throw std::runtime_error("mismatched read names " + name1 + " and " + name2);

   std::string *seq1 = read1.sequence();
   std::string *seq2 = read2.sequence();

   sequence1 = *seq1;
   sequence2 = *seq2;

   delete seq1;
   delete seq2;

   return true;
}
