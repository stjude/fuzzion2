//------------------------------------------------------------------------------------
//
// ubam.cpp - module for reading unaligned Bam files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "ubam.h"
#include <stdexcept>

//------------------------------------------------------------------------------------
// UbamPairReader::getNextPair() gets the next read pair and returns true, or returns
// false when end-of-file is reached

bool UbamPairReader::getNextPair(std::string& name1, std::string& sequence1,
                                 std::string& name2, std::string& sequence2)
{
   BamRead read1, read2;

   if (!reader.getNext(read1))
      return false; // reached EOF

   if (!reader.getNext(read2))
      throw std::runtime_error("odd number of reads in" + filename);

   name1 = read1.constName();
   name2 = read2.constName();

   if (!namesMatch(name1, name2))
      throw std::runtime_error("mismatched read names " + name1 + " and " + name2 +
                               " in " + filename);

   std::string *seq1 = read1.sequence();
   std::string *seq2 = read2.sequence();

   sequence1 = *seq1;
   sequence2 = *seq2;

   delete seq1;
   delete seq2;

   return true;
}

//------------------------------------------------------------------------------------
// isUbamFile() returns true if the named file is an unaligned Bam file

bool isUbamFile(const std::string& filename)
{
   BamReader reader;

   try
   {
      reader.open(filename);
   }
   catch (const std::runtime_error& error) { return false; }

   bool isUbam = false;

   try
   {
      BamRead read1, read2;

      if (reader.getNext(read1) && reader.getNext(read2))
      {
         std::string name1 = read1.constName();
	 std::string name2 = read2.constName();

	 if (namesMatch(name1, name2))
            isUbam = true;
      }
   }
   catch (const std::runtime_error& error) { }

   reader.close();

   return isUbam;
}
