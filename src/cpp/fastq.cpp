//------------------------------------------------------------------------------------
//
// fastq.cpp - module for reading FASTQ files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "fastq.h"
#include <stdexcept>

//------------------------------------------------------------------------------------
// FastqReader::open() opens the named file for input; if the filename ends with
// ".gz", the file is assumed to be gzipped and it is uncompressed

void FastqReader::open(const std::string& filename)
{
   if (f)
      throw std::runtime_error("FASTQ file already open");

   int len = filename.length();

   if (len > 3 && filename.substr(len - 3, 3) == ".gz") // gzipped file
   {
      std::string command = "gunzip -c " + filename;    // we will uncompress it

      f = popen(command.c_str(), "r");
      isPipe = true;
   }
   else // file is uncompressed
   {
      f = std::fopen(filename.c_str(), "r");
      isPipe = false;
   }

   if (!f)
      throw std::runtime_error("unable to open " + filename);
}

//------------------------------------------------------------------------------------
// FastqReader::getNext() gets the next read in the FASTQ file and returns true, or
// returns false when end-of-file is reached

bool FastqReader::getNext(std::string& name, std::string& sequence)
{
   if (!f)
      throw std::runtime_error("attempt to read from unopened FASTQ file");

   const int MAX_LINE_LEN = 100000;
   typedef char Line[MAX_LINE_LEN];

   Line nameLine, seqLine, plusLine, qualLine;

   if (!std::fgets(nameLine, MAX_LINE_LEN, f))
      return false; // reached end-of-file

   if (!std::fgets(seqLine,  MAX_LINE_LEN, f) ||
       !std::fgets(plusLine, MAX_LINE_LEN, f) ||
       !std::fgets(qualLine, MAX_LINE_LEN, f) ||
       nameLine[0] != '@' || plusLine[0] != '+')
      throw std::runtime_error("unexpected format in FASTQ file");

   int i = 1;
   while (nameLine[i] && !isspace(nameLine[i]))
      i++;

   nameLine[i] = 0;

   i = 0;
   while (seqLine[i] && !isspace(seqLine[i]))
      i++;

   seqLine[i] = 0;

   name     = &nameLine[1];
   sequence = seqLine;

   return true;
}

//------------------------------------------------------------------------------------
// FastqReader::close() closes the FASTQ file

void FastqReader::close()
{
   if (!f) // no file is open
      return;

   if (isPipe)
      pclose(f);
   else
      std::fclose(f);

   f      = NULL;
   isPipe = false;
}

//------------------------------------------------------------------------------------
// FastqPairReader::getNextPair() gets the next read pair and returns true, or returns
// false when end-of-file is reached

bool FastqPairReader::getNextPair(std::string& name1, std::string& sequence1,
                                  std::string& name2, std::string& sequence2)
{
   bool gotRead1 = reader1.getNext(name1, sequence1);
   bool gotRead2 = reader2.getNext(name2, sequence2);

   if (!gotRead1 && !gotRead2)
      return false; // reached EOF on both files

   if (gotRead1 && gotRead2)
      if (namesMatch(name1, name2))
         return true;
      else
         throw std::runtime_error("mismatched read names " + name1 + " and " + name2);
   else
      throw std::runtime_error("mismatched number of reads");
}

//------------------------------------------------------------------------------------
// InterleavedFastqPairReader::getNextPair() gets the next read pair and returns true,
// or returns false when end-of-file is reached

bool InterleavedFastqPairReader::getNextPair(
                                  std::string& name1, std::string& sequence1,
                                  std::string& name2, std::string& sequence2)
{
   if (!reader.getNext(name1, sequence1))
      return false; // reached EOF

   if (!reader.getNext(name2, sequence2))
      throw std::runtime_error("mismatched number of reads");

   if (namesMatch(name1, name2))
      return true;
   else
      throw std::runtime_error("mismatched read names " + name1 + " and " + name2);
}
