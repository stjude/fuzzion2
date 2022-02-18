//------------------------------------------------------------------------------------
//
// fastq.cpp - module for reading FASTQ files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "fastq.h"
#include <stdexcept>

//------------------------------------------------------------------------------------
// FastqReader::open() opens a file for input; if the filename ends with ".gz", the
// file is assumed to be gzipped and it is uncompressed

void FastqReader::open()
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
      throw std::runtime_error("attempt to read from unopened FASTQ file " +
                               filename);

   const int MAX_LINE_LEN = 100000;
   typedef char Line[MAX_LINE_LEN];

   Line nameLine, seqLine, plusLine, qualLine;

   if (!std::fgets(nameLine, MAX_LINE_LEN, f))
      return false; // reached end-of-file

   if (!std::fgets(seqLine,  MAX_LINE_LEN, f) ||
       !std::fgets(plusLine, MAX_LINE_LEN, f) ||
       !std::fgets(qualLine, MAX_LINE_LEN, f) ||
       nameLine[0] != '@' || plusLine[0] != '+')
      throw std::runtime_error("unexpected format in FASTQ file " + filename);

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
         throw std::runtime_error("mismatched read names " + name1 + " and " + name2 +
                                  " in "  + reader1.filename +
				  " and " + reader2.filename);
   else
      throw std::runtime_error("different number of reads in " + reader1.filename +
			       " and " + reader2.filename);
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
      throw std::runtime_error("odd number of reads in " + reader.filename);

   if (namesMatch(name1, name2))
      return true;
   else
      throw std::runtime_error("mismatched read names " + name1 + " and " + name2 +
                               " in " + reader.filename);
}

//------------------------------------------------------------------------------------
// isFastqFile() returns true if the named file is a FASTQ file; if so, the names of
// the first two reads in the file are obtained and interleaved is set to true if
// these read names match one another

bool isFastqFile(const std::string& filename, std::string& name1, std::string& name2,
                 bool& interleaved)
{
   FastqReader reader(filename);

   try
   {
      reader.open();
   }
   catch (const std::runtime_error& error) { return false; }

   bool isFastq = false;

   try
   {
      std::string sequence1, sequence2;

      if (reader.getNext(name1, sequence1))
      {
         if (reader.getNext(name2, sequence2))
            interleaved = namesMatch(name1, name2);
	 else
	 {
            name2 = "_NONE_"; // only one read in the file
	    interleaved = false;
	 }

	 isFastq = true;
      }
   }
   catch (const std::runtime_error& error) { }

   reader.close();

   return isFastq;
}
