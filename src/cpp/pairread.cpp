//------------------------------------------------------------------------------------
//
// pairread.cpp - module relating to paired reads
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "pairread.h"
#include <cstring>
#include <stdexcept>

//------------------------------------------------------------------------------------
// InputReader::open() opens the first pair reader

void InputReader::open()
{
   if (readerVector.size() == 0)
      throw std::runtime_error("no input files");

   if (current >= 0) // a pair reader is currently open
      close();       // close it

   readerVector[0]->open();
   current = 0;
}

//------------------------------------------------------------------------------------
// InputReader::getNextPair() gets the next pair of reads in the currently open file;
// if EOF is reached on the currently open file, the next file in the series is
// opened; if there are no more files in the series, false is returned to indicate
// EOF of the series

bool InputReader::getNextPair(std::string& name1, std::string& sequence1,
                              std::string& name2, std::string& sequence2)
{
   int numReaders = readerVector.size();

   if (current < 0 || current >= numReaders)
      return false;

   while (!readerVector[current]->getNextPair(name1, sequence1, name2, sequence2))
      if (current + 1 < numReaders)
      {
         readerVector[current]->close();
	 readerVector[++current]->open();
      }
      else
      {
         close();
	 return false;
      }

   return true;
}

//------------------------------------------------------------------------------------
// InputReader::close() closes the currently open pair reader

void InputReader::close()
{
   if (current >= 0 && current < readerVector.size())
      readerVector[current]->close();

   current = -1;
}

//------------------------------------------------------------------------------------
// namesMatch() returns true if the given read names match

bool namesMatch(const std::string& name1, const std::string& name2)
{
   int len = name1.length();
   if (len != name2.length())
      return false; // names have unequal lengths

   // see if one name ends with "1" and the other name ends with "2"
   int last = len - 1;

   if (last > 0 && (name1[last] == '1' && name2[last] == '2' ||
                    name1[last] == '2' && name2[last] == '1'))
      return (std::strncmp(name1.c_str(), name2.c_str(), last) == 0);

   return (name1 == name2);
}
