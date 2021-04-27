//------------------------------------------------------------------------------------
//
// fuzzort.cpp - this program reads fuzzion2 hits from stdin, sorts them, and writes
//               them to stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "util.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

const std::string VERSION_ID  = "fuzzort v1.0.0, copyright 2021 "
                                "St. Jude Children's Research Hospital";

const std::string PATTERN     = "pattern";
const std::string READ_PAIRS  = "read pairs";

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stdout

void showUsage(const char *progname)
{
   std::cout
      << VERSION_ID << NEWLINE << NEWLINE
      << "Usage: " << progname << " < fuzzion2_hits > sorted_hits" << NEWLINE;
}

//------------------------------------------------------------------------------------

class Hit
{
public:
   Hit(const std::string& inPatternName, const std::string& inRead1Name,
       const std::string& inLines)
      : patternName(inPatternName), read1Name(inRead1Name), lines(inLines) { }

   virtual ~Hit() { }

   std::string patternName, read1Name, lines;
};

typedef std::vector<Hit *> HitVector;

struct HitCompare
{
   bool operator()(Hit* const& a, Hit* const& b) const
   {
      // order by ascending patternName, then by ascending read1Name
      return (a->patternName == b->patternName ?
              (a->read1Name   < b->read1Name) :
	      (a->patternName < b->patternName));
   }
};

void sortHits(HitVector& hitVector)
{
   std::sort(hitVector.begin(), hitVector.end(), HitCompare());
}

//------------------------------------------------------------------------------------
// readHits() reads the hits from stdin and stores them in the vector; the number of
// read pairs is summed and returned to the caller

uint64_t readHits(HitVector& hitVector)
{
   const std::string EOL(1, NEWLINE); // string containing a newline character

   uint64_t numReadPairs = 0;

   std::string line1, line2, line3;

   while (std::getline(std::cin, line1))
   {
      StringVector line1_col, line1a_col, line2_col, line3_col;

      if (splitString(line1, line1_col) != 2)
         throw std::runtime_error("invalid input");

      if (line1_col[1] == READ_PAIRS)
      {
	 std::istringstream stream(line1_col[0]);
	 uint64_t  count;
	 stream >> count;
	 numReadPairs += count;
	 continue;
      }

      // verify that we have encountered a valid hit
      if (splitString(line1_col[0], line1a_col, ' ') != 3 ||
          line1a_col[0] != PATTERN ||
	  !std::getline(std::cin, line2) || splitString(line2, line2_col) != 2 ||
	  !std::getline(std::cin, line3) || splitString(line3, line3_col) != 2)
         throw std::runtime_error("invalid input");

      std::string lines = line1 + EOL + line2 + EOL + line3 + EOL;

      hitVector.push_back(new Hit(line1a_col[1], line2_col[0], lines));
   }

   return numReadPairs;
}

//------------------------------------------------------------------------------------
// writeHits() writes the hits in the vector to stdout, followed by the number of
// read pairs

void writeHits(const HitVector& hitVector, uint64_t numReadPairs)
{
   int numHits = hitVector.size();

   for (int i = 0; i < numHits; i++)
      std::cout << hitVector[i]->lines;

   std::cout << numReadPairs << TAB << READ_PAIRS << NEWLINE;
}

//------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   if (argc > 1)
   {
      showUsage(argv[0]);
      return 1;
   }

   try
   {
      HitVector hitVector;

      uint64_t numReadPairs = readHits(hitVector);

      sortHits(hitVector);

      writeHits(hitVector, numReadPairs);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
