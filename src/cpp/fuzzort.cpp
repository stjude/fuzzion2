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

const std::string VERSION_ID = "fuzzort v1.0.2, copyright 2021 "
                               "St. Jude Children's Research Hospital";

const std::string FUZZION2   = "fuzzion2 ";
const std::string PATTERN    = "pattern ";
const std::string READ       = "read ";
const std::string READ_PAIRS = "read-pairs ";

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_ID << NEWLINE << NEWLINE
      << "Usage: " << progname << " < fuzzion2_hits > sorted_hits" << NEWLINE;
}

//------------------------------------------------------------------------------------

class Hit
{
public:
   Hit(const std::string& inPatternLine, const std::string& inRead1Line,
       const std::string& inRead2Line)
      : patternLine(inPatternLine), read1Line(inRead1Line), read2Line(inRead2Line) { }

   virtual ~Hit() { }

   std::string patternLine, read1Line, read2Line;
};

typedef std::vector<Hit *> HitVector;

struct HitCompare
{
   bool operator()(Hit* const& a, Hit* const& b) const
   {
      // order by ascending patternLine, then by ascending read1Line
      return (a->patternLine == b->patternLine ?
              (a->read1Line   < b->read1Line) :
	      (a->patternLine < b->patternLine));
   }
};

void sortHits(HitVector& hitVector)
{
   std::sort(hitVector.begin(), hitVector.end(), HitCompare());
}

//------------------------------------------------------------------------------------
// readHits() reads the hits from stdin and stores them in the vector; the number of
// read pairs is summed and returned to the caller; the heading line is also provided

uint64_t readHits(HitVector& hitVector, std::string& headingLine)
{
   if (!std::getline(std::cin, headingLine))
      throw std::runtime_error("no input");

   if (!hasPrefix(headingLine, FUZZION2))
      throw std::runtime_error("missing or invalid heading line");

   uint64_t numReadPairs = 0;
   std::string line1, line2, line3;

   while (std::getline(std::cin, line1))
      if (hasPrefix(line1, FUZZION2)) // found another heading line
         if (line1 == headingLine)
            continue; // ignore it
         else
            throw std::runtime_error("inconsistent heading lines");
      else if (hasPrefix(line1, PATTERN)) // found a hit
         if (std::getline(std::cin, line2) && hasPrefix(line2, READ) &&
             std::getline(std::cin, line3) && hasPrefix(line3, READ))
            hitVector.push_back(new Hit(line1, line2, line3));
         else
            throw std::runtime_error("invalid input after " + line1);
      else if (hasPrefix(line1, READ_PAIRS)) // found number of read pairs
      {
         StringVector part;

	 if (splitString(line1, part, ' ') == 2)
	 {
            std::istringstream stream(part[1]);
	    uint64_t count = 0;
	    stream >> count;
	    numReadPairs += count;
	 }
      }
      else
         throw std::runtime_error("invalid input " + line1);

   return numReadPairs;
}

//------------------------------------------------------------------------------------
// writeHits() writes the hits in the vector to stdout, followed by the number of
// read pairs

void writeHits(const HitVector& hitVector, const std::string& headingLine,
               uint64_t numReadPairs)
{
   std::cout << headingLine << NEWLINE;

   int numHits = hitVector.size();

   for (int i = 0; i < numHits; i++)
      std::cout << hitVector[i]->patternLine << NEWLINE
                << hitVector[i]->read1Line   << NEWLINE
		<< hitVector[i]->read2Line   << NEWLINE;

   std::cout << READ_PAIRS << numReadPairs << NEWLINE;
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
      HitVector   hitVector;
      std::string headingLine;

      uint64_t numReadPairs = readHits(hitVector, headingLine);

      sortHits(hitVector);

      writeHits(hitVector, headingLine, numReadPairs);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
