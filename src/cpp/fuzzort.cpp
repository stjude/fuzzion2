//------------------------------------------------------------------------------------
//
// fuzzort.cpp - this program reads fuzzion2 hits from stdin, sorts them, and writes
//               them to stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "hit.h"
#include "version.h"
#include <iostream>
#include <stdexcept>

const std::string VERSION_NAME = "fuzzort " + CURRENT_VERSION;

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname << " < fuzzion2_hits > sorted_hits"
      << NEWLINE;
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
      std::string fuzzion2Version;
      StringVector annotationHeading;
      HitVector hitVector;

      uint64_t numReadPairs = readHits(std::cin, fuzzion2Version, annotationHeading,
                                       hitVector);

      writeHitHeadingLine(fuzzion2Version, annotationHeading);

      int numHits = hitVector.size();

      for (int i = 0; i < numHits; i++)
         hitVector[i]->write();

      writeReadPairLine(numReadPairs);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
