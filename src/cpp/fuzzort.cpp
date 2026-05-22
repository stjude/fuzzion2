//------------------------------------------------------------------------------------
//
// fuzzort.cpp - this program reads fuzzion2 hits from stdin, sorts them, and writes
//               them to stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "hit.h"
#include "version.h"

const String VERSION_NAME = "fuzzort " + CURRENT_VERSION;

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

int main(const int argc, const char *argv[])
{
   if (argc > 1)
   {
      showUsage(argv[0]);
      return 1;
   }

   try
   {
      String       fuzzion2Version;
      StringVector annotationHeading;
      HitVector    hitVector;
      uint64_t     numReads;

      readHits(std::cin, fuzzion2Version, annotationHeading, hitVector, numReads);

      writeHitHeadingLine(std::cout, fuzzion2Version, annotationHeading);

      const int numHits = hitVector.size();

      for (int i = 0; i < numHits; i++)
         hitVector[i]->write(std::cout);

      writeReadCountLine(std::cout, numReads);
   }
   catch (const Error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
