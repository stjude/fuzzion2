//------------------------------------------------------------------------------------
//
// fuzzum.cpp - this programs reads fuzzion2 hits from stdin and writes a summary to
//              stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "hit.h"
#include "summary.h"
#include "version.h"
#include <iostream>
#include <stdexcept>

const std::string VERSION_NAME = FUZZUM + CURRENT_VERSION;

int minStrong = DEFAULT_MIN_STRONG; // minimum overlap for a strong match

std::string id = "";                // identifies the sample

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname << " OPTION ... < fuzzion2_hits > hit_summary"
      << NEWLINE;

   std::cerr
      << NEWLINE
      << "This option is required:" << NEWLINE
      << "  -id=string   "
             << "identifies the sample" << NEWLINE;

   std::cerr
      << NEWLINE
      << "The following is optional:" << NEWLINE
      << "  -strong=N    "
             << "minimum overlap of a strong match in #bases, default is "
	     << DEFAULT_MIN_STRONG << NEWLINE;
}

//------------------------------------------------------------------------------------
// parseArgs() parses the command-line arguments and returns true if all are valid

bool parseArgs(int argc, char *argv[])
{
   for (int i = 1; i < argc; i++)
   {
      std::string arg = argv[i];
      if (arg.length() == 0)
         continue;

      if (arg[0] != '-')
         return false; // not an option

      StringVector opt;

      if (splitString(arg, opt, '=') != 2)
         return false; // incorrect option format

      if (stringOpt(opt, "id",     id) ||
          intOpt   (opt, "strong", minStrong))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   return (id != "" && minStrong > 0);
}

//------------------------------------------------------------------------------------
// summarizePattern() writes to stdout a summary of the hits in hitVector from begin
// (inclusive) to end (exclusive); it is assumed that hitVector has been sorted

void summarizePattern(const HitVector& hitVector, int begin, int end)
{
   std::vector<int> index;
   getDistinctIndices(hitVector, begin, end, index);

   int numDistinct = index.size();
   int numStrong   = 0;

   for (int i = 0; i < numDistinct; i++)
      if (hitVector[index[i]]->isStrong(minStrong))
         numStrong++;

   Summary summary(hitVector[begin]->pattern->name, id, end - begin, numDistinct,
                   numStrong, hitVector[begin]->pattern->annotation);

   summary.write();
}

//------------------------------------------------------------------------------------
// writeSummaries() writes summaries to stdout for the hits in hitVector

void writeSummaries(const StringVector& annotationHeading,
                    const HitVector& hitVector)
{
   writeSummaryHeadingLine(CURRENT_VERSION, annotationHeading);

   std::vector<int> index;
   getPatternIndices(hitVector, index);

   int numPatterns = index.size();

   for (int i = 0; i < numPatterns; i++)
      summarizePattern(hitVector, index[i],
                       (i + 1 < numPatterns ? index[i + 1] : hitVector.size()));
}

//------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   if (!parseArgs(argc, argv))
   {
      showUsage(argv[0]);
      return 1;
   }

   try
   {
      std::string  fuzzion2Version;
      StringVector annotationHeading;
      HitVector    hitVector;

      readHits(fuzzion2Version, annotationHeading, hitVector);

      writeSummaries(annotationHeading, hitVector);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
