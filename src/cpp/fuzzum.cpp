//------------------------------------------------------------------------------------
//
// fuzzum.cpp - this programs reads fuzzion2 hits from stdin and writes a summary to
//              stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "group.h"
#include "version.h"

const String VERSION_NAME = FUZZUM + CURRENT_VERSION;

int minStrong = DEFAULT_MIN_STRONG; // minimum overlap for a strong match

String id = "";           // identifies the sample
String groupColList = ""; // comma-separated list of group column headings

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
      << "  -id=string      "
             << "identifies the sample" << NEWLINE;

   std::cerr
      << NEWLINE
      << "The following are optional:" << NEWLINE
      << "  -group=string   "
             << "comma-separated list of column headings, default is no grouping"
	     << NEWLINE
      << "  -strong=N       "
             << "minimum overlap of a strong match in #bases, default is "
	     << DEFAULT_MIN_STRONG << NEWLINE;
}

//------------------------------------------------------------------------------------
// parseArgs() parses the command-line arguments and returns true if all are valid

bool parseArgs(const int argc, const char *argv[])
{
   for (int i = 1; i < argc; i++)
   {
      const String arg = argv[i];
      if (arg.length() == 0)
         continue;

      if (arg[0] != '-')
         return false; // not an option

      StringVector opt;

      if (splitString(arg, opt, '=') != 2)
         return false; // incorrect option format

      if (stringOpt(opt, "id",     id)           ||
          stringOpt(opt, "group",  groupColList) ||
          intOpt   (opt, "strong", minStrong))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   return (id != "" && minStrong >= TRIMER_LEN);
}

//------------------------------------------------------------------------------------
// writePatternSummaries() writes to stdout a summary of hits for each pattern; it is
// assumed that hitVector has been sorted

void writePatternSummaries(const StringVector& annotationHeading,
                           const HitVector& hitVector)
{
   writeSummaryHeadingLine(std::cout, CURRENT_VERSION, false, annotationHeading);

   IntVector index;
   getPatternIndices(hitVector, index);
   const int numPatterns = index.size();

   for (int i = 0; i < numPatterns; i++)
   {
      const int begin = index[i];
      const int end   = (i + 1 < numPatterns ? index[i + 1] : hitVector.size());

      const Summary *summary = summarizeHits(hitVector, begin, end, minStrong, id);
      summary->write(std::cout);
      delete summary;
   }
}

//------------------------------------------------------------------------------------
// writeGroupSummaries() writes to stdout a summary of hits for each pattern group

void writeGroupSummaries(const GroupManager& groupManager)
{
   writeSummaryHeadingLine(std::cout, CURRENT_VERSION, true,
                           groupManager.annotationHeading);

   const GroupMap& gmap = groupManager.gmap;

   for (GroupMap::const_iterator gpos = gmap.cbegin(); gpos != gmap.cend(); ++gpos)
   {
      const Group& group = gpos->second;

      const Summary *summary = group.summarize(minStrong, id);
      summary->write(std::cout);
      delete summary;
   }
}

//------------------------------------------------------------------------------------

int main(const int argc, const char *argv[])
{
   if (!parseArgs(argc, argv))
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

      if (groupColList == "")
         writePatternSummaries(annotationHeading, hitVector);
      else
      {
         const GroupManager groupManager(groupColList, annotationHeading, hitVector);
	 writeGroupSummaries(groupManager);
      }
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
