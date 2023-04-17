//------------------------------------------------------------------------------------
//
// fuzzall.cpp - this program reads summaries produced by fuzzum and writes an
//               aggregate summary to stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "summary.h"
#include "version.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>

const std::string VERSION_NAME = "fuzzall " + CURRENT_VERSION;

std::string  datasetName = "";
StringVector fuzzumFilename;

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname << " OPTION fuzzum_filename ... > pattern_summary"
      << NEWLINE;

   std::cerr
      << NEWLINE
      << "The following is optional:" << NEWLINE
      << "  -dataset=name   name associated with this dataset" << NEWLINE;
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
      {
         fuzzumFilename.push_back(arg);
	 continue;
      }

      StringVector opt;

      if (splitString(arg, opt, '=') != 2)
         return false; // incorrect option format

      if (stringOpt(opt, "dataset", datasetName))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   return (fuzzumFilename.size() > 0);
}

//------------------------------------------------------------------------------------

class Stats // for collecting statistics
{
public:
   Stats()
      : min(-1), max(-1), sum(0), value() { }

   virtual ~Stats() { }

   void addValue(int newValue);

   double median();
   double mean() const;

   void write();

   int min, max, sum;
   IntVector value;
};

//------------------------------------------------------------------------------------
// Stats::addValue() incorporates the given value into the statistics

void Stats::addValue(int newValue)
{
   if (value.size() == 0)
      min = max = newValue;
   else if (newValue < min)
      min = newValue;
   else if (newValue > max)
      max = newValue;

   sum += newValue;

   value.push_back(newValue);
}

//------------------------------------------------------------------------------------
// Stats::median() returns the median of the set of values

double Stats::median()
{
   int n = value.size();

   if (n == 0)
      return -1.0; // median is undefined for the empty set

   if (n > 1)
      std::sort(value.begin(), value.end());

   if (n % 2 == 0) // even number of values
      return ((value[n / 2 - 1] + value[n / 2]) / 2.0);
   else // odd
      return value[n / 2];
}

//------------------------------------------------------------------------------------
// Stats::mean() returns the mean of the set of values

double Stats::mean() const
{
   int n = value.size();

   if (n == 0)
      return -1.0; // mean is undefined for the empty set

   return (sum / static_cast<double>(n));
}

//------------------------------------------------------------------------------------
// Stats::write() writes to stdout the statistics for a set of values

void Stats::write()
{
   std::cout << TAB << sum
             << TAB << min
	     << TAB << doubleToString(median())
	     << TAB << doubleToString(mean())
	     << TAB << max;
}

//------------------------------------------------------------------------------------
// writeStatsHeadings() writes to stdout the headings for one set of statistics

void writeStatsHeadings(std::string category)
{
   std::cout << TAB << category
             << TAB << "min"
	     << TAB << "median"
	     << TAB << "mean"
	     << TAB << "max";
}

//------------------------------------------------------------------------------------
// writeHeadingLine() writes a heading line to stdout

void writeHeadingLine(const StringVector& annotationHeading)
{
   std::cout << VERSION_NAME;

   if (datasetName != "")
      std::cout << TAB << "dataset";

   writeStatsHeadings("distinct");
   writeStatsHeadings("weak");
   writeStatsHeadings("strong-");
   writeStatsHeadings("strong+");

   std::cout << TAB << "IDs"
             << TAB << "ID list";

   int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotationHeading[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeDataLine() writes one line to stdout with the aggregate summary for one
// pattern or group

void writeDataLine(const std::string& name, Stats& distinctStats, Stats& weakStats,
                   Stats& strongNospanStats, Stats& strongSpanStats, int numIDs,
		   const std::string& idList, const StringVector& annotation)
{
   std::cout << name;

   if (datasetName != "")
      std::cout << TAB << datasetName;

   distinctStats.write();
   weakStats.write();
   strongNospanStats.write();
   strongSpanStats.write();

   std::cout << TAB << numIDs
             << TAB << idList;

   int numAnnotations = annotation.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotation[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// aggregateOne() aggregates the data for one pattern or group and writes it to
// stdout; the first summary is in summaryVector[start]; the return value is the index
// of the first summary of the next pattern or group, or summaryVector.size() if none

int aggregateOne(const SummaryVector& summaryVector, int start)
{
   const std::string& name = summaryVector[start]->name;
   Stats distinctStats, weakStats, strongNospanStats, strongSpanStats;
   int numIDs = 0;
   std::string idList = "";
   const StringVector& annotation = summaryVector[start]->annotation;

   int numSummaries = summaryVector.size();
   int i = start;

   do
   {
      const std::string& sampleID = summaryVector[i]->sampleID;
      int numDistinct = 0, numWeak = 0, numStrongNospan = 0, numStrongSpan = 0;

      do
      {
         numDistinct     += summaryVector[i]->distinct();
	 numWeak         += summaryVector[i]->weak;
	 numStrongNospan += summaryVector[i]->strongNospan;
	 numStrongSpan   += summaryVector[i]->strongSpan;
      }
      while (++i < numSummaries && summaryVector[i]->name == name &&
	     summaryVector[i]->sampleID == sampleID);

      distinctStats.addValue(numDistinct);
      weakStats.addValue(numWeak);
      strongNospanStats.addValue(numStrongNospan);
      strongSpanStats.addValue(numStrongSpan);

      if (++numIDs > 1)
         idList += ", ";

      idList += sampleID + "(" + intToString(numDistinct) + "/" +
                intToString(numStrongSpan) + ")";
   }
   while (i < numSummaries && summaryVector[i]->name == name);

   writeDataLine(name, distinctStats, weakStats, strongNospanStats, strongSpanStats,
                 numIDs, idList, annotation);

   return i;
}

//------------------------------------------------------------------------------------
// aggregateAll() aggregates the data for all patterns or groups and writes it to
// stdout

void aggregateAll(const SummaryVector& summaryVector)
{
   int numSummaries = summaryVector.size();
   int start = 0;

   while (start < numSummaries)
      start = aggregateOne(summaryVector, start);
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
      StringVector  annotationHeading;
      SummaryVector summaryVector;

      readSummaries(fuzzumFilename, annotationHeading, summaryVector);

      writeHeadingLine(annotationHeading);
      aggregateAll(summaryVector);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
