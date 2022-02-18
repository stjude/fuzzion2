//------------------------------------------------------------------------------------
//
// summary.cpp - module for reading and writing fuzzion2 hit summaries
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "summary.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

// column headings
const std::string READ_PAIRS = "read pairs";
const std::string DISTINCT   = "distinct";
const std::string STRONG     = "strong";
const std::string PATTERN    = "pattern";

const int SAMPLE_COL         = 0;
const int RPAIRS_COL         = 1;
const int DISTINCT_COL       = 2;
const int STRONG_COL         = 3;
const int PATTERN_COL        = 4; const int OLD_PATTERN_COL = 2;
const int MIN_HEADING_COLS   = 5;

//------------------------------------------------------------------------------------
// Summary::write() writes one line to stdout summarizing the read pairs of a sample
// that match a pattern

void Summary::write() const
{
   std::cout << sampleID
             << TAB << readPairs
	     << TAB << distinct
	     << TAB << strong
	     << TAB << patternName;

   int numAnnotations = annotation.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotation[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeSummaryHeadingLine() writes a heading line to stdout

void writeSummaryHeadingLine(const std::string& version,
                             const StringVector& annotationHeading)
{
   std::cout << FUZZUM << version
             << TAB << READ_PAIRS
	     << TAB << DISTINCT
	     << TAB << STRONG
	     << TAB << PATTERN;

   int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotationHeading[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// isHeadingLine() returns true if the given line is a heading line

static bool isHeadingLine(const std::string& line)
{
   return hasPrefix(line, FUZZUM);
}

//------------------------------------------------------------------------------------
// validateHeadingLine() verifies that the given line is a valid heading line

static void validateHeadingLine(const std::string& line, const std::string& filename,
                               StringVector& annotationHeading)
{
   StringVector col;
   int numCols = splitString(line, col);

   if (!isHeadingLine(line) || numCols < MIN_HEADING_COLS ||
       col[RPAIRS_COL] != READ_PAIRS || col[DISTINCT_COL] != DISTINCT ||
       col[STRONG_COL] != STRONG     || col[PATTERN_COL]  != PATTERN)
      if (isHeadingLine(line) && numCols > OLD_PATTERN_COL &&
          col[OLD_PATTERN_COL] == PATTERN)
         throw std::runtime_error(FUZZUM + "format is old in " + filename +
                                  "; rerun " + FUZZUM + "using newer version");
      else
         throw std::runtime_error("unexpected heading line in " + filename);

   annotationHeading.clear();

   for (int i = MIN_HEADING_COLS; i < numCols; i++)
      annotationHeading.push_back(col[i]);
}

//------------------------------------------------------------------------------------
// getSummary() returns a pointer to a newly-allocated Summary obtained from the given
// line, or returns NULL if the input is invalid

static Summary *getSummary(const std::string& line)
{
   int readPairs, distinct, strong;
   StringVector annotation, col;

   int numCols = splitString(line, col);

   if (numCols < MIN_HEADING_COLS ||
       (readPairs = stringToNonnegInt(col[RPAIRS_COL]))   <= 0 ||
       (distinct  = stringToNonnegInt(col[DISTINCT_COL])) <= 0 ||
       (strong    = stringToNonnegInt(col[STRONG_COL]))   <  0)
      return NULL;

   for (int i = MIN_HEADING_COLS; i < numCols; i++)
      annotation.push_back(col[i]);

   return new Summary(col[PATTERN_COL], col[SAMPLE_COL], readPairs, distinct, strong,
                      annotation);
}

//------------------------------------------------------------------------------------
// SummaryCompare defines the sort order of hit summaries

struct SummaryCompare
{
   bool operator()(Summary* const& a, Summary* const& b) const
   {
      // sort by ascending pattern name, then by ascending sample ID

      int key1 = std::strcmp(a->patternName.c_str(), b->patternName.c_str());
      if (key1 != 0)
         return (key1 < 0);

      return (a->sampleID < b->sampleID);
   }
};

//------------------------------------------------------------------------------------
// sortSummaries() sorts the summaries in the given vector

static void sortSummaries(SummaryVector& summaryVector)
{
   std::sort(summaryVector.begin(), summaryVector.end(), SummaryCompare());
}

//------------------------------------------------------------------------------------
// readSummaries() reads hit summaries from one or more input files and stores them in
// sorted order in summaryVector

void readSummaries(const StringVector& filename, StringVector& annotationHeading,
                   SummaryVector& summaryVector)
{
   std::string headingLine, line;
   int numFiles = filename.size();

   for (int i = 0; i < numFiles; i++)
   {
      std::ifstream infile(filename[i].c_str());
      if (!infile.is_open())
         throw std::runtime_error("unable to open " + filename[i]);

      if (!getline(infile, line))
         throw std::runtime_error("empty file " + filename[i]);

      validateHeadingLine(line, filename[i], annotationHeading);

      if (i == 0)
         headingLine = line;
      else if (line != headingLine)
         throw std::runtime_error("inconsistent heading lines in input files");

      while (getline(infile, line))
         if (isHeadingLine(line)) // found another heading line in this file
	 {
            if (line != headingLine)
               throw std::runtime_error("inconsistent heading lines in " +
                                        filename[i]);
	 }
         else
	 {
            Summary *summary = getSummary(line);

	    if (summary)
               summaryVector.push_back(summary);
	    else
               throw std::runtime_error("unexpected summary format in " +
                                        filename[i] + ": " + line);
	 }

      infile.close();
   }

   if (summaryVector.size() > 1)
      sortSummaries(summaryVector);
}
