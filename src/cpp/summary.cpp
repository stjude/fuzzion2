//------------------------------------------------------------------------------------
//
// summary.cpp - module for reading and writing fuzzion2 hit summaries
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "summary.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

// column headings
const std::string READPAIR = "read pairs";
const std::string DISTINCT = "distinct";
const std::string WEAK     = "weak";
const std::string STNOSPAN = "strong-";
const std::string STSPAN   = "strong+";
const std::string PATTERN  = "pattern";
const std::string GROUP    = "pattern group";

const int SAMPLE_COL       = 0;
const int READPAIR_COL     = 1;
const int DISTINCT_COL     = 2;
const int WEAK_COL         = 3;
const int STNOSPAN_COL     = 4;
const int STSPAN_COL       = 5;
const int NAME_COL         = 6;

//------------------------------------------------------------------------------------
// Summary::write() writes one line to stdout summarizing the read pairs of a sample
// that match a pattern or group

void Summary::write() const
{
   std::cout << sampleID
             << TAB << readPairs
	     << TAB << distinct()
	     << TAB << weak
	     << TAB << strongNospan
	     << TAB << strongSpan
	     << TAB << name;

   int numAnnotations = annotation.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotation[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeSummaryHeadingLine() writes a heading line to stdout

void writeSummaryHeadingLine(const std::string& version, bool grouping,
                             const StringVector& annotationHeading)
{
   std::cout << FUZZUM << version
             << TAB << READPAIR
	     << TAB << DISTINCT
	     << TAB << WEAK
	     << TAB << STNOSPAN
	     << TAB << STSPAN
	     << TAB << (grouping ? GROUP : PATTERN);

   int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotationHeading[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// summarizeHits() returns a newly-allocated Summary of the hits in the given vector
// from indices begin (inclusive) to end (exclusive)

Summary *summarizeHits(const HitVector& hitVector, int begin, int end, int minStrong,
                       std::string sampleID)
{
   int readPairs = end - begin;
   int weak = 0, strongNospan = 0, strongSpan = 0;

   for (int i = begin; i < end; i++)
   {
      std::string label = hitVector[i]->label(minStrong);

      if (label == HIT_WEAK)
         weak++;
      else if (label == HIT_STRONG_NOSPAN)
         strongNospan++;
      else if (label == HIT_STRONG_SPAN)
         strongSpan++;
   }

   return new Summary(sampleID, readPairs, weak, strongNospan, strongSpan,
                      hitVector[begin]->pattern->name,
		      hitVector[begin]->pattern->annotation);
}

//------------------------------------------------------------------------------------
// isHeadingLine() returns true if the given line is a heading line

static bool isHeadingLine(const std::string& line)
{
   return hasPrefix(line, FUZZUM);
}

//------------------------------------------------------------------------------------
// isValidHeadingLine() returns true if the given line is a valid heading line

static bool isValidHeadingLine(const std::string& line,
                               StringVector& annotationHeading)
{
   StringVector col;
   int numCols = splitString(line, col);

   if (!isHeadingLine(line) || numCols <= NAME_COL ||
       col[READPAIR_COL] != READPAIR || col[DISTINCT_COL] != DISTINCT ||
       col[WEAK_COL]     != WEAK     || col[STNOSPAN_COL] != STNOSPAN ||
       col[STSPAN_COL]   != STSPAN   ||
       col[NAME_COL] != PATTERN && col[NAME_COL] != GROUP)
      return false;

   annotationHeading.clear();

   for (int i = NAME_COL + 1; i < numCols; i++)
      annotationHeading.push_back(col[i]);

   return true;
}

//------------------------------------------------------------------------------------
// getSummary() returns a pointer to a newly-allocated Summary obtained from the given
// line, or returns NULL if the input is invalid

static Summary *getSummary(const std::string& line)
{
   StringVector col;
   int numCols = splitString(line, col);

   int readPairs, weak, strongNospan, strongSpan;

   if (numCols <= NAME_COL ||
       (readPairs    = stringToNonnegInt(col[READPAIR_COL])) <= 0 ||
       (weak         = stringToNonnegInt(col[WEAK_COL]))     <  0 ||
       (strongNospan = stringToNonnegInt(col[STNOSPAN_COL])) <  0 ||
       (strongSpan   = stringToNonnegInt(col[STSPAN_COL]))   <  0)
      return NULL;

   const std::string& sampleID = col[SAMPLE_COL];
   const std::string& name     = col[NAME_COL];

   StringVector annotation;

   for (int i = NAME_COL + 1; i < numCols; i++)
      annotation.push_back(col[i]);

   return new Summary(sampleID, readPairs, weak, strongNospan, strongSpan, name,
                      annotation);
}

//------------------------------------------------------------------------------------
// SummaryCompare defines the sort order of hit summaries

struct SummaryCompare
{
   bool operator()(Summary* const& a, Summary* const& b) const
   {
      // sort by ascending name of pattern or group, then by ascending sample ID

      int key1 = std::strcmp(a->name.c_str(), b->name.c_str());
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

      if (!isValidHeadingLine(line, annotationHeading))
         throw std::runtime_error("unexpected heading line in " + filename[i]);

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
