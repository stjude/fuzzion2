//------------------------------------------------------------------------------------
//
// summary.cpp - module for reading and writing fuzzion2 hit summaries
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "summary.h"

// column headings
const String TOTAL     = "total";
const String DISTINCT  = "distinct";
const String WEAK      = HIT_WEAK;
const String STNOSPAN  = HIT_STRONG_NOSPAN;
const String STSPAN    = HIT_STRONG_SPAN;
const String PATTERN   = "pattern";
const String GROUP     = "pattern group";

const int SAMPLE_COL   = 0;
const int TOTAL_COL    = 1;
const int DISTINCT_COL = 2;
const int WEAK_COL     = 3;
const int STNOSPAN_COL = 4;
const int STSPAN_COL   = 5;
const int NAME_COL     = 6;

//------------------------------------------------------------------------------------
// summarizeHits() returns a newly-allocated Summary of the hits in the given vector
// from indices begin (inclusive) to end (exclusive)

Summary *summarizeHits(const HitVector& hitVector, const int begin, const int end,
                       const int minStrong, const String& sampleID)
{
   const int numMatches = end - begin;
   int weak = 0, strongNospan = 0, strongSpan = 0;

   for (int i = begin; i < end; i++)
   {
      const String label = hitVector[i]->label(minStrong);

      if (label == HIT_WEAK)
         weak++;
      else if (label == HIT_STRONG_NOSPAN)
         strongNospan++;
      else if (label == HIT_STRONG_SPAN)
         strongSpan++;
   }

   return new Summary(sampleID, numMatches, weak, strongNospan, strongSpan,
                      hitVector[begin]->patternName, hitVector[begin]->annotation);
}

//------------------------------------------------------------------------------------
// SummaryCompare defines the sort order of hit summaries

struct SummaryCompare
{
   bool operator()(Summary* const& a, Summary* const& b) const
   {
      // sort by ascending name of pattern or group, then by ascending sample ID

      const int key1 = std::strcmp(a->name.c_str(), b->name.c_str());
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
// writeSummaryHeadingLine() writes a heading line

void writeSummaryHeadingLine(std::ostream& ostream, const String& version,
                             const bool grouping,
                             const StringVector& annotationHeading)
{
   ostream << FUZZUM << version
           << TAB << TOTAL
           << TAB << DISTINCT
           << TAB << WEAK
           << TAB << STNOSPAN
           << TAB << STSPAN
           << TAB << (grouping ? GROUP : PATTERN);

   const int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      ostream << TAB << annotationHeading[i];

   ostream << NEWLINE;
}

//------------------------------------------------------------------------------------
// isHeadingLine() returns true if the given line is a heading line

static bool isHeadingLine(const String& line)
{
   return hasPrefix(line, FUZZUM);
}

//------------------------------------------------------------------------------------
// validHeadingLine() returns true if the given line is a valid heading line

static bool validHeadingLine(const String& line, StringVector& annotationHeading)
{
   if (!isHeadingLine(line))
      return false;

   StringVector col;
   const int numCols = splitString(line, col);

   if (numCols <= NAME_COL || col[TOTAL_COL] != TOTAL ||
       col[DISTINCT_COL] != DISTINCT || col[WEAK_COL]   != WEAK     ||
       col[STNOSPAN_COL] != STNOSPAN || col[STSPAN_COL] != STSPAN   ||
       col[NAME_COL]     != PATTERN  && col[NAME_COL]   != GROUP)
      return false;

   annotationHeading.clear();

   for (int i = NAME_COL + 1; i < numCols; i++)
      annotationHeading.push_back(col[i]);

   return true;
}

//------------------------------------------------------------------------------------
// Summary::write() writes one line summarizing the matches to a pattern or group

void Summary::write(std::ostream& ostream) const
{
   ostream << sampleID
           << TAB << numMatches
           << TAB << distinct()
           << TAB << weak
           << TAB << strongNospan
           << TAB << strongSpan
           << TAB << name;

   const int numAnnotations = annotation.size();

   for (int i = 0; i < numAnnotations; i++)
      ostream << TAB << annotation[i];

   ostream << NEWLINE;
}

//------------------------------------------------------------------------------------
// getSummary() returns a pointer to a newly-allocated Summary obtained from the given
// line, or returns nullptr if the input is invalid

static Summary *getSummary(const String& line)
{
   StringVector col;
   const int numCols = splitString(line, col);

   int numMatches, weak, strongNospan, strongSpan;

   if (numCols <= NAME_COL ||
       (numMatches   = stringToInt(col[TOTAL_COL]))    <= 0 ||
       (weak         = stringToInt(col[WEAK_COL]))     <  0 ||
       (strongNospan = stringToInt(col[STNOSPAN_COL])) <  0 ||
       (strongSpan   = stringToInt(col[STSPAN_COL]))   <  0)
      return nullptr;

   const String& sampleID = col[SAMPLE_COL];
   const String& name     = col[NAME_COL];

   StringVector annotation;

   for (int i = NAME_COL + 1; i < numCols; i++)
      annotation.push_back(col[i]);

   return new Summary(sampleID, numMatches, weak, strongNospan, strongSpan, name,
                      annotation);
}

//------------------------------------------------------------------------------------
// readSummaries() reads hit summaries from one or more input files and stores them in
// sorted order in summaryVector (it is the caller's obligation to de-allocate them)

void readSummaries(const StringVector& filename, StringVector& annotationHeading,
                   SummaryVector& summaryVector)
{
   String headingLine, line;
   const int numFiles = filename.size();

   for (int i = 0; i < numFiles; i++)
   {
      std::ifstream infile(filename[i].c_str());
      if (!infile.is_open())
         throw Error("unable to open " + filename[i]);

      if (!getline(infile, line))
         throw Error("empty file " + filename[i]);

      if (!validHeadingLine(line, annotationHeading))
         throw Error("unexpected heading line in " + filename[i]);

      if (i == 0)
         headingLine = line;
      else if (line != headingLine)
         throw Error("inconsistent heading lines in input files");

      while (getline(infile, line))
         if (isHeadingLine(line)) // found another heading line in this file
         {
            if (line != headingLine)
               throw Error("inconsistent heading lines in " + filename[i]);
         }
         else
         {
            Summary *summary = getSummary(line);

            if (summary)
               summaryVector.push_back(summary);
            else
               throw Error("unexpected summary format in " + filename[i] + ": " +
                           line);
         }

      infile.close();
   }

   if (summaryVector.size() > 1)
      sortSummaries(summaryVector);
}
