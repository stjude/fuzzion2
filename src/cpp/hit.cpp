//------------------------------------------------------------------------------------
//
// hit.cpp - module for reading and writing fuzzion2 hits
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "hit.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

const std::string PATTERN    = "pattern ";
const std::string READ       = "read ";
const std::string READ_PAIRS = "read-pairs ";

// column headings
const std::string SEQUENCE   = "sequence";
const std::string MBASES     = "matching bases";
const std::string POSSIBLE   = "possible";
const std::string PERCENT    = "% match";
const std::string SPANNING   = "junction spanning";
const std::string OVLEFT     = "left overlap";
const std::string OVRIGHT    = "right overlap";
const std::string ISIZE      = "insert size";

const int SEQUENCE_COL       = 1;
const int MBASES_COL         = 2;
const int POSSIBLE_COL       = 3;
const int PERCENT_COL        = 4;
const int SPANNING_COL       = 5;
const int OVLEFT_COL         = 6;
const int OVRIGHT_COL        = 7;
const int ISIZE_COL          = 8;
const int FIRST_ANNOT_COL    = 9; // first annotation column, if any

//------------------------------------------------------------------------------------
// HitPattern::write() writes one line to stdout describing the pattern matched by a
// read pair

void HitPattern::write() const
{
   std::cout << PATTERN << name
             << TAB << displaySequence
	     << TAB << matchingBases
	     << TAB << possible
	     << TAB << doubleToString(percentMatch())
	     << TAB << spanningCount
	     << TAB // no entry in this column
	     << TAB // no entry in this column
	     << TAB << insertSize;

   int numAnnotations = annotation.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotation[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// HitRead::write() writes one line to stdout describing one read of a read pair that
// matches a pattern

void HitRead::write() const
{
   std::cout << READ << name
             << TAB << std::string(leadingBlanks, ' ') << sequence
	     << TAB << matchingBases
	     << TAB << possible()
	     << TAB << doubleToString(percentMatch())
	     << TAB << (isSpanning ? 1 : 0)
	     << TAB << leftOverlap
	     << TAB << rightOverlap
	     << NEWLINE;
}

//------------------------------------------------------------------------------------
// Hit::sameAs() returns true if this hit and the other hit are duplicates

bool Hit::sameAs(const Hit& other) const
{
   return (pattern->name == other.pattern->name &&
           pattern->leftBases  == other.pattern->leftBases &&
	   pattern->rightBases == other.pattern->rightBases);
}

//------------------------------------------------------------------------------------
// Hit::isStrong() returns true if this hit is considered a "strong" hit

bool Hit::isStrong(int minStrong) const
{
   return (std::max(read1->leftOverlap,  read2->leftOverlap)  >= minStrong &&
           std::max(read1->rightOverlap, read2->rightOverlap) >= minStrong);
}

//------------------------------------------------------------------------------------
// Hit::label() returns the hit category

std::string Hit::label(int minStrong) const
{
   if (duplicate)
      return HIT_DUPLICATE;

   if (isStrong(minStrong))
      if (isSpanning())
         return HIT_STRONG_SPAN;
      else
         return HIT_STRONG_NOSPAN;

   return HIT_WEAK;
}

//------------------------------------------------------------------------------------
// writeHitHeadingLine() writes a heading line to stdout

void writeHitHeadingLine(const std::string& version,
                         const StringVector& annotationHeading)
{
   std::cout << FUZZION2 << version
             << TAB << SEQUENCE
	     << TAB << MBASES
	     << TAB << POSSIBLE
	     << TAB << PERCENT
	     << TAB << SPANNING
	     << TAB << OVLEFT
	     << TAB << OVRIGHT
	     << TAB << ISIZE;

   int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotationHeading[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// isHeadingLine() returns true if the given line is a heading line

static bool isHeadingLine(const std::string& line)
{
   return hasPrefix(line, FUZZION2);
}

//------------------------------------------------------------------------------------
// isValidHeadingLine() returns true if the given line is a valid heading line

static bool isValidHeadingLine(const std::string& line, std::string& version,
                               StringVector& annotationHeading)
{
   StringVector col;
   int numCols = splitString(line, col);

   if (!isHeadingLine(line) || numCols <= ISIZE_COL ||
       col[SEQUENCE_COL] != SEQUENCE || col[MBASES_COL]  != MBASES  ||
       col[POSSIBLE_COL] != POSSIBLE || col[PERCENT_COL] != PERCENT ||
       col[SPANNING_COL] != SPANNING || col[OVLEFT_COL]  != OVLEFT  ||
       col[OVRIGHT_COL]  != OVRIGHT  || col[ISIZE_COL]   != ISIZE)
      return false;

   StringVector versionCol;
   if (splitString(col[0], versionCol, ' ') != 2 || versionCol[1] == "")
      return false;

   version = versionCol[1];

   annotationHeading.clear();

   for (int i = FIRST_ANNOT_COL; i < numCols; i++)
      annotationHeading.push_back(col[i]);

   return true;
}

//------------------------------------------------------------------------------------
// writeReadPairLine() writes one line to stdout showing the total number of read
// pairs processed

void writeReadPairLine(uint64_t numReadPairs)
{
   std::cout << READ_PAIRS << numReadPairs << NEWLINE;
}

//------------------------------------------------------------------------------------
// isReadPairLine() returns true if the given line shows the total number of read
// pairs processed

static bool isReadPairLine(const std::string& line)
{
   return hasPrefix(line, READ_PAIRS);
}

//------------------------------------------------------------------------------------
// isValidReadPairLine() returns true if the given line correctly shows the total
// number of read pairs processed

static bool isValidReadPairLine(const std::string& line, uint64_t& numReadPairs)
{
   numReadPairs = UINT64_MAX;

   StringVector col, readPairCol;

   if (isReadPairLine(line) && splitString(line, col) == 1 &&
       splitString(col[0], readPairCol, ' ') == 2)
   {
      std::istringstream stream(readPairCol[1]);
      stream >> numReadPairs;
   }

   return (numReadPairs < UINT64_MAX);
}

//------------------------------------------------------------------------------------
// isPatternLine() returns true if the given line describes a pattern that was matched

static bool isPatternLine(const std::string& line)
{
   return hasPrefix(line, PATTERN);
}

//------------------------------------------------------------------------------------
// getPattern() returns a pointer to a newly-allocated HitPattern describing the
// pattern identified by the given line, or returns NULL if the input is invalid

static HitPattern *getPattern(const std::string& line)
{
   StringVector col;
   int numCols  = splitString(line, col);

   int matchingBases, possible, spanningCount, insertSize;

   if (!isPatternLine(line) || numCols <= ISIZE_COL ||
       (matchingBases  = stringToNonnegInt(col[MBASES_COL]))   <= 0 ||
       (possible       = stringToNonnegInt(col[POSSIBLE_COL])) <= 0 ||
       (spanningCount  = stringToNonnegInt(col[SPANNING_COL])) <  0 ||
       spanningCount > 2 ||
       (insertSize     = stringToNonnegInt(col[ISIZE_COL]))    <= 0)
      return NULL;

   StringVector patternCol;
   if (splitString(col[0], patternCol, ' ') != 2 || patternCol[1] == "")
      return NULL;

   const std::string& name     = patternCol[1];
   const std::string& sequence = col[SEQUENCE_COL];

   StringVector annotation;

   for (int i = FIRST_ANNOT_COL; i < numCols; i++)
      annotation.push_back(col[i]);

   return new HitPattern(name, sequence, annotation, matchingBases, possible,
                         spanningCount, insertSize);
}

//------------------------------------------------------------------------------------
// isReadLine() returns true if the given line describes one read of a read pair that
// matched a pattern

static bool isReadLine(const std::string& line)
{
   return hasPrefix(line, READ);
}

//------------------------------------------------------------------------------------
// getRead() returns a pointer to a newly-allocated HitRead describing the read
// identified by the given line, or returns NULL if the input is invalid

static HitRead *getRead(const std::string& line)
{
   StringVector col;
   int numCols = splitString(line, col);

   int matchingBases, possible, spanningCount, leftOverlap, rightOverlap;

   if (!isReadLine(line) || numCols <= OVRIGHT_COL ||
       (matchingBases  = stringToNonnegInt(col[MBASES_COL]))   <  0 ||
       (possible       = stringToNonnegInt(col[POSSIBLE_COL])) <= 0 ||
       (spanningCount  = stringToNonnegInt(col[SPANNING_COL])) <  0 ||
       spanningCount > 1 ||
       (leftOverlap    = stringToNonnegInt(col[OVLEFT_COL]))   <  0 ||
       (rightOverlap   = stringToNonnegInt(col[OVRIGHT_COL]))  <  0)
      return NULL;

   bool isSpanning = (spanningCount == 1);

   StringVector readCol;
   if (splitString(col[0], readCol, ' ') != 2 || readCol[1] == "")
      return NULL;

   const std::string& name     = readCol[1];
   const std::string& sequence = col[SEQUENCE_COL];

   int leadingBlanks = sequence.find_first_not_of(' '); // find first non-blank
   if (leadingBlanks == std::string::npos ||
       sequence.length() != leadingBlanks + possible)
      return NULL;

   std::string trimmedSequence = sequence.substr(leadingBlanks);

   return new HitRead(name, leadingBlanks, trimmedSequence, matchingBases,
                      isSpanning, leftOverlap, rightOverlap);
}

//------------------------------------------------------------------------------------
// getHit() returns a pointer to a newly-allocated Hit describing the hit identified
// by the given line (the pattern) and the next two lines (the read pair), or returns
// NULL if the input is invalid

static Hit *getHit(std::istream& istream, const std::string& patternLine)
{
   HitPattern *hitPattern = NULL;
   HitRead    *hitRead1   = NULL, *hitRead2;

   std::string readLine1, readLine2;

   if ((hitPattern = getPattern(patternLine)) &&
       getline(istream, readLine1) && (hitRead1 = getRead(readLine1)) &&
       getline(istream, readLine2) && (hitRead2 = getRead(readLine2)))
      return new Hit(hitPattern, hitRead1, hitRead2);

   delete hitPattern;
   delete hitRead1;

   return NULL;
}

//------------------------------------------------------------------------------------
// HitCompare defines the sort order of hits; duplicate hits are placed consecutively
// in the ordering

struct HitCompare
{
   bool operator()(Hit* const& a, Hit* const& b) const
   {
      // sort by ascending pattern name,
      // then by ascending number of left bases,
      // then by ascending number of right bases,
      // then by descending number of junction-spanning reads,
      // then by ascending read1 name

      int key1 = std::strcmp(a->pattern->name.c_str(), b->pattern->name.c_str());
      if (key1 != 0)
         return (key1 < 0);

      int key2 = a->pattern->leftBases - b->pattern->leftBases;
      if (key2 != 0)
         return (key2 < 0);

      int key3 = a->pattern->rightBases - b->pattern->rightBases;
      if (key3 != 0)
         return (key3 < 0);

      int key4 = a->pattern->spanningCount - b->pattern->spanningCount;
      if (key4 != 0)
         return (key4 > 0);

      return (a->read1->name < b->read1->name);
   }
};

//------------------------------------------------------------------------------------
// sortHits() sorts the hits in the given vector

static void sortHits(HitVector& hitVector)
{
   std::sort(hitVector.begin(), hitVector.end(), HitCompare());
}

//------------------------------------------------------------------------------------
// markDuplicates() marks duplicate hits in the given vector, which is assumed to be
// sorted

static void markDuplicates(HitVector& hitVector)
{
   int numHits = hitVector.size();

   for (int i = 1; i < numHits; i++)
      if (hitVector[i]->sameAs(*hitVector[i - 1]))
         hitVector[i]->duplicate = true;
}

//------------------------------------------------------------------------------------
// readHits() reads hits and stores them in sorted order in hitVector; duplicate hits
// are marked; the return value is the total number of read pairs processed

uint64_t readHits(std::istream& istream, std::string& version,
                  StringVector& annotationHeading, HitVector& hitVector)
{
   uint64_t totalReadPairs = 0;
   std::string headingLine, line;

   if (!getline(istream, headingLine))
      throw std::runtime_error("no input");

   if (!isValidHeadingLine(headingLine, version, annotationHeading))
      throw std::runtime_error("unexpected heading line");

   while (getline(istream, line))
      if (isHeadingLine(line)) // found another heading line
      {
         if (line != headingLine)
            throw std::runtime_error("inconsistent heading lines");
      }
      else if (isReadPairLine(line))
      {
         uint64_t numReadPairs;

	 if (isValidReadPairLine(line, numReadPairs))
            totalReadPairs += numReadPairs;
	 else
            throw std::runtime_error("unexpected input line: " + line);
      }
      else // this must be the first of three lines describing a hit
      {
         Hit *hit = getHit(istream, line);

	 if (hit)
            if (hitVector.size() < MAX_HITS)
               hitVector.push_back(hit);
	    else
               throw std::runtime_error("too many hits");
	 else
            throw std::runtime_error("unexpected hit format: " + line);
      }

   if (hitVector.size() > 1)
   {
      sortHits(hitVector);
      markDuplicates(hitVector);
   }

   return totalReadPairs;
}

//------------------------------------------------------------------------------------
// getPatternIndices() determines the index of the first hit of each pattern; it is
// assumed that hitVector has been sorted

void getPatternIndices(const HitVector& hitVector, IntVector& index)
{
   index.clear();

   int numHits = hitVector.size();
   if (numHits == 0)
      return;

   index.push_back(0);

   for (int i = 1; i < numHits; i++)
      if (hitVector[i]->pattern->name != hitVector[i - 1]->pattern->name)
         index.push_back(i);
}

//------------------------------------------------------------------------------------
// maxDisplayLength() returns the maximum length of the pattern display sequence for
// the hits in the given vector, from indices begin (inclusive) to end (exclusive)

int maxDisplayLength(const HitVector& hitVector, int begin, int end)
{
   int maxLength = 0;

   for (int i = begin; i < end; i++)
   {
      int length = hitVector[i]->pattern->displaySequence.length();
      if (length > maxLength)
         maxLength = length;
   }

   return maxLength;
}
