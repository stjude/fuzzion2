//------------------------------------------------------------------------------------
//
// hit.cpp - module for reading and writing fuzzion2 hits
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
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
const std::string ISIZE      = "insert size";

const int SEQUENCE_COL       = 1;
const int MBASES_COL         = 2;
const int POSSIBLE_COL       = 3;
const int PERCENT_COL        = 4;
const int ISIZE_COL          = 5;
const int MIN_HEADING_COLS   = 6;
const int READ_COLS          = 5;

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
// writeHitHeadingLine() writes a heading line to stdout

void writeHitHeadingLine(const std::string& version,
                         const StringVector& annotationHeading)
{
   std::cout << FUZZION2 << version
             << TAB << SEQUENCE
	     << TAB << MBASES
	     << TAB << POSSIBLE
	     << TAB << PERCENT
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

   if (!isHeadingLine(line) || numCols < MIN_HEADING_COLS ||
       col[SEQUENCE_COL] != SEQUENCE || col[MBASES_COL]  != MBASES ||
       col[POSSIBLE_COL] != POSSIBLE || col[PERCENT_COL] != PERCENT ||
       col[ISIZE_COL]    != ISIZE)
      return false;

   StringVector versionCol;
   if (splitString(col[0], versionCol, ' ') != 2 || versionCol[1] == "")
      return false;

   version = versionCol[1];

   annotationHeading.clear();

   for (int i = MIN_HEADING_COLS; i < numCols; i++)
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
// isValidPatternLine() returns true if the given line correctly describes a pattern
// that was matched

static bool isValidPatternLine(const std::string& line, std::string& name,
                               std::string& sequence, StringVector& annotation,
			       int& matchingBases, int& possible, int& insertSize)
{
   if (!isPatternLine(line))
      return false;

   StringVector col;
   int numCols = splitString(line, col);

   if (numCols < MIN_HEADING_COLS ||
       (matchingBases = stringToNonnegInt(col[MBASES_COL]))   <= 0 ||
       (possible      = stringToNonnegInt(col[POSSIBLE_COL])) <= 0 ||
       (insertSize    = stringToNonnegInt(col[ISIZE_COL]))    <= 0)
      return false;

   StringVector patternCol;
   if (splitString(col[0], patternCol, ' ') != 2 || patternCol[1] == "")
      return false;

   name     = patternCol[1];
   sequence = col[SEQUENCE_COL];

   annotation.clear();

   for (int i = MIN_HEADING_COLS; i < numCols; i++)
      annotation.push_back(col[i]);

   return true;
}

//------------------------------------------------------------------------------------
// isReadLine() returns true if the given line describes one read of a read pair that
// matched a pattern

static bool isReadLine(const std::string& line)
{
   return hasPrefix(line, READ);
}

//------------------------------------------------------------------------------------
// isValidReadLine() returns true if the given line correctly describes one read of a
// read pair that matched a pattern

static bool isValidReadLine(const std::string& line, std::string& name,
                            int& leadingBlanks, std::string& sequence,
			    int& matchingBases)
{
   StringVector col, readCol;

   if (!isReadLine(line) || splitString(line, col) != READ_COLS ||
       splitString(col[0], readCol, ' ') != 2 || readCol[1] == "" ||
       (leadingBlanks = col[SEQUENCE_COL].find_first_not_of(' ')) ==
       std::string::npos || col[SEQUENCE_COL].length() !=
       leadingBlanks + stringToNonnegInt(col[POSSIBLE_COL]) ||
       (matchingBases = stringToNonnegInt(col[MBASES_COL])) < 0)
      return false;

   name     = readCol[1];
   sequence = col[SEQUENCE_COL].substr(leadingBlanks);

   return true;
}

//------------------------------------------------------------------------------------
// getHit() returns a pointer to a newly-allocated Hit describing the hit read from
// the given line (the pattern) and the next two lines (the read pair), or returns
// NULL if the input is invalid

static Hit *getHit(std::istream& istream, const std::string& line)
{
   std::string name, sequence, readLine;
   StringVector annotation;
   int matchingBases, possible, insertSize, leadingBlanks;

   if (!isValidPatternLine(line, name, sequence, annotation, matchingBases, possible,
                           insertSize))
      return NULL;

   HitPattern *pattern = new HitPattern(name, sequence, annotation, matchingBases,
		                        possible, insertSize);

   HitRead *read[2];

   for (int i = 0; i < 2; i++)
      if (getline(istream, readLine) &&
          isValidReadLine(readLine, name, leadingBlanks, sequence, matchingBases))
         read[i] = new HitRead(name, leadingBlanks, sequence, matchingBases);
      else
      {
         delete pattern;
	 if (i == 1) delete read[0];
	 return NULL;
      }

   return new Hit(pattern, read[0], read[1]);
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
// readHits() reads hits and stores them in sorted order in hitVector; the return
// value is the total number of read pairs processed by fuzzion2

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
      else
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
      sortHits(hitVector);

   return totalReadPairs;
}

//------------------------------------------------------------------------------------
// getPatternIndices() determines the index of the first hit of each pattern; it is
// assumed that hitVector has been sorted

void getPatternIndices(const HitVector& hitVector, std::vector<int>& index)
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
// getDistinctIndices() determines the index of each distinct hit in hitVector from
// begin (inclusive) to end (exclusive); it is assumed that hitVector has been sorted

void getDistinctIndices(const HitVector& hitVector, int begin, int end,
                        std::vector<int>& index)
{
   index.clear();

   if (begin >= end)
      return;

   index.push_back(begin);

   for (int i = begin + 1; i < end; i++)
      if (!hitVector[i]->sameAs(*hitVector[i - 1]))
         index.push_back(i);
}
