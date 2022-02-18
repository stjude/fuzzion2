//------------------------------------------------------------------------------------
//
// pattern.h - module supporting pattern sequences
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef PATTERN_H
#define PATTERN_H

#include "minimizer.h"
#include "util.h"
#include <unordered_map>

//------------------------------------------------------------------------------------

class Pattern
{
public:
   Pattern(const std::string& inName, const std::string& inSequence,
           const StringVector& inAnnotation);

   virtual ~Pattern() { }

   std::string name;            // pattern name
   std::string sequence;        // sequence w/o brackets or braces
   std::string displaySequence; // sequence w/  brackets or braces

   bool hasBraces;              // false = brackets: left]mid[right
                                // true  = braces  : left}mid{right
   int  delim2;                 // offset of second delimiter: [ or {
   int  leftBases;              // #bases in left side of sequence
   int  middleBases;            // #bases between delimiters
   int  rightBases;             // #bases in right side of sequence

   StringVector annotation;     // zero or more annotations
};

typedef std::vector<Pattern> PatternVector;

//------------------------------------------------------------------------------------

class Location // indicates a location within a pattern sequence
{
public:
   Location(int inIndex, int inOffset)
      : index(inIndex), offset(inOffset) { }

   virtual ~Location() { }

   int index;  // index of pattern in a PatternVector
   int offset; // offset within the pattern sequence
};

typedef std::vector<Location> LocationVector;

// this maps a minimizer to locations of the minimizer in patterns
typedef std::unordered_map<Minimizer, LocationVector> PatternMap;

//------------------------------------------------------------------------------------

bool hasDelimiters(const std::string& sequence, bool& hasBraces, int& delim2,
                   int& leftBases, int& middleBases, int& rightBases);

PatternVector *readPatterns(const std::string& filename,
                            StringVector& annotationHeading);

PatternMap *createPatternMap(const PatternVector *patternVector,
                             MinimizerWindowLength w, const KmerRankTable *rankTable,
			     Minimizer maxMinimizer);

#endif
