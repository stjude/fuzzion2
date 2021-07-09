//------------------------------------------------------------------------------------
//
// pattern.h - module supporting pattern sequences
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
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
   int  leftBases;              // #bases in left  side of sequence
   int  rightBases;             // #bases in right side of sequence
   bool hasBraces;              // true = braces, false = brackets
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

PatternVector *readPatterns(const std::string& filename,
                            StringVector& annotationHeading);

PatternMap *createPatternMap(const PatternVector *patternVector,
                             MinimizerWindowLength w, const KmerRankTable *rankTable,
			     Minimizer maxMinimizer);

//------------------------------------------------------------------------------------
#endif
