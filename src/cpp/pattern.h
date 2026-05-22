//------------------------------------------------------------------------------------
//
// pattern.h - module supporting patterns
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------
// Seven types of patterns are supported:
//
// Description                  Form       left  middle  right  xvector
// --------------------    --------------  ----  ------  -----  -------
// sequence w/o delims        sequence      yes    no     no      no
//
// fusion with middle      left]mid[right   yes    yes    yes     no
//    ITD with middle      left}mid{right   yes    yes    yes     no
//
// fusion with wildcard     left]*[right    yes    no     yes     no
//    ITD with wildcard     left}*{right    yes    no     yes     no
//
// hotspot with inclusion  left<in>right    yes    no     yes     yes
// hotspot with exclusion  left(ex)right    yes    no     yes     yes
//------------------------------------------------------------------------------------

#ifndef PATTERN_H
#define PATTERN_H

#include "minimizer.h"
#include <array>
#include <unordered_map>

const KmerLength TRIMER_LEN  =  3; // k = 3
const int        NUM_TRIMERS = 64; // there are 64 3-mers

typedef std::array<IntVector, NUM_TRIMERS> TrimerOffsets;

void initTrimerOffsets(const char *cstr, int begin, int len,
                       TrimerOffsets& trimerOffsets);

//------------------------------------------------------------------------------------

class Seq
{
public:
   Seq(const String& inStr, bool& invalidBase);
   Seq(const Seq& other);

   virtual ~Seq() { delete[] cstr; }

   String  str;   // original string
   char   *cstr;  // C-style string with all-uppercase ACGT used in comparisons
   const int len; // string length >= 0
};

typedef std::vector<Seq> SeqVector;

//------------------------------------------------------------------------------------

class Pattern
{
public:
   Pattern(const String& inName, const String& inSequence,
           const StringVector& inAnnotation);

   virtual ~Pattern() { delete left; delete middle; delete right; delete xvector; }

   const String name;             // unique pattern name
   const String sequence;         // full sequence with delimiters
   const StringVector annotation; // vector of annotations (possibly empty)

   int delim1, delim2;            // offsets of delimiters in sequence (-1 if N/A)

   Seq *left;                     // left or unsided sequence (always non-null)
   Seq *middle;                   // middle sequence (nullptr if N/A)
   Seq *right;                    // right sequence  (nullptr if N/A)

   SeqVector *xvector;            // extra sequences (nullptr if N/A)
   String xvis;                   // visualization of extra sequences

   TrimerOffsets leftTrimers, rightTrimers; // offsets of 3-mers in left & right
};

typedef std::vector<Pattern *> PatternVector;

//------------------------------------------------------------------------------------

class Ploc // a location in a pattern sequence
{
public:
   Ploc(const int inPindex, const int inPoffset)
      : pindex(inPindex), poffset(inPoffset) { }

   virtual ~Ploc() { }

   const int pindex;  // index of pattern in PatternVector
   const int poffset; // offset in pattern's left/unsided or right sequence
};

typedef std::vector<Ploc> PlocVector;
typedef std::array<PlocVector, 2> PlocDuo; // index 0 is left/unsided; 1 is right

// this maps a minimizer to locations of the minimizer in patterns
typedef std::unordered_map<Minimizer, PlocDuo> PatternMap;

//------------------------------------------------------------------------------------

PatternVector *readPatterns(const String& filename, StringVector& annotationHeading);

PatternMap *createPatternMap(PatternVector& pvector, MinimizerWindowLength w,
                             const KmerRankTable *rankTable, Minimizer maxMinimizer,
                             BoolVector *& inPatternMap);

#endif
