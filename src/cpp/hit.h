//------------------------------------------------------------------------------------
//
// hit.h - module for reading and writing fuzzion2 hits
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef HIT_H
#define HIT_H

#include "pattern.h"
#include <limits>

const std::string FUZZION2 = "fuzzion2 ";

const int MAX_HITS = std::numeric_limits<int>::max();

const int DEFAULT_MIN_STRONG = 15; // default minimum overlap for a strong match

//------------------------------------------------------------------------------------

class HitPattern : public Pattern
{
public:
   HitPattern(const std::string& inName, const std::string& inSequence,
              const StringVector& inAnnotation, int inMatchingBases, int inPossible,
	      int inInsertSize)
      : Pattern(inName, inSequence, inAnnotation), matchingBases(inMatchingBases),
        possible(inPossible), insertSize(inInsertSize) { }

   virtual ~HitPattern() { }

   double percentMatch() const { return (100.0 * matchingBases / possible); }

   bool isStrong(int minStrong) const
   { return (leftBases >= minStrong && rightBases >= minStrong); }

   void write() const;

   int matchingBases; // total number of matching bases in read1 and read2
   int possible;      // possible number of matching bases
   int insertSize;    // insert size of read pair aligned to this pattern
};

//------------------------------------------------------------------------------------

class HitRead
{
public:
   HitRead(const std::string& inName, int inLeadingBlanks,
           const std::string& inSequence, int inMatchingBases)
      : name(inName), leadingBlanks(inLeadingBlanks), sequence(inSequence),
	matchingBases(inMatchingBases) { }

   virtual ~HitRead() { }

   int possible() const { return sequence.length(); }

   double percentMatch() const { return (100.0 * matchingBases / possible()); }

   void write() const;

   std::string name;     // read name
   int leadingBlanks;    // number of blanks preceding read sequence in the display
   std::string sequence; // read sequence
   int matchingBases;    // number of matching bases (zero for an unmatched mate)
};

//------------------------------------------------------------------------------------

class Hit
{
public:
   Hit(HitPattern *inPattern, HitRead *inRead1, HitRead *inRead2)
      : pattern(inPattern), read1(inRead1), read2(inRead2) { }

   virtual ~Hit() { delete pattern; delete read1; delete read2; }

   bool sameAs(const Hit& other) const;
   bool isStrong(int minStrong)  const { return pattern->isStrong(minStrong); }

   void write() const { pattern->write(); read1->write(); read2->write(); }

   HitPattern *pattern; // describes the pattern that was matched by a read pair
   HitRead    *read1;   // describes the first  read of the read pair
   HitRead    *read2;   // describes the second read of the read pair
};

typedef std::vector<Hit *> HitVector;

//------------------------------------------------------------------------------------

void writeHitHeadingLine(const std::string& version,
                         const StringVector& annotationHeading);

void writeReadPairLine(uint64_t numReadPairs);

uint64_t readHits(std::string& version, StringVector& annotationHeading,
                  HitVector& hitVector);

void getPatternIndices(const HitVector& hitVector, std::vector<int>& index);

void getDistinctIndices(const HitVector& hitVector, int begin, int end,
                        std::vector<int>& index);

#endif
