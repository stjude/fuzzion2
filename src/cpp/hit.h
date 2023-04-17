//------------------------------------------------------------------------------------
//
// hit.h - module for reading and writing fuzzion2 hits
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef HIT_H
#define HIT_H

#include "pattern.h"
#include <istream>
#include <limits>

const std::string FUZZION2 = "fuzzion2 ";

const std::string HIT_DUPLICATE     = "dup";
const std::string HIT_WEAK          = "weak";
const std::string HIT_STRONG_NOSPAN = "strong-";
const std::string HIT_STRONG_SPAN   = "strong+";

const int MAX_HITS = std::numeric_limits<int>::max();

const int DEFAULT_MIN_STRONG = 15; // default minimum overlap for a strong match

//------------------------------------------------------------------------------------

class HitPattern : public Pattern
{
public:
   HitPattern(const std::string& inName, const std::string& inSequence,
              const StringVector& inAnnotation, int inMatchingBases, int inPossible,
	      int inSpanningCount, int inInsertSize)
      : Pattern(inName, inSequence, inAnnotation), matchingBases(inMatchingBases),
        possible(inPossible), spanningCount(inSpanningCount),
	insertSize(inInsertSize) { }

   virtual ~HitPattern() { }

   double percentMatch() const { return (100.0 * matchingBases / possible); }

   bool isStrong(int minStrong) const
   { return (leftBases >= minStrong && rightBases >= minStrong); }

   bool isSpanning() const { return (spanningCount > 0); }

   void write() const;

   int matchingBases; // total number of matching bases in read1 and read2
   int possible;      // possible number of matching bases
   int spanningCount; // number of junction-spanning reads (0, 1 or 2)
   int insertSize;    // insert size of read pair aligned to this pattern
};

//------------------------------------------------------------------------------------

class HitRead
{
public:
   HitRead(const std::string& inName, int inLeadingBlanks,
           const std::string& inSequence, int inMatchingBases, bool inIsSpanning)
      : name(inName), leadingBlanks(inLeadingBlanks), sequence(inSequence),
	matchingBases(inMatchingBases), isSpanning(inIsSpanning) { }

   virtual ~HitRead() { }

   int possible() const { return sequence.length(); }

   double percentMatch() const { return (100.0 * matchingBases / possible()); }

   void write() const;

   std::string name;     // read name
   int  leadingBlanks;   // number of blanks preceding read sequence in the display
   std::string sequence; // read sequence
   int  matchingBases;   // number of matching bases (zero for an unmatched mate)
   bool isSpanning;      // true if this is a junction-spanning read
};

//------------------------------------------------------------------------------------

class Hit
{
public:
   Hit(HitPattern *inPattern, HitRead *inRead1, HitRead *inRead2)
      : pattern(inPattern), read1(inRead1), read2(inRead2), duplicate(false) { }

   virtual ~Hit() { delete pattern; delete read1; delete read2; }

   bool sameAs(const Hit& other) const;
   bool isStrong(int minStrong)  const { return pattern->isStrong(minStrong); }
   bool isSpanning()             const { return pattern->isSpanning(); }

   std::string label(int minStrong) const;

   void write() const { pattern->write(); read1->write(); read2->write(); }

   HitPattern *pattern; // describes the pattern that was matched by a read pair
   HitRead    *read1;   // describes the first  read of the read pair
   HitRead    *read2;   // describes the second read of the read pair
   bool duplicate;      // true if this hit is a duplicate of another hit
};

typedef std::vector<Hit *> HitVector;

//------------------------------------------------------------------------------------

void writeHitHeadingLine(const std::string& version,
                         const StringVector& annotationHeading);

void writeReadPairLine(uint64_t numReadPairs);

uint64_t readHits(std::istream& istream, std::string& version,
		  StringVector& annotationHeading, HitVector& hitVector);

void getPatternIndices(const HitVector& hitVector, IntVector& index);

int maxDisplayLength(const HitVector& hitVector, int begin, int end);

#endif
