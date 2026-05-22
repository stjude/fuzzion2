//------------------------------------------------------------------------------------
//
// hit.h - module for reading and writing fuzzion2 hits
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef HIT_H
#define HIT_H

#include "match.h"

const String FUZZION2 = "fuzzion2 ";

const String HIT_DUPLICATE     = "dup";
const String HIT_WEAK          = "weak";
const String HIT_STRONG_NOSPAN = "strong-";
const String HIT_STRONG_SPAN   = "strong+";

const int MAX_HITS = std::numeric_limits<int>::max();

const int DEFAULT_MIN_STRONG = 15; // default minimum overlap for a strong match

//------------------------------------------------------------------------------------

class HitRead // represents one read of a fuzzion2 hit
{
public:
   HitRead(const String& inName, const String& inVis, const int inPossible,
           const int inMatches, const int inSpanning, const int inLoverlap,
           const int inRoverlap, const int inLen, const int inHash)
      : name(inName), vis(inVis), possible(inPossible), matches(inMatches),
        spanning(inSpanning), loverlap(inLoverlap), roverlap(inRoverlap), len(inLen),
        hash(inHash) { }

   virtual ~HitRead() { }

   void write(std::ostream& ostream) const;

   const String name;   // read name
   const String vis;    // visualization of read sequence
   const int possible;  // possible #matching bases
   const int matches;   // actual   #matching bases
   const int spanning;  // 1 if spanning, 0 if not
   const int loverlap;  // length of left/unsided alignment
   const int roverlap;  // length of right alignment
   const int len;       // read length
   const uint32_t hash; // hash of read sequence (for duplicate detection)
};

//------------------------------------------------------------------------------------

class Hit // represents a fuzzion2 hit
{
public:
   Hit(const String& inPatternName, const String& inPatternVis,
       int inPossible, int inMatches, int inSpanning, int inInsertSize,
       const StringVector& inAnnotation, const HitRead *inRead1,
       const HitRead *inRead2=nullptr);

   virtual ~Hit() { delete read1; delete read2; }

   bool sameAs(const Hit& other) const;
   bool isStrong(int minStrong)  const;
   String label (int minStrong)  const;

   void write(std::ostream& ostream) const;

   const String patternName;      // pattern name
   const String patternVis;       // visualization of pattern sequence
   const int delim1;              // offset of first delimiter in patternVis
   const int possible;            // possible #matching bases
   const int matches;             // actual   #matching bases
   const int spanning;            // #spanning reads (0, 1 or 2)
   const int insertSize;          // insert size of read pair or 0 if N/A
   const StringVector annotation; // pattern annotations
   const HitRead *read1;          // required read (cannot be nullptr)
   const HitRead *read2;          // optional second read (may be nullptr)
   uint64_t hash;                 // hash of read sequences (for duplicate detection)
   bool duplicate;                // true if this hit is a duplicate of another hit
};

typedef std::vector<Hit *> HitVector;

//------------------------------------------------------------------------------------

Hit *createHitFromSingleMatch(const Pattern& pattern, int maxmidlen,
                              double minPercentAgreement, int minOverlap,
                              const String& readName, const Seq& readSeq,
                              const SingleMatch& match);

Hit *createHitFromPairMatch(const Pattern& pattern, int maxmidlen,
                            double minPercentAgreement, int minOverlap,
                            const String& readName1, const Seq& readSeq1,
                            const String& readName2, const Seq& readSeq2,
                            const PairMatch& pairMatch);

void writeHitHeadingLine(std::ostream& ostream, const String& version,
                         const StringVector& annotationHeading);

void writeReadCountLine(std::ostream& ostream, uint64_t numReads);

void readHits(std::istream& istream, String& version, StringVector& annotationHeading,
              HitVector& hitVector, uint64_t& numReads);

void getPatternIndices(const HitVector& hitVector, IntVector& index);

int maxVisLength(const HitVector& hitVector, int begin, int end);

#endif
