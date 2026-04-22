//------------------------------------------------------------------------------------
//
// align.h - module for aligning substrings of sequences
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef ALIGN_H
#define ALIGN_H

#include "pattern.h"
#include "window.h"

const char SPACER = '-'; // spacer symbol used in alignment visualization

//------------------------------------------------------------------------------------
// An "origin" is an offset in sequence1 that corresponds to the first base of
// sequence2 in an alignment.  It can be negative.  See examples below.
//
// Adding the origin to a sequence2 offset gives the corresponding sequence1 offset.
// Subtracting the origin from a sequence1 offset gives the corresponding sequence2
// offset.
//                  0123456789...
//     sequence1    ACCGTATGAA...    origin is 4
//     sequence2        TATAAA...
//
//                  0123456789...
//     sequence1    TCAGAGAAGG...    origin is 0
//     sequence2    TCCGAGATGG...
//
//                     0123456...
//     sequence1       TGCAAGA...    origin is -3
//     sequence2    CGATCCAAGA...

typedef int Origin;

inline Origin getOrigin(const int offset1, const int offset2)
{ return (offset1 - offset2); } // returns origin given sequence offsets that align

typedef std::vector<Origin> OriginVector;
typedef std::unordered_map<Origin, int> OriginCountMap; // origin, count

class OriginCounter // counts occurrences of each origin
{
public:
   OriginCounter()
      : ocmap() { }

   virtual ~OriginCounter() { }

   void incrementCount(Origin origin);
   int  maxCount() const;
   int  retrieveOrigins(int minCount, OriginVector& ovector) const;

   OriginCountMap ocmap;
};

typedef std::array<OriginCounter, 2> OriginCounterDuo;    // index 0 is left/unsided;
                                                          // 1 is right
typedef std::map<int, OriginCounterDuo> PatternOriginMap; // key is pattern index

//------------------------------------------------------------------------------------

class Sub // a substring of a sequence
{
public:
   Sub(const Seq& inSeq)
      : seq(inSeq), begin(0), len(inSeq.len) { }
   Sub(const Seq& inSeq, const int inBegin, const int inLen)
      : seq(inSeq), begin(inBegin), len(inLen) { }
   Sub(const Sub& other, bool before);

   virtual ~Sub() { }

   String str()       const { return seq.str.substr(begin, len); }
   const char *cstr() const { return &seq.cstr[begin]; }

   bool hasBefore()   const { return (begin > 0); }
   bool hasAfter()    const { return (seq.len > begin + len); }

   const Seq& seq; // the sequence

   int begin; // first offset of substring in sequence
   int len;   // substring length >= 0
   int end() const { return (begin + len); } // exclusive last offset
};

//------------------------------------------------------------------------------------

class Align // an alignment of two substrings
{
public:
   Align(const Sub& inSub1, const Sub& inSub2,
         const String& inVis1, const String& inVis2)
      : sub1(inSub1), sub2(inSub2), vis1(inVis1), vis2(inVis2) { }

   virtual ~Align() { }

   Origin origin() const { return getOrigin(sub1.begin, sub2.begin); }

   void adjust(int lbases1, int lbases2, int rbases1, int rbases2);
   void getStats(int& possible, int& matches) const;

   Sub    sub1, sub2; // the two substrings
   String vis1, vis2; // visualization of alignment
};

typedef std::vector<Align *> AlignVector;
typedef std::array<AlignVector, 2> AlignDuo; // index 0 is left/unsided; 1 is right
typedef std::map<int, AlignDuo> PatternAlignMap; // key is pattern index

//------------------------------------------------------------------------------------

void connectReadToPatterns(const WindowVector& readWindowVector,
                           Minimizer maxMinimizer, const PatternMap& pmap,
			   const BoolVector& inPmap, PatternOriginMap& omap,
			   const BoolVector *eligible=nullptr);

void getOverlap(const Sub& sub1, const Sub& sub2, Origin origin,
                int& begin1, int& begin2, int& len, int& lbases1, int& lbases2,
		int& rbases1, int& rbases2);

inline int computeMinMatches(const int len, const double minPercentAgreement)
{ return std::ceil((minPercentAgreement / 100) * len); }

Align *alignSubstrings(const Sub& sub1, const Sub& sub2, int minMatches);

void alignReadToPatterns(const Seq& readSeq, const PatternVector& pvector,
                         double minPercentAgreement, const PatternOriginMap& omap,
			 PatternAlignMap& amap);

#endif
