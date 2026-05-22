//------------------------------------------------------------------------------------
//
// summary.h - module for reading and writing fuzzion2 hit summaries
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef SUMMARY_H
#define SUMMARY_H

#include "hit.h"

const String FUZZUM = "fuzzum ";

//------------------------------------------------------------------------------------

class Summary
{
public:
   Summary(const String& inSampleID, const int inNumMatches, const int inWeak,
           const int inStrongNospan, const int inStrongSpan, const String& inName,
           const StringVector& inAnnotation)
      : sampleID(inSampleID), numMatches(inNumMatches), weak(inWeak),
        strongNospan(inStrongNospan), strongSpan(inStrongSpan), name(inName),
        annotation(inAnnotation) { }

   virtual ~Summary() { }

   int distinct() const { return weak + strongNospan + strongSpan; }

   void write(std::ostream& ostream) const;

   const String sampleID;         // identifies the sample
   const int numMatches;          // total #matches to the pattern or group
   const int weak;                // #distinct matches that are weak matches
   const int strongNospan;        // #distinct matches that are strong/non-spanning
   const int strongSpan;          // #distinct matches that are strong/spanning
   const String name;             // name of pattern or group
   const StringVector annotation; // pattern or group annotations
};

typedef std::vector<Summary *> SummaryVector;

//------------------------------------------------------------------------------------

Summary *summarizeHits(const HitVector& hitVector, int begin, int end, int minStrong,
                       const String& sampleID="");

void writeSummaryHeadingLine(std::ostream& ostream, const String& version,
                             bool grouping, const StringVector& annotationHeading);

void readSummaries(const StringVector& filename, StringVector& annotationHeading,
                   SummaryVector& summaryVector);

#endif
