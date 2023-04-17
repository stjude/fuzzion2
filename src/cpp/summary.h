//------------------------------------------------------------------------------------
//
// summary.h - module for reading and writing fuzzion2 hit summaries
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef SUMMARY_H
#define SUMMARY_H

#include "hit.h"

const std::string FUZZUM = "fuzzum ";

//------------------------------------------------------------------------------------

class Summary
{
public:
   Summary(const std::string& inSampleID, int inReadPairs, int inWeak,
           int inStrongNospan, int inStrongSpan, const std::string& inName,
	   const StringVector& inAnnotation)
      : sampleID(inSampleID), readPairs(inReadPairs), weak(inWeak),
	strongNospan(inStrongNospan), strongSpan(inStrongSpan), name(inName),
	annotation(inAnnotation) { }

   virtual ~Summary() { }

   int distinct() const { return weak + strongNospan + strongSpan; }

   void write() const;

   std::string sampleID;    // identifies the sample
   int readPairs;           // #read pairs matching the pattern or group
   int weak;                // #distinct read pairs that are weak matches
   int strongNospan;        // #distinct read pairs that are strong/non-spanning
   int strongSpan;          // #distinct read pairs that are strong/spanning
   std::string name;        // name of pattern or group
   StringVector annotation; // pattern or group annotations
};

typedef std::vector<Summary *> SummaryVector;

//------------------------------------------------------------------------------------

void writeSummaryHeadingLine(const std::string& version, bool grouping,
                             const StringVector& annotationHeading);

Summary *summarizeHits(const HitVector& hitVector, int begin, int end, int minStrong,
                       std::string sampleID="");

void readSummaries(const StringVector& filename, StringVector& annotationHeading,
                   SummaryVector& summaryVector);

#endif
