//------------------------------------------------------------------------------------
//
// summary.h - module for reading and writing fuzzion2 hit summaries
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef SUMMARY_H
#define SUMMARY_H

#include "util.h"

const std::string FUZZUM = "fuzzum ";

//------------------------------------------------------------------------------------

class Summary
{
public:
   Summary(const std::string& inPatternName, const std::string& inSampleID,
           int inReadPairs, int inDistinct, int inStrong,
	   const StringVector& inAnnotation)
      : patternName(inPatternName), sampleID(inSampleID), readPairs(inReadPairs),
	distinct(inDistinct), strong(inStrong), annotation(inAnnotation) { }

   virtual ~Summary() { }

   void write() const;

   std::string patternName; // name of the pattern
   std::string sampleID;    // identifies the sample
   int readPairs;           // #read pairs from this sample matching the pattern
   int distinct;            // #distinct read pairs that match
   int strong;              // #distinct read pairs that are strong matches
   StringVector annotation; // pattern annotations
};

typedef std::vector<Summary *> SummaryVector;

//------------------------------------------------------------------------------------

void writeSummaryHeadingLine(const std::string& version,
                             const StringVector& annotationHeading);

void readSummaries(const StringVector& filename, StringVector& annotationHeading,
                   SummaryVector& summaryVector);

#endif
