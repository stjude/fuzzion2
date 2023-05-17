//------------------------------------------------------------------------------------
//
// fuzzhop.cpp - this program examines the hits in two or more fuzzion2 output files
//               and writes to stdout possible instances of index hopping
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "hit.h"
#include "version.h"
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

const std::string VERSION_NAME = "fuzzhop " + CURRENT_VERSION;

StringVector fuzzion2Filename;
int numFuzzion2Files;

//------------------------------------------------------------------------------------

typedef std::map<std::string, IntVector> FlowcellLaneMap; // flowcell lane,
                                                          // #hits in each file

class PatternHitCount
{
public:
   PatternHitCount(const StringVector& inAnnotation, int numFiles)
      : annotation(inAnnotation), hitCount(numFiles, 0), fmap() { }

   virtual ~PatternHitCount() { }

   void addHit(const std::string& readName, int fileIndex);

   StringVector annotation; // pattern annotations
   IntVector hitCount;      // holds #hits of a pattern in each file
   FlowcellLaneMap fmap;    // for each flowcell lane, holds #hits in each file
};

typedef std::map<std::string, PatternHitCount> PatternHitMap; // key is pattern name

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname
      << " fuzzion2_filename1 fuzzion2_filename2 ... > possible_index_hops"
      << NEWLINE;
}

//------------------------------------------------------------------------------------
// parseArgs() parses the command-line arguments and returns true if all are valid

bool parseArgs(int argc, char *argv[])
{
   for (int i = 1; i < argc; i++)
   {
      std::string arg = argv[i];
      if (arg.length() == 0)
         continue;

      if (arg[0] != '-')
      {
         fuzzion2Filename.push_back(arg);
	 continue;
      }

      return false; // unrecognized option
   }

   numFuzzion2Files = fuzzion2Filename.size();

   return (numFuzzion2Files > 1);
}

//------------------------------------------------------------------------------------
// getFlowcellLane() returns the flowcell lane that is embedded in the given read
// name, or returns the empty string if the flowcell lane cannot be determined; it is
// assumed that a read name contains a series of values separated by colons, and that
// the flowcell and lane are identified by all but the last three values of the series

std::string getFlowcellLane(const std::string& readName)
{
   StringVector part;

   int numParts = splitString(readName, part, ':');
   if (numParts < 4)
      return "";

   int numKeep = numParts - 3;

   std::string flowcellLane = part[0];

   for (int i = 1; i < numKeep; i++)
      flowcellLane += ":" + part[i];

   return flowcellLane;
}

//------------------------------------------------------------------------------------
// PatternHitCount::addHit() increments the #hits for the file having the specified
// index from 0 to (numFuzzion2Files - 1)

void PatternHitCount::addHit(const std::string& readName, int fileIndex)
{
   std::string flowcellLane = getFlowcellLane(readName);
   if (flowcellLane == "")
      throw std::runtime_error("unable to obtain flowcell lane from read name " +
                               readName + " in " + fuzzion2Filename[fileIndex]);

   hitCount[fileIndex]++; // increment the overall hit count

   // now increment the hit count for this flowcell and lane

   FlowcellLaneMap::iterator fpos = fmap.find(flowcellLane);
   if (fpos == fmap.end())
   {
      fmap.insert(std::make_pair(flowcellLane, IntVector(numFuzzion2Files, 0)));
      fpos = fmap.find(flowcellLane);
   }

   IntVector& flowcellLaneCount = fpos->second;
   flowcellLaneCount[fileIndex]++;
}

//------------------------------------------------------------------------------------
// addPatternHit() increments the #hits for the file having the specified index from
// 0 to (numFuzzion2Files - 1); this is done for the pattern and read name of the
// given hit

void addPatternHit(PatternHitMap& pmap, const Hit *hit, int fileIndex)
{
   const std::string& patternName = hit->pattern->name;

   PatternHitMap::iterator ppos = pmap.find(patternName);
   if (ppos == pmap.end())
   {
      pmap.insert(std::make_pair(patternName,
                  PatternHitCount(hit->pattern->annotation, numFuzzion2Files)));
      ppos = pmap.find(patternName);
   }

   PatternHitCount& patternHitCount = ppos->second;
   patternHitCount.addHit(hit->read1->name, fileIndex);
}

//------------------------------------------------------------------------------------
// initializePatternHitMap() reads the hits from the input files and updates the
// PatternHitMap

void initializePatternHitMap(PatternHitMap& pmap, StringVector& annotationHeading)
{
   for (int i = 0; i < numFuzzion2Files; i++) // i is the file index
   {
      std::ifstream infile(fuzzion2Filename[i].c_str());
      if (!infile.is_open())
         throw std::runtime_error("unable to open " + fuzzion2Filename[i]);

      std::string version;
      HitVector   hitVector;

      readHits(infile, version, annotationHeading, hitVector);

      infile.close();

      int numHits = hitVector.size();

      for (int j = 0; j < numHits; j++)
      {
         Hit *hit = hitVector[j];
	 addPatternHit(pmap, hit, i);
	 delete hit;
      }
   }
}

//------------------------------------------------------------------------------------
// writeHeadingLine() writes a heading line to stdout

void writeHeadingLine(const StringVector& annotationHeading)
{
   std::cout << VERSION_NAME
             << TAB << "flowcell lane"
	     << TAB << "read pairs"
	     << TAB << "other read pairs"
	     << TAB << "file name";

   int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotationHeading[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// writePossibleHops() writes to stdout each detected possible instance of index
// hopping

void writePossibleHops(PatternHitMap& pmap)
{
   for (PatternHitMap::iterator ppos = pmap.begin(); ppos != pmap.end(); ++ppos)
   {
      const std::string& patternName     = ppos->first;
      PatternHitCount&   patternHitCount = ppos->second;
      FlowcellLaneMap&   fmap            = patternHitCount.fmap;

      for (FlowcellLaneMap::iterator fpos = fmap.begin(); fpos != fmap.end(); ++fpos)
      {
         const std::string& flowcellLane      = fpos->first;
	 const IntVector&   flowcellLaneCount = fpos->second;

	 IntVector positiveFileIndex;  // indices of files having hits

	 for (int i = 0; i < numFuzzion2Files; i++)
            if (flowcellLaneCount[i] > 0)
               positiveFileIndex.push_back(i);

	 int numPositive = positiveFileIndex.size();
	 if (numPositive < 2)
            continue;

	 // found possible index hopping
	 for (int j = 0; j < numPositive; j++)
	 {
            int i = positiveFileIndex[j];

	    std::cout << patternName
                      << TAB << flowcellLane
		      << TAB << flowcellLaneCount[i]
		      << TAB << patternHitCount.hitCount[i] - flowcellLaneCount[i]
		      << TAB << fuzzion2Filename[i];

	    int numAnnotations = patternHitCount.annotation.size();

	    for (int k = 0; k < numAnnotations; k++)
               std::cout << TAB << patternHitCount.annotation[k];

	    std::cout << NEWLINE;
	 }
      }
   }
}

//------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   if (!parseArgs(argc, argv))
   {
      showUsage(argv[0]);
      return 1;
   }

   try
   {
      PatternHitMap pmap;
      StringVector  annotationHeading;

      initializePatternHitMap(pmap, annotationHeading);

      writeHeadingLine(annotationHeading);
      writePossibleHops(pmap);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
