//------------------------------------------------------------------------------------
//
// fuzzhop.cpp - this program examines the hits in two or more fuzzion2 output files
//               and writes to stdout possible instances of index hopping
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "hit.h"
#include "version.h"

const String VERSION_NAME = "fuzzhop " + CURRENT_VERSION;

const int DEFAULT_BY_LANE = 1;
int byLane    = DEFAULT_BY_LANE;    // setting of -bylane option
int minStrong = DEFAULT_MIN_STRONG; // minimum overlap for a strong match

StringVector fuzzion2Filename;
int numFuzzion2Files;

//------------------------------------------------------------------------------------

// this maps flowcell or lane to the number of strong hits in each file
typedef std::map<String, IntVector> FlowcellMap;

class PatternHitCount
{
public:
   PatternHitCount(const StringVector& inAnnotation, const int numFiles)
      : annotation(inAnnotation), hitCount(numFiles, 0), fmap() { }

   virtual ~PatternHitCount() { }

   void addStrongHit(const String& readName, int fileIndex);

   const StringVector annotation; // pattern annotations
   IntVector   hitCount;          // holds #strong hits of a pattern in each file
   FlowcellMap fmap;              // for each flowcell or lane, #strong hits per file
}; 

typedef std::map<String, PatternHitCount> PatternHitMap; // key is pattern name

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname << " OPTION ..."
      << " fuzzion2_filename1 fuzzion2_filename2 ... > possible_index_hops"
      << NEWLINE;

   std::cerr
      << NEWLINE
      << "The following are optional:" << NEWLINE
      << "  -bylane=N   "
             << "group by flowcell lane (1) or by flowcell (0), default is "
	     << DEFAULT_BY_LANE << NEWLINE
      << "  -strong=N   "
             << "minimum overlap of a strong match in #bases,   default is "
	     << DEFAULT_MIN_STRONG << NEWLINE;
}

//------------------------------------------------------------------------------------
// parseArgs() parses the command-line arguments and returns true if all are valid

bool parseArgs(const int argc, const char *argv[])
{
   for (int i = 1; i < argc; i++)
   {
      const String arg = argv[i];
      if (arg.length() == 0)
         continue;

      if (arg[0] != '-')
      {
         fuzzion2Filename.push_back(arg);
	 continue;
      }

      StringVector opt;

      if (splitString(arg, opt, '=') != 2)
         return false; // incorrect option format

      if (intOpt(opt, "bylane", byLane) ||
          intOpt(opt, "strong", minStrong))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   numFuzzion2Files = fuzzion2Filename.size();

   return (numFuzzion2Files > 1 && (byLane == 0 || byLane == 1) &&
           minStrong >= TRIMER_LEN);
}

//------------------------------------------------------------------------------------
// getFlowcell() returns the flowcell or lane that is embedded in the given read name,
// or returns the empty string if it cannot be determined; we assume that a read name
// follows Illumina conventions and contains a series of values separated by colons,
// and that the flowcell and lane are identified by all but the last three values of
// the series

String getFlowcell(const String& readName)
{
   StringVector part;

   const int numParts = splitString(readName, part, ':');
   const int numKeep  = numParts - (byLane ? 3 : 4);

   if (numKeep < 1)
      return "";

   String flowcell = part[0];

   for (int i = 1; i < numKeep; i++)
      flowcell += ":" + part[i];

   return flowcell;
}

//------------------------------------------------------------------------------------
// PatternHitCount::addStrongHit() increments the #strong hits for the file having the
// specified index from 0 to (numFuzzion2Files - 1)

void PatternHitCount::addStrongHit(const String& readName, const int fileIndex)
{
   const String flowcell = getFlowcell(readName);
   if (flowcell == "")
      throw std::runtime_error("unable to obtain flowcell information from read name "
                               + readName + " in " + fuzzion2Filename[fileIndex]);

   hitCount[fileIndex]++; // increment the overall hit count

   // now increment the hit count for this flowcell

   FlowcellMap::iterator fpos = fmap.find(flowcell);
   if (fpos == fmap.end())
      fpos = fmap.insert(
          std::make_pair(flowcell, IntVector(numFuzzion2Files, 0))).first;

   IntVector& flowcellCount = fpos->second;
   flowcellCount[fileIndex]++;
}

//------------------------------------------------------------------------------------
// addPatternHit() increments the #strong hits for the file having the specified index
// from 0 to (numFuzzion2Files - 1); this is done for the pattern and read1 name of
// the given hit

void addPatternHit(PatternHitMap& pmap, const Hit& hit, const int fileIndex)
{
   const String& patternName = hit.patternName;

   PatternHitMap::iterator ppos = pmap.find(patternName);
   if (ppos == pmap.end())
      ppos = pmap.insert(std::make_pair(patternName,
                         PatternHitCount(hit.annotation, numFuzzion2Files))).first;

   PatternHitCount& patternHitCount = ppos->second;
   patternHitCount.addStrongHit(hit.read1->name, fileIndex);
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

      String    version;
      HitVector hitVector;
      uint64_t  numReads;

      readHits(infile, version, annotationHeading, hitVector, numReads);

      infile.close();

      const int numHits = hitVector.size();

      for (int j = 0; j < numHits; j++)
      {
         const Hit *hit = hitVector[j];

         if (hit->isStrong(minStrong))
            addPatternHit(pmap, *hit, i);

	 delete hit;
      }
   }
}

//------------------------------------------------------------------------------------
// writeHeadingLine() writes a heading line to stdout

void writeHeadingLine(const StringVector& annotationHeading)
{
   std::cout << VERSION_NAME
             << TAB << (byLane ? "flowcell lane" : "flowcell")
	     << TAB << "strong hits here"
	     << TAB << "strong hits elsewhere"
	     << TAB << "file name";

   const int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotationHeading[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// writePossibleHops() writes to stdout each detected possible instance of index
// hopping

void writePossibleHops(const PatternHitMap& pmap)
{
   for (PatternHitMap::const_iterator ppos = pmap.cbegin();
        ppos != pmap.cend(); ++ppos)
   {
      const String& patternName = ppos->first;
      const PatternHitCount& patternHitCount = ppos->second;
      const FlowcellMap& fmap = patternHitCount.fmap;

      for (FlowcellMap::const_iterator fpos = fmap.cbegin();
           fpos != fmap.cend(); ++fpos)
      {
         const String&    flowcell      = fpos->first;
	 const IntVector& flowcellCount = fpos->second;

	 IntVector positiveFileIndex;  // indices of files having strong hits

	 for (int i = 0; i < numFuzzion2Files; i++)
            if (flowcellCount[i] > 0)
               positiveFileIndex.push_back(i);

	 const int numPositive = positiveFileIndex.size();
	 if (numPositive < 2)
            continue;

	 // found possible index hopping
	 for (int j = 0; j < numPositive; j++)
	 {
            const int i = positiveFileIndex[j];

	    std::cout << patternName
                      << TAB << flowcell
		      << TAB << flowcellCount[i]
		      << TAB << patternHitCount.hitCount[i] - flowcellCount[i]
		      << TAB << fuzzion2Filename[i];

	    const int numAnnotations = patternHitCount.annotation.size();

	    for (int k = 0; k < numAnnotations; k++)
               std::cout << TAB << patternHitCount.annotation[k];

	    std::cout << NEWLINE;
	 }
      }
   }
}

//------------------------------------------------------------------------------------

int main(const int argc, const char *argv[])
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
