//------------------------------------------------------------------------------------
//
// fuzzall.cpp - this program reads summaries produced by fuzzum and writes a pattern
//               summary to stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "util.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

const std::string VERSION_NAME = "fuzzall v1.0.2";
const std::string COPYRIGHT    =
   "copyright 2021 St. Jude Children's Research Hospital";

const std::string FUZZUM       = "fuzzum ";
const std::string READ_PAIRS   = "read pairs"; const int RDPAIR_COL  = 1;
const std::string PATTERN      = "pattern";    const int PATTERN_COL = 2;

const int MIN_HEADING_COLS     = 3;

StringVector fuzzumFilename;

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname << " fuzzum_filename ... > pattern_summary"
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
         fuzzumFilename.push_back(arg);
	 continue;
      }

      return false; // unrecognized option
   }

   return (fuzzumFilename.size() > 0);
}

//------------------------------------------------------------------------------------

typedef std::map<std::string, uint64_t> IdMap; // id, #read pairs

class Pattern
{
public:
   Pattern(const std::string& inName, const StringVector& inAnnotation)
      : name(inName), annotation(inAnnotation), idmap() { }

   virtual ~Pattern() { }

   void write();

   std::string  name;
   StringVector annotation;
   IdMap        idmap;
};

typedef std::map<std::string, Pattern> PatternMap; // key is pattern name

//------------------------------------------------------------------------------------
// headingsMatch() returns true if two vectors of headings are identical

bool headingsMatch(const StringVector& headingA, const StringVector& headingB)
{
   int numHeadings = headingA.size();

   if (headingB.size() != numHeadings)
      return false;

   for (int i = 0; i < numHeadings; i++)
      if (headingA[i] != headingB[i])
         return false;

   return true;
}

//------------------------------------------------------------------------------------
// increment() updates the pattern map

void increment(PatternMap& pmap, const std::string& id, uint64_t readPairs,
               const std::string& patternName, const StringVector& annotation)
{
   PatternMap::iterator ppos = pmap.find(patternName);

   if (ppos == pmap.end())
   {
      pmap.insert(std::make_pair(patternName, Pattern(patternName, annotation)));
      ppos = pmap.find(patternName);
   }

   Pattern& pattern = ppos->second;
   IdMap&   idmap   = pattern.idmap;

   IdMap::iterator idpos = idmap.find(id);

   if (idpos == idmap.end())
      idmap.insert(std::make_pair(id, readPairs));
   else
      idpos->second += readPairs;
}

//------------------------------------------------------------------------------------
// readFuzzumFile() reads the named fuzzum file and updates the map; annotation
// headings, if any, are stored in the vector

void readFuzzumFile(const std::string& filename, PatternMap& pmap,
                    StringVector& annotationHeading)
{
   std::ifstream infile(filename.c_str());
   if (!infile.is_open())
      throw std::runtime_error("unable to open " + filename);

   std::string  headingLine, line;
   StringVector headingCol;

   if (!getline(infile, headingLine))
      throw std::runtime_error("empty file " + filename);

   int numCols = splitString(headingLine, headingCol);

   if (numCols < MIN_HEADING_COLS ||
       !hasPrefix(headingCol[0], FUZZUM) ||
       headingCol[RDPAIR_COL]  != READ_PAIRS ||
       headingCol[PATTERN_COL] != PATTERN)
      throw std::runtime_error("unexpected heading line in " + filename);

   for (int i = MIN_HEADING_COLS; i < numCols; i++)
      annotationHeading.push_back(headingCol[i]);

   while (getline(infile, line))
   {
      if (hasPrefix(line, FUZZUM)) // found another heading line
         if (line == headingLine)
            continue; // ignore it
         else
            throw std::runtime_error("inconsistent heading lines in " + filename);

      StringVector col;

      if (splitString(line, col) != numCols)
         throw std::runtime_error("inconsistent #columns in " + filename);

      std::istringstream stream(col[RDPAIR_COL]);
      uint64_t readPairs = 0;
      stream >> readPairs;

      if (readPairs == 0)
         throw std::runtime_error("invalid #read pairs in " + filename);

      StringVector annotation;

      for (int i = MIN_HEADING_COLS; i < numCols; i++)
         annotation.push_back(col[i]);

      increment(pmap, col[0], readPairs, col[PATTERN_COL], annotation);
   }

   infile.close();
}

//------------------------------------------------------------------------------------
// readAllFuzzumFiles() reads all of the fuzzum files and updates the map; annotation
// headings, if any, are stored in the vector

void readAllFuzzumFiles(PatternMap& pmap, StringVector& annotationHeading)
{
   int numFiles = fuzzumFilename.size();

   readFuzzumFile(fuzzumFilename[0], pmap, annotationHeading);

   for (int i = 1; i < numFiles; i++)
   {
      StringVector heading;

      readFuzzumFile(fuzzumFilename[i], pmap, heading);

      if (!headingsMatch(heading, annotationHeading))
         throw std::runtime_error("inconsistent heading lines in input files");
   }
}

//------------------------------------------------------------------------------------
// writeHeadingLine() writes a heading line to stdout

void writeHeadingLine(const StringVector& annotationHeading)
{
   std::cout << VERSION_NAME
             << TAB << "min"
             << TAB << "median"
             << TAB << "mean"
             << TAB << "max"
             << TAB << "total"
             << TAB << "IDs"
             << TAB << "ID list";

   int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotationHeading[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// Pattern::write() writes information about a pattern to stdout

void Pattern::write()
{
   std::vector<uint64_t> readPairs;
   uint64_t total = 0;

   for (IdMap::iterator idpos = idmap.begin(); idpos != idmap.end(); ++idpos)
   {
      readPairs.push_back(idpos->second);
      total += idpos->second;
   }

   std::sort(readPairs.begin(), readPairs.end());

   int n = readPairs.size();

   uint64_t min  = readPairs[0];
   uint64_t max  = readPairs[n - 1];

   double median = (n % 2 == 1 ? readPairs[n / 2] :
                                (readPairs[n / 2] + readPairs[n / 2 - 1]) / 2.0);

   double mean   = total / static_cast<double>(n);

   std::cout << name
             << TAB << min
	     << TAB << std::fixed << std::setprecision(1) << median
	     << TAB << std::fixed << std::setprecision(1) << mean
             << TAB << max
             << TAB << total
             << TAB << n
	     << TAB;

   bool first = true;

   for (IdMap::iterator idpos = idmap.begin(); idpos != idmap.end(); ++idpos)
   {
      if (first)
         first = false;
      else
         std::cout << ", ";

      std::cout << idpos->first << "(" << idpos->second << ")";
   }

   int numAnnotations = annotation.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotation[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// writePatterns() writes information about all patterns to stdout

void writePatterns(PatternMap& pmap)
{
   for (PatternMap::iterator ppos = pmap.begin(); ppos != pmap.end(); ++ppos)
   {
      Pattern& pattern = ppos->second;
      pattern.write();
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
      PatternMap   pmap;
      StringVector annotationHeading;

      readAllFuzzumFiles(pmap, annotationHeading);

      writeHeadingLine(annotationHeading);
      writePatterns(pmap);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
