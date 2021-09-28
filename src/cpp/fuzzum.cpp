//------------------------------------------------------------------------------------
//
// fuzzum.cpp - this programs reads fuzzion2 hits from stdin and writes a summary to
//              stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "util.h"
#include <iostream>
#include <map>
#include <stdexcept>

const std::string VERSION_NAME = "fuzzum v1.0.2";
const std::string COPYRIGHT    =
   "copyright 2021 St. Jude Children's Research Hospital";

const std::string FUZZION2     = "fuzzion2 ";
const std::string PATTERN      = "pattern ";
const std::string READ         = "read ";
const std::string READ_PAIRS   = "read-pairs ";

const int MIN_HEADING_COLS     = 6;

std::string id = ""; // identifies the sample

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname << " OPTION ... < fuzzion2_hits > hit_summary"
      << NEWLINE;

   std::cerr
      << NEWLINE
      << "This option is required:" << NEWLINE
      << "  -id=string   "
             << "identifies the sample" << NEWLINE;
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
         return false; // not an option

      StringVector opt;

      if (splitString(arg, opt, '=') != 2)
         return false; // incorrect option format

      if (stringOpt(opt, "id", id))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   return (id != "");
}

//------------------------------------------------------------------------------------

class Pattern
{
public:
   Pattern(const std::string& inName, const StringVector& inAnnotation,
           uint64_t inReadPairs)
      : name(inName), annotation(inAnnotation), readPairs(inReadPairs) { }

   virtual ~Pattern() { }

   void write() const;

   std::string  name;
   StringVector annotation;
   uint64_t readPairs;
};

typedef std::map<std::string, Pattern> PatternMap; // key is pattern name

//------------------------------------------------------------------------------------
// readHits() reads the hits from stdin and stores them in the map; annotation
// headings, if any, are stored in the vector

void readHits(PatternMap& pmap, StringVector& annotationHeading)
{
   std::string headingLine, line1, line2, line3;
   StringVector headingCol;

   if (!getline(std::cin, headingLine))
      throw std::runtime_error("no input");

   int numCols = splitString(headingLine, headingCol);

   if (numCols < MIN_HEADING_COLS || !hasPrefix(headingCol[0], FUZZION2))
      throw std::runtime_error("unexpected heading line");

   for (int i = MIN_HEADING_COLS; i < numCols; i++)
      annotationHeading.push_back(headingCol[i]);

   while (getline(std::cin, line1))
   {
      if (hasPrefix(line1, FUZZION2)) // found another heading line
         if (line1 == headingLine)
            continue; // ignore it
         else
            throw std::runtime_error("inconsistent heading lines");

      if (hasPrefix(line1, READ_PAIRS)) // found total number of read pairs
         continue; // ignore it

      StringVector col, part;

      if (!hasPrefix(line1, PATTERN) ||
          splitString(line1, col) != numCols ||
          splitString(col[0], part, ' ') != 2 || part[1] == "" ||
          !getline(std::cin, line2) || !hasPrefix(line2, READ) ||
          !getline(std::cin, line3) || !hasPrefix(line3, READ))
         throw std::runtime_error("invalid format in: " + line1);

      const std::string& patternName = part[1];

      PatternMap::iterator ppos = pmap.find(patternName);

      if (ppos == pmap.end())
      {
         StringVector annotation;

	 for (int i = MIN_HEADING_COLS; i < numCols; i++)
            annotation.push_back(col[i]);

	 pmap.insert(std::make_pair(patternName,
                                    Pattern(patternName, annotation, 1)));
      }
      else
      {
         Pattern& pattern = ppos->second;
	 pattern.readPairs++;
      }
   }
}

//------------------------------------------------------------------------------------
// writeHeadingLine() writes a heading line to stdout

void writeHeadingLine(const StringVector& annotationHeading)
{
   std::cout << VERSION_NAME
             << TAB << "read pairs"
	     << TAB << "pattern";

   int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotationHeading[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// Pattern::write() writes information about a pattern to stdout

void Pattern::write() const
{
   std::cout << id
             << TAB << readPairs
	     << TAB << name;

   int numAnnotations = annotation.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << annotation[i];

   std::cout << NEWLINE;
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

      readHits(pmap, annotationHeading);

      writeHeadingLine(annotationHeading);

      for (PatternMap::iterator ppos = pmap.begin(); ppos != pmap.end(); ++ppos)
      {
         const Pattern& pattern = ppos->second;
	 pattern.write();
      }
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
