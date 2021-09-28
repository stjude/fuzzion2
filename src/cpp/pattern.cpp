//------------------------------------------------------------------------------------
//
// pattern.cpp - module supporting pattern sequences
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "pattern.h"
#include "window.h"
#include <fstream>
#include <stdexcept>

const std::string PATTERN_HEADING  = "pattern";
const std::string SEQUENCE_HEADING = "sequence";

//------------------------------------------------------------------------------------
// Pattern::Pattern() constructs an object representing the named pattern

Pattern::Pattern(const std::string& inName, const std::string& inSequence,
                 const StringVector& inAnnotation)
   : name(inName), displaySequence(inSequence), annotation(inAnnotation)
{
   if (name.length() == 0)
      throw std::runtime_error("zero-length pattern name");

   if (name.find(' ') != std::string::npos)
      throw std::runtime_error("space not allowed in pattern name \"" + name + "\"");

   bool foundDelimiters = false;
   int  delim1, delim2;

   if ((delim1 = inSequence.find(']')) != std::string::npos &&
       (delim2 = inSequence.find('[')) != std::string::npos)
   {
      foundDelimiters = true;
      hasBraces       = false;
   }
   else if ((delim1 = inSequence.find('}')) != std::string::npos &&
            (delim2 = inSequence.find('{')) != std::string::npos)
   {
      foundDelimiters = true;
      hasBraces       = true;
   }

   if (!foundDelimiters || delim1 == 0 || delim2 == inSequence.length() - 1 ||
       delim1 > delim2)
      throw std::runtime_error("invalid pattern " + inSequence);

   int middleBases = delim2 - delim1 - 1;

   leftBases  = delim1;
   rightBases = inSequence.length() - 1 - delim2;

   sequence   = inSequence.substr(0, leftBases) +
                inSequence.substr(delim1 + 1, middleBases) +
		inSequence.substr(delim2 + 1, rightBases);
}

//------------------------------------------------------------------------------------
// readPatterns() returns a vector of patterns read from the named file; it is the
// caller's obligation to de-allocate it

PatternVector *readPatterns(const std::string& filename,
                            StringVector& annotationHeading)
{
   PatternVector *patternVector = new PatternVector();

   std::ifstream infile(filename.c_str());
   if (!infile.is_open())
      throw std::runtime_error("unable to open " + filename);

   std::string line;

   if (!getline(infile, line))
      throw std::runtime_error("empty file " + filename);

   StringVector heading;

   int numColumns = splitString(line, heading);

   if (numColumns < 2 || heading[0] != PATTERN_HEADING ||
       heading[1] != SEQUENCE_HEADING)
      throw std::runtime_error("invalid format in " + filename);

   for (int i = 2; i < numColumns; i++)
      annotationHeading.push_back(heading[i]);

   while (getline(infile, line))
   {
      StringVector column, annotation;

      if (splitString(line, column) != numColumns)
         throw std::runtime_error("inconsistent #columns in " + filename);

      for (int i = 2; i < numColumns; i++)
         annotation.push_back(column[i]);

      patternVector->push_back(Pattern(column[0], column[1], annotation));
   }

   infile.close();

   return patternVector;
}

//------------------------------------------------------------------------------------
// createPatternMap() returns a map constructed from the given pattern vector; it is
// the caller's obligation to de-allocate it

PatternMap *createPatternMap(const PatternVector *patternVector,
                             MinimizerWindowLength w,
			     const KmerRankTable *rankTable,
			     Minimizer maxMinimizer)
{
   PatternMap *patternMap = new PatternMap();

   int numPatterns = patternVector->size();

   for (int i = 0; i < numPatterns; i++)
   {
      WindowVector windowVector;

      getWindows((*patternVector)[i].sequence, w, rankTable, windowVector);

      int numWindows = windowVector.size();

      for (int j = 0; j < numWindows; j++)
      {
         Minimizer minimizer = windowVector[j].minimizer;
	 if (minimizer > maxMinimizer)
            continue; // don't put common minimizer in map

	 PatternMap::iterator ppos = patternMap->find(minimizer);

	 if (ppos == patternMap->end())
	 {
            patternMap->insert(std::make_pair(minimizer, LocationVector()));
	    ppos = patternMap->find(minimizer);
	 }

	 LocationVector& locationVector = ppos->second;

	 locationVector.push_back(Location(i, windowVector[j].offset));
      }
   }

   return patternMap;
}
