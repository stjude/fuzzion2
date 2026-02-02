//------------------------------------------------------------------------------------
//
// pattern.cpp - module supporting patterns
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "pattern.h"
#include "window.h"
#include <cstring>
#include <fstream>
#include <set>
#include <stdexcept>

const String PATTERN_HEADING  = "pattern";
const String SEQUENCE_HEADING = "sequence";

//------------------------------------------------------------------------------------

class TrimerFinder : public KmerFinder // for finding 3-mers in a sequence
{
public:
   TrimerFinder(const char *cstr, const int begin, const int len,
                TrimerOffsets& trimerOffsets)
      : KmerFinder(&cstr[begin], len, TRIMER_LEN), cstrOffset(begin),
	tOffsets(trimerOffsets) { }

   virtual ~TrimerFinder() { }

   virtual bool reportKmer(const Kmer trimer, const int startIndex) override
   {
      // save the offset of this 3-mer
      tOffsets[trimer].push_back(cstrOffset + startIndex);
      return true;
   }

   const int cstrOffset;    // starting offset of sequence
   TrimerOffsets& tOffsets; // records the offsets of each 3-mer
};

//------------------------------------------------------------------------------------
// initTrimerOffsets() finds 3-mers in the given sequence and saves their offsets

void initTrimerOffsets(const char *cstr, const int begin, const int len,
                       TrimerOffsets& trimerOffsets)
{
   TrimerFinder finder(cstr, begin, len, trimerOffsets);
   finder.find();
}

//------------------------------------------------------------------------------------
// Seq::Seq() allocates and initializes a C-style string, converting lowercase acgt
// to uppercase ACGT; invalidBase is set to true if the string contains an
// unrecognized base

Seq::Seq(const String& inStr, bool& invalidBase)
   : str(inStr), len(inStr.length())
{
   cstr = new char[len + 1];
   invalidBase = false;

   for (int i = 0; i < len; i++)
      switch (str[i])
      {
         case 'A': case 'a': cstr[i] = 'A'; break;
         case 'C': case 'c': cstr[i] = 'C'; break;
         case 'G': case 'g': cstr[i] = 'G'; break;
         case 'T': case 't': cstr[i] = 'T'; break;
         default : cstr[i] = str[i]; invalidBase = true;
      }

   cstr[len] = 0;
}

//------------------------------------------------------------------------------------
// Seq::Seq() copy constructor

Seq::Seq(const Seq& other)
   : str(other.str), len(other.len)
{
   cstr = new char[len + 1];
   std::strcpy(cstr, other.cstr);
}

//------------------------------------------------------------------------------------
// findChar() returns the offset of the first occurrence of a character in a string,
// or returns -1 if the character is not found

static int findChar(const String& s, const char c)
{
   const std::size_t offset = s.find(c);
   return (offset == String::npos ? -1 : offset);
}

//------------------------------------------------------------------------------------
// validateNoDelims() validates a sequence containing no delimiters

static bool validateNoDelims(const String& sequence, Seq *& left)
{
   bool invalidBase;
   left = new Seq(sequence, invalidBase);

   if (!invalidBase)
      return true;

   delete left; left = nullptr;
   return false;
}

//------------------------------------------------------------------------------------
// validateJunction() validates a sequence describing a junction

static bool validateJunction(const String& sequence,
                             const int delim1, const int delim2,
                             Seq *& left, Seq *& middle, Seq *& right)
{
   if (delim1 == 0 || delim1 > delim2 || delim2 == sequence.length() - 1)
      return false;

   const String lsequence = sequence.substr(0, delim1);
   const String msequence = sequence.substr(delim1 + 1, delim2 - delim1 - 1);
   const String rsequence = sequence.substr(delim2 + 1);

   bool invalidLeftBase, invalidMiddleBase, invalidRightBase;
   left   = new Seq(lsequence, invalidLeftBase);
   middle = new Seq(msequence, invalidMiddleBase);
   right  = new Seq(rsequence, invalidRightBase);

   if (!invalidLeftBase && !invalidRightBase)
      if (msequence == "*") // wildcard middle
      {
         delete middle; middle = nullptr;
	 return true;
      }
      else if (!invalidMiddleBase)
         return true;

   delete left;   left   = nullptr;
   delete middle; middle = nullptr;
   delete right;  right  = nullptr;
   return false;
}

//------------------------------------------------------------------------------------
// validateExclusion() validates a sequence containing exclusions

static bool validateExclusion(const String& sequence,
                              const int lparen, const int rparen,
                              Seq *& left, Seq *& right, SeqVector *& xvector,
			      String& xvis)
{
   if (lparen == 0 || lparen >= rparen - 1 || rparen == sequence.length() - 1)
      return false;

   // left(ex)right
   const String lsequence = sequence.substr(0, lparen);
   const String xsequence = sequence.substr(lparen + 1, rparen - lparen - 1);
   const String rsequence = sequence.substr(rparen + 1);

   bool invalidLeftBase, invalidRightBase;
   left  = new Seq(lsequence, invalidLeftBase);
   right = new Seq(rsequence, invalidRightBase);

   if (!invalidLeftBase && !invalidRightBase)
   {
      StringVector v;
      const int n = splitString(xsequence, v, '|'); // exclusions separated by '|'

      xvector = new SeqVector();
      bool invalid = false;

      for (int i = 0; i < n && !invalid; i++)
         if (v[i].length() == 0)
            invalid = true;
         else
            xvector->emplace_back(v[i], invalid);

      if (!invalid)
      {
         xvis = xsequence;
         return true;
      }

      delete xvector; xvector = nullptr;
   }

   delete left;  left  = nullptr;
   delete right; right = nullptr;
   return false;
}

//------------------------------------------------------------------------------------
// Pattern::Pattern() constructs an object representing the named pattern

Pattern::Pattern(const String& inName, const String& inSequence,
                 const StringVector& inAnnotation)
   : name(inName), sequence(inSequence), annotation(inAnnotation),
     delim1(-1), delim2(-1), left(nullptr), middle(nullptr), right(nullptr),
     xvector(nullptr), xvis(""), leftTrimers(), rightTrimers()
{
   if (name.length() == 0)
      throw std::runtime_error("missing pattern name");

   if (findChar(name, ' ') >= 0)
      throw std::runtime_error("space not allowed in pattern name \"" + name + "\"");

   if (sequence.length() == 0)
      throw std::runtime_error("missing sequence for pattern \"" + name + "\"");

   if (findChar(sequence, ' ') >= 0)
      throw std::runtime_error("space not allowed in sequence for pattern \"" + name +
                               "\"");

   const int lbracket = findChar(sequence, '['), rbracket = findChar(sequence, ']');
   const int lbrace   = findChar(sequence, '{'), rbrace   = findChar(sequence, '}');
   const int lparen   = findChar(sequence, '('), rparen   = findChar(sequence, ')');

   if (lbracket < 0 && rbracket < 0 && lbrace < 0 && rbrace < 0 &&
       lparen < 0 && rparen < 0) // no delimiters in sequence
   {
      if (validateNoDelims(sequence, left))
         return;
   }
   else if (lbracket >= 0 && rbracket >= 0)
   {
      delim1 = rbracket; delim2 = lbracket;

      if (validateJunction(sequence, rbracket, lbracket, left, middle, right))
         return;
   }
   else if (lbrace >= 0 && rbrace >= 0)
   {
      delim1 = rbrace; delim2 = lbrace;

      if (validateJunction(sequence, rbrace, lbrace, left, middle, right))
         return;
   }
   else if (lparen >= 0 && rparen >= 0)
   {
      delim1 = lparen; delim2 = rparen;

      if (validateExclusion(sequence, lparen, rparen, left, right, xvector, xvis))
         return;
   }

   throw std::runtime_error("invalid sequence for pattern \"" + name + "\"");
}

//------------------------------------------------------------------------------------
// readPatterns() returns a vector of patterns read from the named file; it is the
// caller's obligation to de-allocate it

PatternVector *readPatterns(const String& filename, StringVector& annotationHeading)
{
   PatternVector *pvector = new PatternVector();

   std::ifstream infile(filename.c_str());
   if (!infile.is_open())
      throw std::runtime_error("unable to open " + filename);

   String line;

   if (!getline(infile, line))
      throw std::runtime_error("empty file " + filename);

   StringVector heading;
   const int numCols = splitString(line, heading);

   if (numCols < 2 || heading[0] != PATTERN_HEADING || heading[1] != SEQUENCE_HEADING)
      throw std::runtime_error("invalid format in " + filename);

   for (int i = 2; i < numCols; i++)
      annotationHeading.push_back(heading[i]);

   std::set<String> nameSet; // set of pattern names

   while (getline(infile, line))
   {
      StringVector col, annotation;

      if (splitString(line, col) != numCols)
         throw std::runtime_error("inconsistent #columns in " + filename);

      const String& name     = col[0];
      const String& sequence = col[1];

      if (nameSet.find(name) != nameSet.end())
         throw std::runtime_error("duplicate pattern name \"" + name + "\"");

      nameSet.insert(name);

      for (int i = 2; i < numCols; i++)
         annotation.push_back(col[i]);

      pvector->push_back(new Pattern(name, sequence, annotation));
   }

   infile.close();

   if (pvector->size() == 0)
      throw std::runtime_error("no patterns specified in " + filename);

   return pvector;
}

//------------------------------------------------------------------------------------
// getMinimizers() gets minimizers from a pattern sequence and stores them in a map;
// the number of minimizers stored in the map is returned

static int getMinimizers(const Seq *seq, const int side, const int pindex,
                         const MinimizerWindowLength w,
			 const KmerRankTable *rankTable,
			 const Minimizer maxMinimizer, PatternMap& pmap,
			 BoolVector& inPatternMap)
{
   WindowVector wvector;
   getWindows(seq->cstr, seq->len, w, rankTable, wvector);
   const int numWindows = wvector.size();
   int numMinimizers = 0;

   for (int i = 0; i < numWindows; i++)
   {
      const Minimizer minimizer = wvector[i].minimizer;
      if (minimizer > maxMinimizer)
         continue; // don't put common minimizer in map

      inPatternMap[minimizer] = true;

      PatternMap::iterator ppos = pmap.find(minimizer);
      if (ppos == pmap.end())
         ppos = pmap.insert(std::make_pair(minimizer, PlocDuo())).first;

      PlocDuo& duo = ppos->second;
      duo[side].emplace_back(pindex, wvector[i].offset);

      numMinimizers++;
   }

   return numMinimizers;
}

//------------------------------------------------------------------------------------
// createPatternMap() returns a map constructed from the given pattern vector and
// creates a Boolean vector indicating which minimizers are in the map; it is the
// caller's obligation to de-allocate these

PatternMap *createPatternMap(PatternVector& pvector, const MinimizerWindowLength w,
			     const KmerRankTable *rankTable,
			     const Minimizer maxMinimizer, BoolVector *& inPatternMap)
{
   PatternMap *pmap = new PatternMap();
   inPatternMap     = new BoolVector(maxMinimizer + 1, false);

   const int numPatterns = pvector.size();

   for (int pindex = 0; pindex < numPatterns; pindex++)
   {
      Pattern *pattern = pvector[pindex];

      int numMinimizers = getMinimizers(pattern->left,  0, pindex, w, rankTable,
                                        maxMinimizer, *pmap, *inPatternMap);
      if (pattern->right)
         numMinimizers += getMinimizers(pattern->right, 1, pindex, w, rankTable,
                                        maxMinimizer, *pmap, *inPatternMap);

      if (numMinimizers == 0)
         throw std::runtime_error("the sequence for pattern \"" + pattern->name +
                               "\" has no minimizers; it is too short or too common");

      if (pattern->right)
      {
         initTrimerOffsets(pattern->left->cstr,  0, pattern->left->len,
                           pattern->leftTrimers);
	 initTrimerOffsets(pattern->right->cstr, 0, pattern->right->len,
                           pattern->rightTrimers);
      }
   }

   return pmap;
}
