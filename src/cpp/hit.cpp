//------------------------------------------------------------------------------------
//
// hit.cpp - module for reading and writing fuzzion2 hits
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "hit.h"

const String PATTERN   = "pattern ";
const String READ      = "read ";
const String RDCOUNT   = "#reads ";

// column headings
const String VIS       = "sequence";
const String MATCHES   = "matching bases";
const String POSSIBLE  = "possible";
const String PERCENT   = "% match";
const String SPANNING  = "spanning";
const String OVLEFT    = "left overlap";
const String OVRIGHT   = "right overlap";
const String ISIZE     = "isize|rdlen";
const String HASH      = "hash";

const int VIS_COL         =  1;
const int MATCHES_COL     =  2;
const int POSSIBLE_COL    =  3;
const int PERCENT_COL     =  4;
const int SPANNING_COL    =  5;
const int OVLEFT_COL      =  6;
const int OVRIGHT_COL     =  7;
const int ISIZE_COL       =  8;
const int HASH_COL        =  9;
const int FIRST_ANNOT_COL = 10; // first annotation column, if any

//------------------------------------------------------------------------------------
// computeHash() computes a 31-bit integer to represent the given sequence; it is used
// to detect duplicate read sequences

static int computeHash(const Seq& seq)
{
   const char *cstr = seq.cstr;
   const int   len  = seq.len;
   uint32_t hash = 0;

   for (int i = 0; i < len; i++)
   {
      const uint8_t ch = cstr[i];
      hash ^= (ch << (i % 24));
   }

   return hash;
}

//------------------------------------------------------------------------------------
// findDelimiter() returns the offset of the first delimiter in the pattern
// visualization

static int findDelimiter(const String& patternVis)
{
   const int len = patternVis.length();

   int i = 0;
   while (i < len && patternVis[i] != ']' && patternVis[i] != '}' &&
          patternVis[i] != '(')
      i++;

   return i;
}

//------------------------------------------------------------------------------------
// Hit::Hit() performs initialization; a 63-bit hash is computed to represent a read
// pair

Hit::Hit(const String& inPatternName, const String& inPatternVis,
         const int inPossible, const int inMatches, const int inSpanning,
	 const int inInsertSize, const StringVector& inAnnotation,
	 const HitRead *inRead1, const HitRead *inRead2)
   : patternName(inPatternName), patternVis(inPatternVis),
     delim1(findDelimiter(inPatternVis)), possible(inPossible), matches(inMatches),
     spanning(inSpanning), insertSize(inInsertSize), annotation(inAnnotation),
     read1(inRead1), read2(inRead2), duplicate(false)
{
   // compute the hash value representing the read(s)
   hash = (read2 ? read2->hash : 0);
   hash = (hash << 32) + read1->hash;
}

//------------------------------------------------------------------------------------
// createHitFromSingleMatch() returns a Hit to represent the given single-read match;
// it is the caller's obligation to de-allocate it

Hit *createHitFromSingleMatch(const Pattern& pattern, const int maxmidlen,
                              const double minPercentAgreement, const int minOverlap,
			      const String& readName, const Seq& readSeq,
			      const SingleMatch& match)
{
   String patternVis, readVis;
   match.getVis(pattern, maxmidlen, patternVis, readVis);

   // the alignment of a read is considered to be spanning if it spans the left and
   // right sides or the pattern is unsided
   int spanning = 0;
   if (!pattern.right || match.spanning(minPercentAgreement, minOverlap))
      spanning = 1;

   HitRead *read = new HitRead(readName, readVis, match.possible, match.matches,
                               spanning, match.loverlap, match.roverlap, readSeq.len,
			       computeHash(readSeq));

   return new Hit(pattern.name, patternVis, match.possible, match.matches, spanning,
                  0, pattern.annotation, read);
}

//------------------------------------------------------------------------------------
// createHitFromPairMatch() returns a Hit to represent the given read-pair match; it
// is the caller's obligation to de-allocate it

Hit *createHitFromPairMatch(const Pattern& pattern, const int maxmidlen,
                            const double minPercentAgreement, const int minOverlap,
			    const String& readName1, const Seq& readSeq1,
			    const String& readName2, const Seq& readSeq2,
			    const PairMatch& pairMatch)
{
   StringVector vis;
   pairMatch.getVis(pattern, maxmidlen, vis);

   // the alignment of a read is considered to be spanning if it spans the left and
   // right sides or the pattern is unsided
   int spanning1 = 0;
   if (!pattern.right || pairMatch.match1->spanning(minPercentAgreement, minOverlap))
      spanning1 = 1;

   int spanning2 = 0;
   if (!pattern.right || pairMatch.match2->spanning(minPercentAgreement, minOverlap))
      spanning2 = 1;

   HitRead *read1 = new HitRead(readName1, vis[1], pairMatch.match1->possible,
                                pairMatch.match1->matches, spanning1,
				pairMatch.match1->loverlap,
				pairMatch.match1->roverlap, readSeq1.len,
				computeHash(readSeq1));

   HitRead *read2 = new HitRead(readName2, vis[2], pairMatch.match2->possible,
                                pairMatch.match2->matches, spanning2,
				pairMatch.match2->loverlap,
				pairMatch.match2->roverlap, readSeq2.len,
				computeHash(readSeq2));

   return new Hit(pattern.name, vis[0], pairMatch.possible, pairMatch.matches,
                  spanning1 + spanning2, pairMatch.insertSize, pattern.annotation,
		  read1, read2);
}

//------------------------------------------------------------------------------------
// Hit::sameAs() returns true if this hit and the other hit are duplicates

bool Hit::sameAs(const Hit& other) const
{
   return (patternName == other.patternName && delim1 == other.delim1 &&
           hash == other.hash);
}

//------------------------------------------------------------------------------------
// HitCompare defines the sort order of hits; duplicate hits are placed consecutively
// in the ordering

struct HitCompare
{
   bool operator()(Hit* const& a, Hit* const& b) const
   {
      // sort by ascending pattern name,
      // then by ascending offset of first pattern delimiter,
      // then by ascending hash,
      // then by ascending read1 name

      const int key1 = std::strcmp(a->patternName.c_str(), b->patternName.c_str());
      if (key1 != 0)
         return (key1 < 0);

      const int key2 = a->delim1 - b->delim1;
      if (key2 != 0)
         return (key2 < 0);

      if (a->hash != b->hash)
         return (a->hash < b->hash);

      return (a->read1->name < b->read1->name);
   }
};

//------------------------------------------------------------------------------------
// sortHits() sorts the hits in the given vector

static void sortHits(HitVector& hitVector)
{
   std::sort(hitVector.begin(), hitVector.end(), HitCompare());
}

//------------------------------------------------------------------------------------
// markDuplicates() marks duplicate hits in the given vector, which is assumed to be
// sorted

static void markDuplicates(HitVector& hitVector)
{
   const int numHits = hitVector.size();

   for (int i = 1; i < numHits; i++)
      if (hitVector[i]->sameAs(*hitVector[i - 1]))
         hitVector[i]->duplicate = true;
}

//------------------------------------------------------------------------------------
// Hit::isStrong() returns true if this hit is considered a "strong" hit

bool Hit::isStrong(const int minStrong) const
{
   const int lmax =
      (read2 ? std::max(read1->loverlap, read2->loverlap) : read1->loverlap);
   const int rmax =
      (read2 ? std::max(read1->roverlap, read2->roverlap) : read1->roverlap);

   // lmax will be zero if the pattern does not have a left side;
   // rmax will be zero if the pattern does not have a right side

   return ((lmax == 0 || lmax >= minStrong) &&
           (rmax == 0 || rmax >= minStrong));
}

//------------------------------------------------------------------------------------
// Hit::label() returns the hit category

String Hit::label(const int minStrong) const
{
   if (duplicate)
      return HIT_DUPLICATE;

   if (isStrong(minStrong))
      if (spanning > 0)
         return HIT_STRONG_SPAN;
      else
         return HIT_STRONG_NOSPAN;

   return HIT_WEAK;
}

//------------------------------------------------------------------------------------
// writeHitHeadingLine() writes a heading line

void writeHitHeadingLine(std::ostream& ostream, const String& version,
                         const StringVector& annotationHeading)
{
   ostream << FUZZION2 << version
           << TAB << VIS
           << TAB << MATCHES
	   << TAB << POSSIBLE
	   << TAB << PERCENT
	   << TAB << SPANNING
	   << TAB << OVLEFT
	   << TAB << OVRIGHT
	   << TAB << ISIZE
	   << TAB << HASH;

   const int numAnnotations = annotationHeading.size();

   for (int i = 0; i < numAnnotations; i++)
      ostream << TAB << annotationHeading[i];

   ostream << NEWLINE;
}

//------------------------------------------------------------------------------------
// isHeadingLine() returns true if the given line is a heading line

static bool isHeadingLine(const String& line)
{
   return hasPrefix(line, FUZZION2);
}

//------------------------------------------------------------------------------------
// validHeadingLine() returns true if the given line is a valid heading line

static bool validHeadingLine(const String& line, String& version,
                             StringVector& annotationHeading)
{
   if (!isHeadingLine(line))
      return false;

   StringVector col;
   const int numCols = splitString(line, col);

   if (numCols <= HASH_COL ||
       col[VIS_COL]      != VIS      || col[MATCHES_COL] != MATCHES ||
       col[POSSIBLE_COL] != POSSIBLE || col[PERCENT_COL] != PERCENT ||
       col[SPANNING_COL] != SPANNING || col[OVLEFT_COL]  != OVLEFT  ||
       col[OVRIGHT_COL]  != OVRIGHT  || col[ISIZE_COL]   != ISIZE   ||
       col[HASH_COL]     != HASH)
      return false;

   StringVector versionCol;
   if (splitString(col[0], versionCol, ' ') != 2 || versionCol[1] == "")
      return false;

   version = versionCol[1];

   annotationHeading.clear();

   for (int i = FIRST_ANNOT_COL; i < numCols; i++)
      annotationHeading.push_back(col[i]);

   return true;
}

//------------------------------------------------------------------------------------
// writeReadCountLine() writes one line showing the total number of reads processed

void writeReadCountLine(std::ostream& ostream, const uint64_t numReads)
{
   ostream << RDCOUNT << numReads << NEWLINE;
}

//------------------------------------------------------------------------------------
// isReadCountLine() returns true if the given line shows the total number of reads
// processed

static bool isReadCountLine(const String& line)
{
   return hasPrefix(line, RDCOUNT);
}

//------------------------------------------------------------------------------------
// validReadCountLine() returns true if the given line is a valid read count line

static bool validReadCountLine(const String& line, uint64_t& numReads)
{
   if (!isReadCountLine(line))
      return false;

   StringVector col;
   if (splitString(line, col, ' ') != 2 || col[1] == "")
      return false;

   numReads = UINT64_MAX;
   std::istringstream stream(col[1]);
   stream >> numReads;

   return (numReads != UINT64_MAX);
}

//------------------------------------------------------------------------------------
// HitRead::write() writes one line describing a read that matches the pattern

void HitRead::write(std::ostream& ostream) const
{
   ostream << READ << name
           << TAB << vis
	   << TAB << matches
	   << TAB << possible
	   << TAB << doubleToString(100.0 * matches / possible)
	   << TAB << spanning
	   << TAB << loverlap
	   << TAB << roverlap
	   << TAB << len
	   << TAB << hash
	   << NEWLINE;
}

//------------------------------------------------------------------------------------
// Hit::write() writes a hit, first a line describing the pattern, and then one or two
// lines describing matching reads

void Hit::write(std::ostream& ostream) const
{
   ostream << PATTERN << patternName
           << TAB << patternVis
	   << TAB << matches
	   << TAB << possible
	   << TAB << doubleToString(100.0 * matches / possible)
	   << TAB << spanning
	   << TAB // no entry in this column
	   << TAB // no entry in this column
	   << TAB << insertSize
	   << TAB;// no entry in this column

   const int numAnnotations = annotation.size();

   for (int i = 0; i < numAnnotations; i++)
      ostream << TAB << annotation[i];

   ostream << NEWLINE;

   read1->write(ostream);

   if (read2)
      read2->write(ostream);
}

//------------------------------------------------------------------------------------
// isReadLine() returns true if the given line describes a read that matched a pattern

static bool isReadLine(const String& line)
{
   return hasPrefix(line, READ);
}

//------------------------------------------------------------------------------------
// getRead() returns a pointer to a newly-allocated HitRead describing the read
// identified by the given line, or returns nullptr if the input is invalid

static HitRead *getRead(const String& line)
{
   if (!isReadLine(line))
      return nullptr;

   StringVector col;
   const int numCols = splitString(line, col);

   int matches, possible, spanning, loverlap, roverlap, len, hash;

   if (numCols <= HASH_COL ||
       (matches  = stringToInt(col[MATCHES_COL]))  <  0 ||
       (possible = stringToInt(col[POSSIBLE_COL])) <= 0 ||
       (spanning = stringToInt(col[SPANNING_COL])) <  0 || spanning > 1 ||
       (loverlap = stringToInt(col[OVLEFT_COL]))   <  0 ||
       (roverlap = stringToInt(col[OVRIGHT_COL]))  <  0 ||
       (len      = stringToInt(col[ISIZE_COL]))    <= 0 ||
       (hash     = stringToInt(col[HASH_COL]))     <  0)
      return nullptr;

   StringVector nameCol;
   if (splitString(col[0], nameCol, ' ') != 2 || nameCol[1] == "")
      return nullptr;

   const String& name = nameCol[1];
   const String& vis  = col[VIS_COL];

   return new HitRead(name, vis, possible, matches, spanning, loverlap, roverlap, len,
                      hash);
}

//------------------------------------------------------------------------------------
// isPatternLine() returns true if the given line describes a pattern that was matched

static bool isPatternLine(const String& line)
{
   return hasPrefix(line, PATTERN);
}

//------------------------------------------------------------------------------------
// getHit() returns a pointer to a newly-allocated Hit describing the pattern
// identified by the given line, or returns nullptr if the input is invalid

static Hit *getHit(const String& line, const HitRead *read1,
                   const HitRead *read2=nullptr)
{
   if (!isPatternLine(line))
      return nullptr;

   StringVector col;
   const int numCols  = splitString(line, col);

   int matches, possible, spanning, insertSize;

   if (numCols <= HASH_COL ||
       (matches    = stringToInt(col[MATCHES_COL]))  <= 0 ||
       (possible   = stringToInt(col[POSSIBLE_COL])) <= 0 ||
       (spanning   = stringToInt(col[SPANNING_COL])) <  0 || spanning > 2 ||
       (insertSize = stringToInt(col[ISIZE_COL]))    <  0)
      return nullptr;

   StringVector nameCol;
   if (splitString(col[0], nameCol, ' ') != 2 || nameCol[1] == "")
      return nullptr;

   const String& name = nameCol[1];
   const String& vis  = col[VIS_COL];

   StringVector annotation;

   for (int i = FIRST_ANNOT_COL; i < numCols; i++)
      annotation.push_back(col[i]);

   return new Hit(name, vis, possible, matches, spanning, insertSize, annotation,
                  read1, read2);
}

//------------------------------------------------------------------------------------
// readHits() reads hits and stores them in sorted order in hitVector (it is the
// caller's obligation to de-allocate them); duplicate hits are marked

void readHits(std::istream& istream, String& version, StringVector& annotationHeading,
	      HitVector& hitVector, uint64_t& numReads)
{
   numReads = 0;
   String headingLine, line;

   if (!getline(istream, headingLine))
      throw std::runtime_error("no input");

   if (!validHeadingLine(headingLine, version, annotationHeading))
      throw std::runtime_error("unexpected heading line");

   if (!getline(istream, line))
      return; // no hits

   while (true)
      if (isHeadingLine(line)) // found another heading line
      {
         if (line != headingLine)
            throw std::runtime_error("inconsistent heading lines");

	 if (!getline(istream, line))
            break;
      }
      else if (isReadCountLine(line))
      {
         uint64_t readCount;

	 if (validReadCountLine(line, readCount))
            numReads += readCount;
	 else
            throw std::runtime_error("unexpected input line: " + line);

	 if (!getline(istream, line))
            break;
      }
      else // this must be a hit expressed on two or three lines
      {
         String rdline1, rdline2;
	 HitRead *read1, *read2;
	 Hit *hit = nullptr;

	 if (isPatternLine(line)        && getline(istream, rdline1) &&
             (read1 = getRead(rdline1)) && getline(istream, rdline2))
            if (isReadLine(rdline2)) // found a second read
	    {
               if (read2 = getRead(rdline2))
                  hit = getHit(line, read1, read2);
	    }
	    else // only one read
               hit = getHit(line, read1);

	 if (hit)
            if (hitVector.size() < MAX_HITS)
	    {
               hitVector.push_back(hit);

	       if (!hit->read2)
                  line = rdline2;
	       else if (!getline(istream, line))
                  break;
	    }
	    else
               throw std::runtime_error("too many hits");
	 else
            throw std::runtime_error("unexpected hit format: " + line);
      }

   if (hitVector.size() > 1)
   {
      sortHits(hitVector);
      markDuplicates(hitVector);
   }
}

//------------------------------------------------------------------------------------
// getPatternIndices() determines the index of the first hit of each pattern; it is
// assumed that hitVector has been sorted

void getPatternIndices(const HitVector& hitVector, IntVector& index)
{
   index.clear();

   const int numHits = hitVector.size();
   if (numHits == 0)
      return;

   index.push_back(0);

   for (int i = 1; i < numHits; i++)
      if (hitVector[i]->patternName != hitVector[i - 1]->patternName)
         index.push_back(i);
}

//------------------------------------------------------------------------------------
// maxVisLength() returns the maximum length of the visualization sequence for
// the hits in the given vector, from indices begin (inclusive) to end (exclusive)

int maxVisLength(const HitVector& hitVector, const int begin, const int end)
{
   int maxlen = 0;

   for (int i = begin; i < end; i++)
   {
      const int len = hitVector[i]->patternVis.length();
      if (len > maxlen)
         maxlen = len;
   }

   return maxlen;
}
