//------------------------------------------------------------------------------------
//
// align.cpp - module for aligning substrings of sequences
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "align.h"
#include <algorithm>

//------------------------------------------------------------------------------------
// OriginCounter::incrementCount() increments the number of occurrences of the given
// origin

void OriginCounter::incrementCount(const Origin origin)
{
   OriginCountMap::iterator ocpos = ocmap.find(origin);

   if (ocpos == ocmap.end())
      ocmap[origin] = 1;
   else
      ocpos->second++;
}

//------------------------------------------------------------------------------------
// OriginCounter::maxCount() returns the maximum number of occurrences of an origin

int OriginCounter::maxCount() const
{
   int max = 0;

   for (OriginCountMap::const_iterator ocpos = ocmap.cbegin();
        ocpos != ocmap.cend(); ++ocpos)
      if (ocpos->second > max)
         max = ocpos->second;

   return max;
}

//------------------------------------------------------------------------------------
// OriginCounter::retrieveOrigins() puts all origins having at least minCount
// occurrences in the given vector, sorts the vector by ascending origin, and returns
// the number of these origins

int OriginCounter::retrieveOrigins(const int minCount, OriginVector& ovector) const
{
   for (OriginCountMap::const_iterator ocpos = ocmap.cbegin();
        ocpos != ocmap.cend(); ++ocpos)
      if (ocpos->second >= minCount)
         ovector.push_back(ocpos->first);

   int numOrigins = ovector.size();
   if (numOrigins > 1)
      std::sort(ovector.begin(), ovector.end());

   return numOrigins;
}

//------------------------------------------------------------------------------------
// Sub::Sub() constructs a substring that is before or after the given substring

Sub::Sub(const Sub& other, const bool before)
   : seq(other.seq)
{
   if (before)
   {
      begin = 0;
      len   = other.begin;
   }
   else // after
   {
      begin = other.end();
      len   = seq.len - begin;
   }
}

//------------------------------------------------------------------------------------
// initialSpacers() returns the number of initial spacer characters in the given
// visualization string

static int initialSpacers(const String& vis)
{
   const int len = vis.length();
   int count = 0;

   for (int i = 0; i < len && vis[i] == SPACER; i++)
      count++;

   return count;
}

//------------------------------------------------------------------------------------
// trailingSpacers() returns the number of trailing spacer characters in the given
// visualization string

static int trailingSpacers(const String& vis)
{
   const int len = vis.length();
   int count = 0;

   for (int i = len - 1; i >= 0 && vis[i] == SPACER; i--)
      count++;

   return count;
}

//------------------------------------------------------------------------------------
// adjustInitial() replaces each initial spacer in visA with a matching base if there
// is one, and excess initial bases are trimmed

static void adjustInitial(const int initial, const int lbases, Sub& subA, Sub& subB,
                          String& visA, String& visB)
{
   int i = initial - 1, bases = 0;
   while (i >= 0 && bases < lbases)
   {
      if (subA.seq.cstr[subA.begin - 1] == subB.seq.cstr[subB.begin + i])
      {
	 visA[i] = subA.seq.str[--subA.begin];
	 subA.len++;

	 bases++;
      }

      i--;
   }

   if (++i > 0)
   {
      // trim excess initial bases
      visA = visA.substr(i);
      visB = visB.substr(i);
      subB.begin += i;
      subB.len   -= i;
   }
}

//------------------------------------------------------------------------------------
// adjustTrailing() replaces each trailing spacer in visA with a matching base if
// there is one, and excess trailing bases are trimmed

static void adjustTrailing(const int trailing, const int rbases, Sub& subA, Sub& subB,
                           String& visA, String& visB)
{
   int i = 0, bases = 0;
   while (i < trailing && bases < rbases)
   {
      if (subA.seq.cstr[subA.end()] == subB.seq.cstr[subB.end() - trailing + i])
      {
	 visA[visA.length() - trailing + i] = subA.seq.str[subA.end()];
	 subA.len++;

	 bases++;
      }

      i++;
   }

   if (i < trailing)
   {
      // trim excess trailing bases
      visA = visA.substr(0, visA.length() - trailing + i);
      visB = visB.substr(0, visB.length() - trailing + i);
      subB.len -= trailing - i;
   }
}

//------------------------------------------------------------------------------------
// Align::adjust() adjusts the alignment when there are initial or trailing spacers;
// these spacers represent insertions and deletions at the left or right end of the
// alignment

void Align::adjust(const int lbases1, const int lbases2,
                   const int rbases1, const int rbases2)
{
   int initial = initialSpacers(vis1);
   if (initial > 0)
      adjustInitial(initial, lbases1, sub1, sub2, vis1, vis2);
   else
   {
      initial = initialSpacers(vis2);
      if (initial > 0)
         adjustInitial(initial, lbases2, sub2, sub1, vis2, vis1);
   }

   int trailing = trailingSpacers(vis1);
   if (trailing > 0)
      adjustTrailing(trailing, rbases1, sub1, sub2, vis1, vis2);
   else
   {
      trailing = trailingSpacers(vis2);
      if (trailing > 0)
         adjustTrailing(trailing, rbases2, sub2, sub1, vis2, vis1);
   }
}

//------------------------------------------------------------------------------------
// Align::getStats() gets the possible and actual number of matching bases in the
// alignment

void Align::getStats(int& possible, int& matches) const
{
   possible = vis1.length();
   if (possible != vis2.length()) // this should never happen
      throw std::runtime_error("internal vis error");

   matches = 0;
   int j1 = sub1.begin, j2 = sub2.begin;
   const int end1 = sub1.end(), end2 = sub2.end();

   for (int i = 0; i < possible; i++)
      if (vis1[i] == SPACER)
         j2++;
      else if (vis2[i] == SPACER)
         j1++;
      else if (j1 < end1 && j2 < end2)
      {
         if (sub1.seq.cstr[j1++] == sub2.seq.cstr[j2++])
            matches++;
      }
      else // this should never happen
         throw std::runtime_error("internal vis error");

   if (j1 != end1 || j2 != end2) // this should never happen
      throw std::runtime_error("internal vis error");
}

//------------------------------------------------------------------------------------
// connectReadToPatterns() identifies patterns that share minimizers with a read;
// origins suggested by shared minimizers are tabulated in the resulting
// PatternOriginMap; the eligible vector is indexed by pattern index and indicates
// which patterns are to be considered; if this vector is not provided, all patterns
// are considered

void connectReadToPatterns(const WindowVector& readWindowVector,
                           const Minimizer maxMinimizer, const PatternMap& pmap,
			   const BoolVector& inPmap, PatternOriginMap& omap,
			   const BoolVector *eligible)
{
   const int numReadWindows = readWindowVector.size();

   for (int i = 0; i < numReadWindows; i++)
   {
      const Window&   readWindow = readWindowVector[i];
      const Minimizer minimizer  = readWindow.minimizer;

      if (minimizer > maxMinimizer)
         continue; // ignore common minimizer

      // a large percentage of read minimizers will not be found in the pattern map;
      // it is much faster to check for this using the Boolean vector than to call
      // pmap.find()
      if (!inPmap[minimizer])
         continue; // this read minimizer is not found in any pattern

      PatternMap::const_iterator ppos = pmap.find(minimizer);
      if (ppos == pmap.cend()) // this should never happen
         continue;

      const PlocDuo& pduo = ppos->second;

      for (int side = 0; side < 2; side++) // 0 for left/unsided; 1 for right side
      {
         const int numPloc = pduo[side].size();

	 for (int j = 0; j < numPloc; j++)
	 {
            const Ploc& ploc = pduo[side][j];

	    if (eligible && !(*eligible)[ploc.pindex])
               continue;

	    PatternOriginMap::iterator opos = omap.find(ploc.pindex);
	    if (opos == omap.end())
               opos =
                  omap.insert(std::make_pair(ploc.pindex, OriginCounterDuo())).first;

	    OriginCounterDuo& oduo = opos->second;
	    oduo[side].incrementCount(getOrigin(ploc.poffset, readWindow.offset));
	 }
      }
   }
}

//------------------------------------------------------------------------------------
// getOverlap() determines the longest overlapping substrings (of equal length) of
// sub1 and sub2 for the given origin

void getOverlap(const Sub& sub1, const Sub& sub2, const Origin origin,
                int& begin1, int& begin2, int& len,
		int& lbases1, int& lbases2, int& rbases1, int& rbases2)
{
   begin1 = std::max(sub1.begin, sub2.begin + origin);
   begin2 = begin1 - origin;

   len = std::min(sub1.end() - begin1, sub2.end() - begin2);
   if (len < 0)
      len = 0;

   // compute number of unused bases before the overlapping substrings
   lbases1 = begin1 - sub1.begin;
   lbases2 = begin2 - sub2.begin;

   // compute number of unused bases after the overlapping substrings
   rbases1 = sub1.end() - (begin1 + len);
   rbases2 = sub2.end() - (begin2 + len);
}

//------------------------------------------------------------------------------------

class BaseCounts // number of A's, C'S, G's, T's in substring
{
public:
   BaseCounts(const Sub& sub)
      : acount(0), ccount(0), gcount(0), tcount(0)
   {
      const char *cstr = sub.cstr();
      const int   len  = sub.len;

      for (int i = 0; i < len; i++)
         switch (cstr[i])
	 {
            case 'A': acount++; break;
            case 'C': ccount++; break;
            case 'G': gcount++; break;
            case 'T': tcount++; break;
            default : break;
	 }
   }

   int maxMatches(const BaseCounts& other) const
   {
      return std::min(acount, other.acount) +
             std::min(ccount, other.ccount) +
	     std::min(gcount, other.gcount) +
	     std::min(tcount, other.tcount);
   }

   int acount, ccount, gcount, tcount;
};

//------------------------------------------------------------------------------------
// generateVis() generates a visualization of the alignment of two substrings; the
// matrix from the LCS computation is traversed from lower right to upper left, and
// the visualization strings are constructed from right to left

static void generateVis(const Sub& sub1, const Sub& sub2, int **matrix,
                        String& vis1, String& vis2)
{
   const String str1 = sub1.str();
   const String str2 = sub2.str();

   const char *cstr1 = sub1.cstr();
   const char *cstr2 = sub2.cstr();

   const int len1    = sub1.len;
   const int len2    = sub2.len;
   const int buflen  = len1 + len2;

   char *buf1  = new char[buflen + 1];
   char *buf2  = new char[buflen + 1];

   int k = buflen;

   buf1[k] = 0; // null terminating character
   buf2[k] = 0;

   int i = len1, j = len2;

   while (i > 0 || j > 0)
   {
      k--;

      if (i > 0 && j > 0 && (cstr1[i - 1] == cstr2[j - 1] ||
          matrix[i - 1][j - 1] >= std::max(matrix[i][j - 1], matrix[i - 1][j])))
      {
         // match or substitution
	 buf1[k] = str1[--i];
	 buf2[k] = str2[--j];
      }
      else if (i == 0 || j > 0 && matrix[i][j - 1] >= matrix[i - 1][j])
      {
         buf1[k] = SPACER;
	 buf2[k] = str2[--j];
      }
      else
      {
         buf1[k] = str1[--i];
	 buf2[k] = SPACER;
      }
   }

   vis1 = &buf1[k];
   vis2 = &buf2[k];

   delete[] buf1;
   delete[] buf2;
}

//------------------------------------------------------------------------------------
// alignSubstrings() computes a Longest Common Subsequence (LCS) of two substrings
// using a classic dynamic programming algorithm; if the length of the LCS is at least
// minMatches, an Align object is returned (it is the caller's obligation to
// de-allocate it); otherwise, nullptr is returned

Align *alignSubstrings(const Sub& sub1, const Sub& sub2, const int minMatches)
{
   if (minMatches > 0)
   {
      const BaseCounts baseCounts1(sub1), baseCounts2(sub2);
      if (baseCounts1.maxMatches(baseCounts2) < minMatches)
         return nullptr;
   }

   const char *cstr1 = sub1.cstr(); const int len1 = sub1.len;
   const char *cstr2 = sub2.cstr(); const int len2 = sub2.len;

   int **matrix  = new int *[len1 + 1]();
   int *previous = matrix[0] = new int[len2 + 1](); // init first row to all zeros

   for (int i = 1; i <= len1; i++) // for each row
   {
      int *current = matrix[i] = new int[len2 + 1];
      current[0] = 0;

      const int irem = len1 - i;
      bool  mayAbort = (std::min(irem, len2) < minMatches);

      int last = 0, diag = 0;

      for (int j = 1; j <= len2; j++) // for each column
      {
         const int above = previous[j];

	 current[j] = last =
            (cstr1[i - 1] == cstr2[j - 1] ? diag + 1 : std::max(last, above));

	 diag = above;

	 if (mayAbort && last + std::min(irem, len2 - j) >= minMatches)
            mayAbort = false;
      }

      if (mayAbort)
         break; // unable to obtain at least minMatches, so abort the process

      previous = current;
   }

   Align *align = nullptr;

   if (matrix[len1] && matrix[len1][len2] >= minMatches) // success
   {
      String vis1, vis2;
      generateVis(sub1, sub2, matrix, vis1, vis2);

      align = new Align(sub1, sub2, vis1, vis2);
   }

   // cleanup
   int i = 0;
   while (i <= len1 && matrix[i])
      delete[] matrix[i++];

   delete[] matrix;

   return align;
}

//------------------------------------------------------------------------------------
// alignByOrigin() returns the alignment of the longest overlapping substrings (of
// equal length) of sub1 and sub2 for the given origin (it is the caller's obligation
// to de-allocate it), but nullptr is returned if there are not enough matching bases
// in the alignment

static Align *alignByOrigin(const Sub& sub1, const Sub& sub2, const Origin origin,
                            const double minPercentAgreement)
{
   int begin1, begin2, len, lbases1, lbases2, rbases1, rbases2;
   getOverlap(sub1, sub2, origin, begin1, begin2, len, lbases1, lbases2,
              rbases1, rbases2);

   const Sub ovSub1(sub1.seq, begin1, len);
   const Sub ovSub2(sub2.seq, begin2, len);

   int minMatches = computeMinMatches(len, minPercentAgreement);

   // the number of matches might be increased by the adjustment, so let us not be too
   // strict when aligning the substrings
   minMatches = std::max(static_cast<int>(TRIMER_LEN), minMatches - 2);

   Align *align = alignSubstrings(ovSub1, ovSub2, minMatches);

   if (align)
      align->adjust(lbases1, lbases2, rbases1, rbases2);

   return align;
}

//------------------------------------------------------------------------------------
// alignByEachOrigin() invokes alignByOrigin() for each origin in the given vector;
// alignments are saved in the AlignVector (it is the caller's obligation to
// de-allocate them) and the number of alignments is returned

static int alignByEachOrigin(const Sub& sub1, const Sub& sub2,
                             const OriginVector& ovector,
			     const double minPercentAgreement, AlignVector& avector)
{
   int numAligns  = 0;
   const int numOrigins = ovector.size();

   for (int i = 0; i < numOrigins; i++)
   {
      Align *align = alignByOrigin(sub1, sub2, ovector[i], minPercentAgreement);

      if (align)
      {
         avector.push_back(align);
	 numAligns++;
      }
   }

   return numAligns;
}

//------------------------------------------------------------------------------------
// alignReadToPatterns() aligns a read to patterns that share minimizers with the read
// using the most frequently occurring origins, which are likely to produce the best
// alignments; the alignments, if any, are stored in the PatternAlignMap; it is the
// caller's obligation to de-allocate them

void alignReadToPatterns(const Seq& readSeq, const PatternVector& pvector,
                         const double minPercentAgreement,
			 const PatternOriginMap& omap, PatternAlignMap& amap)
{
   for (PatternOriginMap::const_iterator opos = omap.cbegin();
        opos != omap.cend(); ++opos)
   {
      const int pindex = opos->first;
      const Pattern& pattern = *pvector[pindex];
      const OriginCounterDuo& oduo = opos->second;

      for (int side = 0; side < 2; side++) // 0 for left/unsided; 1 for right side
      {
         const int maxCount = oduo[side].maxCount(); // max #occurrences of an origin
	 if (maxCount == 0)
            continue; // nothing on this side

	 OriginVector ovector;
	 oduo[side].retrieveOrigins(maxCount / 2, ovector);

	 const Seq& patternSeq = (side == 0 ? *pattern.left : *pattern.right);
	 const Sub patternSub(patternSeq);
	 const Sub readSub(readSeq);
	 AlignVector avector;

	 if (alignByEachOrigin(patternSub, readSub, ovector, minPercentAgreement,
                               avector) > 0) // got one or more alignments
	 {
            PatternAlignMap::iterator apos = amap.find(pindex);
	    if (apos == amap.end())
               apos = amap.insert(std::make_pair(pindex, AlignDuo())).first;

	    AlignDuo& aduo = apos->second;
	    aduo[side] = avector; // save the alignments here
	 }
      }
   }
}
