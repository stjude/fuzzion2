//------------------------------------------------------------------------------------
//
// match.cpp - module supporting pattern matching
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "match.h"

const char BLANK      = ' ';
const int  PENALTY    = 3;    // penalty charged for each mismatched base in scoring
const int  MISALIGNED = 9999; // value to indicate read2 is misaligned ahead of read1

//------------------------------------------------------------------------------------
// computeScore() returns the score for the given number of possible and actual
// matching bases

static int computeScore(const int possible, const int matches)
{
   return (matches - PENALTY * (possible - matches));
}

//------------------------------------------------------------------------------------
// SingleMatch::SingleMatch() makes copies of the given alignments

SingleMatch::SingleMatch(const Align *inLeft, const Align *inMiddle,
                         const Align *inRight, const bool includeMiddleStats)
   : left(nullptr), middle(nullptr), right(nullptr), possible(0), matches(0),
     loverlap(0), roverlap(0), qual(MatchQual::UNDETERMINED)
{
   if (inLeft)
   {
      left = new Align(*inLeft);
      int lpossible, lmatches;
      left->getStats(lpossible, lmatches);
      possible += lpossible;
      matches  += lmatches;
      loverlap  = left->sub1.len;
   }

   if (inMiddle)
   {
      middle = new Align(*inMiddle);

      if (includeMiddleStats)
      {
         int mpossible, mmatches;
         middle->getStats(mpossible, mmatches);
         possible += mpossible;
         matches  += mmatches;
      }
   }

   if (inRight)
   {
      right = new Align(*inRight);
      int rpossible, rmatches;
      right->getStats(rpossible, rmatches);
      possible += rpossible;
      matches  += rmatches;
      roverlap  = right->sub1.len;
   }

   score = computeScore(possible, matches);
}

//------------------------------------------------------------------------------------
// SingleMatch::SingleMatch() copy constructor

SingleMatch::SingleMatch(const SingleMatch& other)
   : possible(other.possible), matches(other.matches), score(other.score),
     loverlap(other.loverlap), roverlap(other.roverlap), qual(other.qual)
{
   left   = (other.left   ? new Align(*other.left)   : nullptr);
   middle = (other.middle ? new Align(*other.middle) : nullptr);
   right  = (other.right  ? new Align(*other.right)  : nullptr);
}

//------------------------------------------------------------------------------------
// SingleMatch::sufficientAgeement() returns true if there is a sufficient number of
// matching bases

bool SingleMatch::sufficientAgreement(const double minPercentAgreement) const
{
   return (matches >= computeMinMatches(possible, minPercentAgreement));
}

//------------------------------------------------------------------------------------
// SingleMatch::spanning() returns true if this match qualifies as "spanning," i.e.,
// both sides are sufficiently matched

bool SingleMatch::spanning(const double minPercentAgreement,
                           const int minOverlap) const
{
   if (left && right && loverlap >= minOverlap && roverlap >= minOverlap)
   {
      int lpossible, lmatches, rpossible, rmatches;
      left ->getStats(lpossible, lmatches);
      right->getStats(rpossible, rmatches);

      return (lmatches >= computeMinMatches(lpossible, minPercentAgreement) &&
              rmatches >= computeMinMatches(rpossible, minPercentAgreement));
   }
   else
      return false;
}

//------------------------------------------------------------------------------------
// SingleMatch::betterThan() returns true if this match is a better match than the
// given match

bool SingleMatch::betterThan(const SingleMatch& other) const
{
   if (score == other.score)
      return (matches > other.matches);
   else
      return (score > other.score);
}

//------------------------------------------------------------------------------------
// leftCheck() returns true if one of these cases is detected:
//              insertion   deletion
//    pattern:  -b...b(bX   bb...b(X
//    read:     bb...b(X    -b...b(bX

static bool leftCheck(const Align& left, char base, bool insertion)
{
   const char *cstr1 = left.sub1.cstr();
   const char *cstr2 = left.sub2.cstr();

   int i1 = left.sub1.len - 1;
   int i2 = left.sub2.len - 1;
   int v  = left.vis1.length() - 1;

   while (i1 >= 0 && cstr1[i1] == base &&
          i2 >= 0 && cstr2[i2] == base)
   {
      i1--;
      i2--;
      v--;
   }

   if (insertion)
      return (i2 >= 0 && cstr2[i2] == base && v >= 0 && left.vis1[v] == SPACER);
   else // deletion
      return (i1 >= 0 && cstr1[i1] == base && v >= 0 && left.vis2[v] == SPACER);
}

//------------------------------------------------------------------------------------
// SingleMatch::leftRightQual() returns the qualification indicator for a spanning
// match of a pattern of this type: left(ex)right

MatchQual SingleMatch::leftRightQual(const Pattern& pattern)
{
   const char *m    = middle->sub2.cstr();
   const int   mlen = middle->sub2.len;
   const int   n    = pattern.xvector->size();

   // look for an exact match of an excluded middle
   for (int i = 0; i < n; i++)
   {
      const Seq&  xseq = (*pattern.xvector)[i];
      const char *x    = xseq.cstr;
      const int   xlen = xseq.len;

      if (xlen == mlen     && std::strncmp(x, m, xlen) == 0 ||
          xlen == mlen + 1 && std::strncmp(&x[1], m, mlen) == 0 &&
	  leftCheck(*left, x[0], true) ||
	  xlen == mlen - 1 && std::strncmp(x, &m[1], xlen) == 0 &&
	  leftCheck(*left, m[0], false))
         return MatchQual::DISQUALIFIED;
   }

   // full match not found
   return MatchQual::QUALIFIED;
}

//------------------------------------------------------------------------------------
// SingleMatch::setQual() sets an indicator of whether the match qualifies or
// disqualifies as a hit, or whether it is undetermined

void SingleMatch::setQual(const Pattern& pattern, const double minPercentAgreement,
                          const int minOverlap)
{
   if (!pattern.xvector)
      if (loverlap >= minOverlap && (!pattern.right ||
          spanning(minPercentAgreement, minOverlap)))
         qual = MatchQual::QUALIFIED;
      else
         qual = MatchQual::UNDETERMINED;
   else
      if (left && middle && right)
      {
         qual = leftRightQual(pattern);

	 if (qual == MatchQual::QUALIFIED &&
             !spanning(minPercentAgreement, minOverlap))
            qual = MatchQual::UNDETERMINED;
      }
      else
         qual = MatchQual::UNDETERMINED;
}

//------------------------------------------------------------------------------------
// truncate() truncates the given string on the left or right to a new length; the
// length of the string after truncation is returned

static int truncate(const int newlen, const bool truncateRight, String& s)
{
   const int oldlen = s.length();
   if (oldlen <= newlen) // truncation is not needed
      return oldlen;

   if (truncateRight)
      s = s.substr(0, newlen - 1) + ELLIPSIS_MARK;
   else // truncate on the left
      s = ELLIPSIS_MARK + s.substr(oldlen - newlen + 1);

   return newlen;
}

//------------------------------------------------------------------------------------
// pad() pads the given string with blanks on the left or right to a new length; the
// length of the string after padding is returned

static int pad(const int newlen, const bool padRight, String& s)
{
   const int oldlen = s.length();
   if (oldlen >= newlen) // padding is not needed
      return oldlen;

   if (padRight)
      s += String(newlen - oldlen, BLANK);
   else // pad on the left
      s  = String(newlen - oldlen, BLANK) + s;

   return newlen;
}

//------------------------------------------------------------------------------------
// ellide() returns the ELLIPSIS_MARK if the given expression is true; otherwise,
// BLANK is returned

static char ellide(const bool elliding)
{
   return (elliding ? ELLIPSIS_MARK : BLANK);
}

//------------------------------------------------------------------------------------
// SingleMatch::getMiddleVis() prepares a visualization of the middle of a match; if
// the middle is unbounded on one side, it is limited to maxmidlen

void SingleMatch::getMiddleVis(const Pattern& pattern, const int maxmidlen,
                               String& mvis1, String& mvis2) const
{
   if (pattern.middle)
      mvis1 = middle->vis1;
   else if (pattern.xvector)
      mvis1 = pattern.xvis;
   else
      mvis1 = "";

   mvis2 = middle->vis2;

   int mlen1 = mvis1.length();
   int mlen2 = mvis2.length();

   if (mlen1 != mlen2)
   {
      if (!left || !right)
         mlen2 = truncate(std::max(mlen1, maxmidlen), (left ? true : false), mvis2);

      if (mlen1 < mlen2)
         pad(mlen2, (left ? true : false), mvis1);
      else if (mlen1 > mlen2)
         pad(mlen1, (left ? true : false), mvis2);
   }
}

//------------------------------------------------------------------------------------
// SingleMatch::getVis() prepares a visualization of the match; if the middle is
// unbounded on one side, it is limited to maxmidlen

void SingleMatch::getVis(const Pattern& pattern, const int maxmidlen,
                         String& vis1, String& vis2) const
{
   String lvis1 = "", mvis1 = "", rvis1 = "";
   String lvis2 = "", mvis2 = "", rvis2 = "";

   if (left)
   {
      if (left->sub1.hasBefore() || left->sub2.hasBefore())
      {
         lvis1 = ellide(left->sub1.hasBefore()) + left->vis1;
         lvis2 = ellide(left->sub2.hasBefore()) + left->vis2;
      }
      else
      {
         lvis1 = left->vis1;
         lvis2 = left->vis2;
      }

      if (middle)
      {
         // append first delimiter to left visualization
         const char delim = pattern.sequence[pattern.delim1];
	 lvis1 += delim;
	 lvis2 += delim;
      }

      // may need to indicate ellide on the right if this is an unsided match
      if (!middle && !right && (left->sub1.hasAfter() || left->sub2.hasAfter()))
      {
         lvis1 += ellide(left->sub1.hasAfter());
         lvis2 += ellide(left->sub2.hasAfter());
      }
   }

   if (middle)
      getMiddleVis(pattern, maxmidlen, mvis1, mvis2);

   if (right)
   {
      if (right->sub1.hasAfter() || right->sub2.hasAfter())
      {
         rvis1 = right->vis1 + ellide(right->sub1.hasAfter());
         rvis2 = right->vis2 + ellide(right->sub2.hasAfter());
      }
      else
      {
         rvis1 = right->vis1;
         rvis2 = right->vis2;
      }

      if (middle)
      {
         // prepend second delimiter to right visualization
         const char delim = pattern.sequence[pattern.delim2];
	 rvis1 = delim + rvis1;
	 rvis2 = delim + rvis2;
      }
   }

   vis1 = lvis1 + mvis1 + rvis1;
   vis2 = lvis2 + mvis2 + rvis2;
}

//------------------------------------------------------------------------------------
// save() constructs a match from the given alignments and saves it in the given
// vector

static void save(const Align *left, const Align *middle, const Align *right,
                 const bool includeMiddleStats, SingleMatchVector& mvector)
{
   mvector.push_back(new SingleMatch(left, middle, right, includeMiddleStats));
}

//------------------------------------------------------------------------------------
// saveLeftRight() constructs a match by connecting the given left and right
// alignments and saves it in the given vector; however, false is returned if the left
// and right alignments are incompatible

static bool saveLeftRight(const Align *left, const Align *right,
                          const Seq *patternMiddle, SingleMatchVector& mvector)
{
   const int mbegin = left->sub2.end();
   const int mlen   = right->sub2.begin - mbegin;

   if (mlen < 0) // left and right alignments are incompatible
      return false;

   const Sub sub2(left->sub2.seq, mbegin, mlen);
   Align *middle;

   if (patternMiddle)
   {
      const Sub sub1(*patternMiddle);
      middle = alignSubstrings(sub1, sub2, 0);
   }
   else // the pattern has no middle
   {
      const Sub sub1(left->sub1.seq, left->sub1.len, 0); // zero-length substring
      middle = new Align(sub1, sub2, "", sub2.str());
   }

   save(left, middle, right, (patternMiddle ? true : false), mvector);
   delete middle;
   return true;
}

//------------------------------------------------------------------------------------
// countOrigins() counts origins obtained from shared 3-mers

static void countOrigins(const TrimerOffsets& offsets1,
                         const TrimerOffsets& offsets2, OriginCounter& ocounter)
{
   int n1, n2;

   for (int t = 0; t < NUM_TRIMERS; t++) // for each 3-mer
      if ((n1 = offsets1[t].size()) > 0 &&
          (n2 = offsets2[t].size()) > 0)
         for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
	    {
               const Origin origin = getOrigin(offsets1[t][i], offsets2[t][j]);
	       ocounter.incrementCount(origin);
	    }
}

//------------------------------------------------------------------------------------
// saveLeftExtend() constructs matches by extending the given left alignment to the
// middle and right and saves them in the given vector

static void saveLeftExtend(const Align *left, const Pattern& pattern,
                           SingleMatchVector& mvector)
{
   const Sub afterSub(left->sub2, false);

   if (afterSub.len >= TRIMER_LEN && pattern.right)
   {
      // try extending to the middle and right based on shared 3-mers
      TrimerOffsets afterTrimers;
      initTrimerOffsets(afterSub.seq.cstr, afterSub.begin, afterSub.len,
                        afterTrimers);

      OriginCounter ocounter;
      countOrigins(pattern.rightTrimers, afterTrimers, ocounter);

      const int maxCount = ocounter.maxCount();
      if (maxCount > 0)
      {
         const Sub patternSub(*pattern.right);
	 OriginVector ovector;
	 const int numOrigins = ocounter.retrieveOrigins(maxCount - 2, ovector);

	 for (int i = 0; i < numOrigins; i++)
	 {
	    int begin1, begin2, len, lbases1, lbases2, rbases1, rbases2;
	    getOverlap(patternSub, afterSub, ovector[i], begin1, begin2, len,
                       lbases1, lbases2, rbases1, rbases2);
	    if (begin1 > len)
               continue; // too many deletions at left end of pattern right

	    const Sub ovSub1(patternSub.seq, begin1, len);
	    const Sub ovSub2(afterSub.seq,   begin2, len);

	    const int minMatches = std::max(static_cast<int>(TRIMER_LEN), len - 5);

	    Align *right = alignSubstrings(ovSub1, ovSub2, minMatches);
	    if (right)
	    {
               right->adjust(lbases1, lbases2, rbases1, rbases2);

               const int deletions = right->sub1.begin;
	       if (deletions > 0)
	       {
                  // prepend deletions to alignment
		  right->sub1.begin = 0;
		  right->sub1.len  += deletions;
		  right->vis1 = patternSub.seq.str.substr(0, deletions) + right->vis1;
		  right->vis2 = String(deletions, SPACER) + right->vis2;
	       }

	       saveLeftRight(left, right, pattern.middle, mvector);
	       delete right;
	    }
	 }
      }
   }

   // try extending to the middle only
   Align *middle;

   if (pattern.middle)
   {
      const int len = std::min(pattern.middle->len, afterSub.len);
      const Sub sub1(*pattern.middle, 0, len);
      const int minMatches = std::max(0, len - 5);
      middle = alignSubstrings(sub1, afterSub, minMatches);
   }
   else // the pattern has no middle
   {
      const Sub sub1(left->sub1.seq, left->sub1.len, 0); // zero-length substring
      middle = new Align(sub1, afterSub, "", afterSub.str());
   }

   if (middle)
   {
      save(left, middle, nullptr, (pattern.middle ? true : false), mvector);
      delete middle;
   }
}

//------------------------------------------------------------------------------------
// saveRightExtend() constructs matches by extending the given right alignment to the
// middle and left and saves them in the given vector

static void saveRightExtend(const Align *right, const Pattern& pattern,
                            SingleMatchVector& mvector)
{
   const Sub beforeSub(right->sub2, true);

   if (beforeSub.len >= TRIMER_LEN)
   {
      // try extending to the middle and left based on shared 3-mers
      TrimerOffsets beforeTrimers;
      initTrimerOffsets(beforeSub.seq.cstr, beforeSub.begin, beforeSub.len,
                        beforeTrimers);

      OriginCounter ocounter;
      countOrigins(pattern.leftTrimers, beforeTrimers, ocounter);

      const int maxCount = ocounter.maxCount();
      if (maxCount > 0)
      {
         const Sub patternSub(*pattern.left);
	 OriginVector ovector;
	 const int numOrigins = ocounter.retrieveOrigins(maxCount - 2, ovector);

	 for (int i = 0; i < numOrigins; i++)
	 {
	    int begin1, begin2, len, lbases1, lbases2, rbases1, rbases2;
	    getOverlap(patternSub, beforeSub, ovector[i], begin1, begin2, len,
                       lbases1, lbases2, rbases1, rbases2);
	    if (patternSub.seq.len - (begin1 + len) > len)
               continue; // too many deletions at right end of pattern left

	    const Sub ovSub1(patternSub.seq, begin1, len);
	    const Sub ovSub2(beforeSub.seq,  begin2, len);

	    const int minMatches = std::max(static_cast<int>(TRIMER_LEN), len - 5);

	    Align *left = alignSubstrings(ovSub1, ovSub2, minMatches);
	    if (left)
	    {
               left->adjust(lbases1, lbases2, rbases1, rbases2);

               const int deletions = patternSub.seq.len - left->sub1.end();
	       if (deletions > 0)
	       {
                  // append deletions to alignment
		  left->sub1.len += deletions;
		  left->vis1 +=
                     patternSub.seq.str.substr(patternSub.seq.len - deletions);
		  left->vis2 += String(deletions, SPACER);
	       }

	       saveLeftRight(left, right, pattern.middle, mvector);
	       delete left;
	    }
	 }
      }
   }

   // try extending to the middle only
   Align *middle;

   if (pattern.middle)
   {
      const int len   = std::min(pattern.middle->len, beforeSub.len);
      const int begin = pattern.middle->len - len;
      const Sub sub1(*pattern.middle, begin, len);
      const int minMatches = std::max(0, len - 5);
      middle = alignSubstrings(sub1, beforeSub, minMatches);
   }
   else // the pattern has no middle
   {
      const Sub sub1(right->sub1.seq, 0, 0); // zero-length substring
      middle = new Align(sub1, beforeSub, "", beforeSub.str());
   }

   if (middle)
   {
      save(nullptr, middle, right, (pattern.middle ? true : false), mvector);
      delete middle;
   }
}

//------------------------------------------------------------------------------------
// findMatches() examines the given alignments and finds matches

static void findMatches(AlignVector& lvector, AlignVector& rvector,
                        const Pattern& pattern, SingleMatchVector& mvector)
{
   const int lcount = lvector.size();
   const int rcount = rvector.size();

   if (!pattern.right) // the pattern is unsided
      for (int i = 0; i < lcount; i++)
         save(lvector[i], nullptr, nullptr, false, mvector);
   else
   {
      // look for one-sided matches
      for (int i = 0; i < lcount; i++)
         if (!lvector[i]->sub2.hasAfter())
            save(lvector[i], nullptr, nullptr, false, mvector);

      for (int i = 0; i < rcount; i++)
         if (!rvector[i]->sub2.hasBefore())
            save(nullptr, nullptr, rvector[i], false, mvector);

      // now look for two-sided matches
      bool foundCompatible = false;

      for (int i = 0; i < lcount; i++)
         if (lvector[i]->sub2.hasAfter())
            for (int j = 0; j < rcount; j++)
               if (rvector[j]->sub2.hasBefore() &&
                   saveLeftRight(lvector[i], rvector[j], pattern.middle, mvector))
                  foundCompatible = true;

      if (!foundCompatible)
      {
         // try to extend alignments
	 for (int i = 0; i < lcount; i++)
            if (lvector[i]->sub2.hasAfter())
               saveLeftExtend(lvector[i], pattern, mvector);

	 for (int i = 0; i < rcount; i++)
            if (rvector[i]->sub2.hasBefore())
               saveRightExtend(rvector[i], pattern, mvector);
      }
   }

   //cleanup
   for (int i = 0; i < lcount; i++)
      delete lvector[i];

   for (int i = 0; i < rcount; i++)
      delete rvector[i];
}

//------------------------------------------------------------------------------------
// getSingleMatches() finds matches of each pattern; the matches may be designated as
// QUALIFIED or UNDETERMINED; the number of these matches is returned and if greater
// than zero, readSeq and the matches in mmap must be de-allocated by the caller; if
// mateMap is non-null, it identifies a subset of patterns to consider (otherwise,
// all patterns are considered)

int getSingleMatches(const String& readStr, const MinimizerWindowLength w,
                     const KmerRankTable *rankTable, const Minimizer maxMinimizer,
		     const PatternMap& pmap, const BoolVector& inPmap,
		     const PatternVector& pvector, const double minPercentAgreement,
		     const int minOverlap, Seq *& readSeq, SingleMatchMap& mmap,
		     const SingleMatchMap *mateMap)
{
   WindowVector readWindowVector;
   getWindows(readStr.c_str(), readStr.length(), w, rankTable, readWindowVector);

   PatternOriginMap omap;

   if (mateMap)
   {
      BoolVector *eligible = new BoolVector(pvector.size(), false);

      for (SingleMatchMap::const_iterator mpos = mateMap->cbegin();
           mpos != mateMap->cend(); ++mpos)
         (*eligible)[mpos->first] = true;

      connectReadToPatterns(readWindowVector, maxMinimizer, pmap, inPmap, omap,
                            eligible);

      delete eligible;
   }
   else
      connectReadToPatterns(readWindowVector, maxMinimizer, pmap, inPmap, omap);

   if (omap.size() == 0) // found no connection to any eligible pattern
   {
      readSeq = nullptr;
      return 0;
   }

   bool invalidBase;
   readSeq = new Seq(readStr, invalidBase);

   PatternAlignMap amap;
   alignReadToPatterns(*readSeq, pvector, minPercentAgreement, omap, amap);

   int numMatches =  0;

   for (PatternAlignMap::iterator apos = amap.begin(); apos != amap.end(); ++apos)
   {
      const int pindex = apos->first;
      const Pattern& pattern = *pvector[pindex];
      AlignDuo& aduo = apos->second;
      SingleMatchVector mvector, mvectorAgree;

      findMatches(aduo[0], aduo[1], pattern, mvector);

      const int count = mvector.size();
      if (count == 0) // no matches
         continue;

      // delete matches with insufficient agreement

      for (int i = 0; i < count; i++)
         if (!mvector[i]->sufficientAgreement(minPercentAgreement))
	 {
            delete mvector[i];
	    mvector[i] = nullptr;
	 }

      // determine whether the remaining matches are qualified, undetermined, or
      // disqualifed

      bool foundDisqualified = false;

      for (int i = 0; i < count; i++)
         if (mvector[i])
	 {
            mvector[i]->setQual(pattern, minPercentAgreement, minOverlap);

	    if (mvector[i]->qual == MatchQual::DISQUALIFIED)
	    {
               foundDisqualified = true;
	       break;
	    }
	 }

      // if one is disqualified, all are disqualified; retain only those matches that
      // are qualified or undetermined

      if (foundDisqualified)
         for (int i = 0; i < count; i++)
            delete mvector[i];
      else
      {
         SingleMatchVector mvectorAgree;
	 for (int i = 0; i < count; i++)
            if (mvector[i])
               mvectorAgree.push_back(mvector[i]);

	 const int numAgree = mvectorAgree.size();
	 if (numAgree > 0)
	 {
            mmap.insert(std::make_pair(pindex, mvectorAgree));
	    numMatches += numAgree;
	 }
      }
   }

   if (numMatches == 0)
   {
      delete readSeq;
      readSeq = nullptr;
   }

   return numMatches;
}

//------------------------------------------------------------------------------------
// selectBestSingleMatches() identifies the best qualified match overall or for each
// pattern, depending on the value of the Boolean parameter

void selectBestSingleMatches(const bool bestOverall, SingleMatchMap& mmap,
                             SingleMatchMap& bestMap)
{
   int pindex;
   SingleMatch *best = nullptr;

   for (SingleMatchMap::iterator mpos = mmap.begin(); mpos != mmap.end(); ++mpos)
   {
      SingleMatchVector& mvector = mpos->second;
      const int count = mvector.size();

      for (int i = 0; i < count; i++)
         if (mvector[i]->qual == MatchQual::QUALIFIED &&
             (!best || mvector[i]->betterThan(*best)))
	 {
            pindex = mpos->first;
	    best   = mvector[i];
	 }

      if (!bestOverall && best)
      {
         SingleMatchVector selected;
	 selected.push_back(new SingleMatch(*best));
	 bestMap.insert(std::make_pair(pindex, selected));
	 best = nullptr; // reset for next pattern
      }
   }

   if (bestOverall && best)
   {
      SingleMatchVector selected;
      selected.push_back(new SingleMatch(*best));
      bestMap.insert(std::make_pair(pindex, selected));
   }
}

//------------------------------------------------------------------------------------
// freeSingleMatchMap() de-allocates the matches in the given map

void freeSingleMatchMap(SingleMatchMap& mmap)
{
   for (SingleMatchMap::iterator mpos = mmap.begin(); mpos != mmap.end(); ++mpos)
   {
      SingleMatchVector& mvector = mpos->second;
      const int count = mvector.size();

      for (int i = 0; i < count; i++)
         delete mvector[i];
   }
}

//------------------------------------------------------------------------------------
// computePairStats() computes for the given matches the insert size and the number of
// bases read2 is aligned ahead of read1 (called misalignment)

static void computePairStats(const Pattern& pattern, const SingleMatch& match1,
                             const SingleMatch& match2, int& insertSize,
			     int& misalignment)
{
   if (!match1.left && match2.left)
   {
      insertSize   = 0;
      misalignment = MISALIGNED;
      return;
   }

   // these variables will hold the inclusive beginning and exclusive ending positions
   // of read1 and read2 with respect to the pattern's left, middle, and right
   int begin1, end1; // for read1
   int begin2, end2; // for read2

   if (match1.left)
   {
      begin1 = match1.left->origin();
      end1   = begin1 + match1.left->sub2.seq.len;

      if (match2.left)
      {
         begin2 = match2.left->origin();
	 end2   = begin2 + match2.left->sub2.seq.len;
      }
      else // match2.right
      {
         if (match1.right)
            begin2 = begin1 + match2.right->origin() - match1.right->origin();
	 else
	 {
            const int mid1   = (match1.middle ? match1.middle->sub2.len : 0);
            const int mid2   = (match2.middle ? match2.middle->sub2.len : 0);
	    const int maxmid = std::max(mid1, mid2);
	    begin2 = pattern.left->len + maxmid + match2.right->origin();
	 }

	 end2 = begin2 + match2.right->sub2.seq.len;
      }
   }
   else // there are right alignments and no left alignments
   {
      begin1 = match1.right->origin();
      end1   = begin1 + match1.right->sub2.seq.len;

      begin2 = match2.right->origin();
      end2   = begin2 + match2.right->sub2.seq.len;
   }

   insertSize   = std::max(end1, end2) - std::min(begin1, begin2);
   misalignment = std::max(0, begin1 - begin2);
}

//------------------------------------------------------------------------------------
// PairMatch::PairMatch() makes copies of the given matches

PairMatch::PairMatch(const Pattern& pattern, const SingleMatch *inMatch1,
                     const SingleMatch *inMatch2)
   : qual(MatchQual::UNDETERMINED)
{
   match1   = new SingleMatch(*inMatch1);
   match2   = new SingleMatch(*inMatch2);

   possible = match1->possible + match2->possible;
   matches  = match1->matches  + match2->matches;
   score    = computeScore(possible, matches);

   computePairStats(pattern, *match1, *match2, insertSize, misalignment);
}

//------------------------------------------------------------------------------------
// PairMatch::PairMatch() copy constructor

PairMatch::PairMatch(const PairMatch& other)
   : possible(other.possible), matches(other.matches), score(other.score),
     insertSize(other.insertSize), misalignment(other.misalignment), qual(other.qual)
{
   match1 = new SingleMatch(*other.match1);
   match2 = new SingleMatch(*other.match2);
}

//------------------------------------------------------------------------------------
// PairMatch::betterThan() returns true if this match is a better match than the
// given match

bool PairMatch::betterThan(const PairMatch& other) const
{
   if (score == other.score)
      if (matches == other.matches)
         return (insertSize < other.insertSize);
      else
         return (matches > other.matches);
   else
      return (score > other.score);
}

//------------------------------------------------------------------------------------
// sideAgreement() returns true if the weighted percentage of matches on one side is
// sufficient agreement

static bool sideAgreement(const Align *side1, const Align *side2,
                          const double minPercentAgreement)
{
   int possible1 = 0, matches1 = 0, possible2 = 0, matches2 = 0;

   if (side1) side1->getStats(possible1, matches1);
   if (side2) side2->getStats(possible2, matches2);

   return (matches1 + matches2 >=
           computeMinMatches(possible1 + possible2, minPercentAgreement));
}

//------------------------------------------------------------------------------------
// PairMatch::setQual() sets an indicator of whether the match qualifies or
// disqualifies as a hit

void PairMatch::setQual(const Pattern& pattern, const double minPercentAgreement,
                        const int minOverlap)
{
   if (pattern.right)
      if (std::max(match1->loverlap, match2->loverlap) >= minOverlap &&
          std::max(match1->roverlap, match2->roverlap) >= minOverlap &&
	  sideAgreement(match1->left,  match2->left,  minPercentAgreement) &&
	  sideAgreement(match1->right, match2->right, minPercentAgreement) &&
	  (match1->qual == MatchQual::QUALIFIED ||
	   match2->qual == MatchQual::QUALIFIED ||
	   pattern.sequence[pattern.delim1] == ']'))
         qual = MatchQual::QUALIFIED;
      else
         qual = MatchQual::DISQUALIFIED;
   else // unsided pattern
      if (match1->qual == MatchQual::QUALIFIED &&
          match2->qual == MatchQual::QUALIFIED)
         qual = MatchQual::QUALIFIED;
      else
         qual = MatchQual::DISQUALIFIED;
}

//------------------------------------------------------------------------------------
// appendVis() is given the alignments of a read pair to a pattern substring and
// appends their visualization to the given string vector; an object is returned
// describing the pattern substring spanned by the visualization (it is the caller's
// obligation to de-allocate it) or nullptr is returned if an error was detected

static Sub *appendVis(const Align *align1, const Align *align2, StringVector& vis)
{
   if (!align1 && !align2 ||
        align1 && align1->vis1.length() != align1->vis2.length() ||
	align2 && align2->vis1.length() != align2->vis2.length())
      return nullptr; // this should never happen

   if (align1 && !align2) // only the first mate is aligned to the pattern
   {
      vis[0] += align1->vis1;
      vis[1] += align1->vis2;
      vis[2] += String(align1->vis1.length(), BLANK);
      return new Sub(align1->sub1);
   }

   if (!align1 && align2) // only the second mate is aligned to the pattern
   {
      vis[0] += align2->vis1;
      vis[1] += String(align2->vis1.length(), BLANK);
      vis[2] += align2->vis2;
      return new Sub(align2->sub1);
   }

   // both mates are aligned to the pattern; we must merge their visualizations

   const int pbegin1 = align1->sub1.begin;
   const int pbegin2 = align2->sub1.begin;

   const int pend1   = align1->sub1.end();
   const int pend2   = align2->sub1.end();

   const int mergedBegin = std::min(pbegin1, pbegin2); // inclusive
   const int mergedEnd   = std::max(pend1, pend2);     // exclusive

   const String& pstr = align1->sub1.seq.str;

   const int len1 = align1->vis1.length();
   const int len2 = align2->vis1.length();

   int i1 = 0, i2 = 0, poffset = mergedBegin;

   while (i1 < len1 || i2 < len2)
   {
      char p1 = (i1 < len1 && poffset >= pbegin1 ? align1->vis1[i1] : BLANK);
      char p2 = (i2 < len2 && poffset >= pbegin2 ? align2->vis1[i2] : BLANK);

      if (p1 == SPACER || p2 == SPACER)
      {
         vis[0] += SPACER;
	 vis[1] += (p1 == SPACER ? align1->vis2[i1++] : BLANK);
	 vis[2] += (p2 == SPACER ? align2->vis2[i2++] : BLANK);
      }
      else
      {
         const char p = pstr[poffset++];

	 if (p1 != BLANK && p1 != p ||
             p2 != BLANK && p2 != p)
            return nullptr; // this should never happen

	 vis[0] += p;
	 vis[1] += (p1 == BLANK ? BLANK : align1->vis2[i1++]);
	 vis[2] += (p2 == BLANK ? BLANK : align2->vis2[i2++]);
      }
   }

   if (poffset != mergedEnd)
      return nullptr; // this should never happen

   return new Sub(align1->sub1.seq, mergedBegin, mergedEnd - mergedBegin);
}

//------------------------------------------------------------------------------------
// appendMiddleVis() is given the matches of a read pair to a pattern and appends
// their visualization of the middle alignment to the given string vector

static void appendMiddleVis(const Pattern& pattern, const int maxmidlen,
                            const SingleMatch& match1, const SingleMatch& match2,
			    StringVector& vis)
{
   if (match1.middle && !match2.middle) // only the first mate's middle is aligned
   {
      String mvis1, mvis2;
      match1.getMiddleVis(pattern, maxmidlen, mvis1, mvis2);

      vis[0] += mvis1;
      vis[1] += mvis2;
      vis[2] += String(mvis1.length(), BLANK);
      return;
   }

   if (!match1.middle && match2.middle) // only the second mate's middle is aligned
   {
      String mvis1, mvis2;
      match2.getMiddleVis(pattern, maxmidlen, mvis1, mvis2);

      vis[0] += mvis1;
      vis[1] += String(mvis1.length(), BLANK);
      vis[2] += mvis2;
      return;
   }

   // each mate's middle is aligned; we must merge their visualizations

   if (pattern.middle)
   {
      const Sub *psub = appendVis(match1.middle, match2.middle, vis);
      if (!psub)
         throw std::runtime_error("internal vis error");

      delete psub;
      return;
   }

   String mvis1 = (pattern.xvector ? pattern.xvis : "");
   String mvis2 = match1.middle->vis2;
   String mvis3 = match2.middle->vis2;

   int mlen1 = mvis1.length();
   int mlen2 = mvis2.length();
   int mlen3 = mvis3.length();

   if (mlen1 != mlen2 || mlen1 != mlen3)
   {
      int maxlen = std::max(mlen1, maxmidlen);

      if (match1.left && match1.right && mlen2 > maxlen)
         maxlen = mlen2;

      if (match2.left && match2.right && mlen3 > maxlen)
         maxlen = mlen3;

      mlen2 = truncate(maxlen, (match1.left ? true : false), mvis2);
      mlen3 = truncate(maxlen, (match2.left ? true : false), mvis3);

      const int eqlen = std::max(mlen1, std::max(mlen2, mlen3));

      pad(eqlen, (match1.left || match2.left ? true : false), mvis1);
      pad(eqlen, (match1.left ? true : false), mvis2);
      pad(eqlen, (match2.left ? true : false), mvis3);
   }

   vis[0] += mvis1;
   vis[1] += mvis2;
   vis[2] += mvis3;
}

//------------------------------------------------------------------------------------
// PairMatch::getVis() prepares a visualization of the match; if the middle is
// unbounded on one side, it is limited to maxmidlen

void PairMatch::getVis(const Pattern& pattern, const int maxmidlen,
                       StringVector& vis) const
{
   vis.clear(); vis.push_back(""); vis.push_back(""); vis.push_back("");

   if (match1->left || match2->left)
   {
      const Sub *psub = appendVis(match1->left, match2->left, vis);
      if (!psub)
         throw std::runtime_error("internal vis error");

      if (match1->middle || match1->right || match2->middle || match2->right)
      {
         if (psub->hasAfter())
	 {
            const String after = psub->seq.str.substr(psub->end());
	    const String blanks(after.length(), BLANK);

	    vis[0] += after;
	    vis[1] += blanks;
	    vis[2] += blanks;
	 }

	 const char delim = pattern.sequence[pattern.delim1];

	 vis[0] += delim;
	 vis[1] += (match1->left && match1->middle ? delim : BLANK);
	 vis[2] += (match2->left && match2->middle ? delim : BLANK);
      }

      if (psub->hasBefore() ||
          match1->left && match1->left->sub2.hasBefore() ||
          match2->left && match2->left->sub2.hasBefore())
      {
         vis[0] = ellide(psub->hasBefore()) + vis[0];
	 vis[1] = ellide(match1->left && match1->left->sub2.hasBefore()) + vis[1];
	 vis[2] = ellide(match2->left && match2->left->sub2.hasBefore()) + vis[2];
      }

      // may need to indicate ellide on the right if this is an unsided match
      if (!match1->middle && !match1->right && !match2->middle && !match2->right &&
          (psub->hasAfter() ||
	   match1->left->sub2.hasAfter() ||
	   match2->left->sub2.hasAfter()))
      {
         vis[0] += ellide(psub->hasAfter());
	 vis[1] += ellide(match1->left->sub2.hasAfter());
	 vis[2] += ellide(match2->left->sub2.hasAfter());
      }

      delete psub;
   }

   if (match1->middle || match2->middle)
      appendMiddleVis(pattern, maxmidlen, *match1, *match2, vis);
   else if (match1->left && match2->right ||
            match2->left && match1->right)
   {
      String mvis = "";

      if (pattern.middle)
         mvis = pattern.middle->str;
      else if (pattern.xvector)
         mvis = pattern.xvis;

      if (mvis.length() > 0)
      {
         const String blanks(mvis.length(), BLANK);

	 vis[0] += mvis;
	 vis[1] += blanks;
	 vis[2] += blanks;
      }
   }

   if (match1->right || match2->right)
   {
      StringVector rvis(3, "");

      const Sub *psub = appendVis(match1->right, match2->right, rvis);
      if (!psub)
         throw std::runtime_error("internal vis error");

      if (match1->left || match1->middle || match2->left || match2->middle)
      {
         const char delim = pattern.sequence[pattern.delim2];

	 vis[0] += delim;
	 vis[1] += (match1->middle && match1->right ? delim : BLANK);
	 vis[2] += (match2->middle && match2->right ? delim : BLANK);

	 if (psub->hasBefore())
	 {
            const String before = psub->seq.str.substr(0, psub->begin);
	    const String blanks(before.length(), BLANK);

	    vis[0] += before;
	    vis[1] += blanks;
	    vis[2] += blanks;
	 }
      }

      if (psub->hasAfter() ||
          match1->right && match1->right->sub2.hasAfter() ||
          match2->right && match2->right->sub2.hasAfter())
      {
         rvis[0] += ellide(psub->hasAfter());
	 rvis[1] += ellide(match1->right && match1->right->sub2.hasAfter());
	 rvis[2] += ellide(match2->right && match2->right->sub2.hasAfter());
      }

      delete psub;

      vis[0] += rvis[0];
      vis[1] += rvis[1];
      vis[2] += rvis[2];
   }
}

//------------------------------------------------------------------------------------
// bestPairMatch() returns the best qualified match for the given pattern, or nullptr
// if none

static PairMatch *bestPairMatch(const Pattern& pattern,
                                const double minPercentAgreement,
				const int minOverlap, const int maxInsert,
				const int maxTrim,
				const SingleMatchVector& mvector1,
				const SingleMatchVector& mvector2)
{
   PairMatch *best = nullptr;
   const int count1 = mvector1.size();
   const int count2 = mvector2.size();

   for (int i = 0; i < count1; i++)
      for (int j = 0; j < count2; j++)
      {
         PairMatch *candidate = new PairMatch(pattern, mvector1[i], mvector2[j]);

	 if (candidate->isPlausible(maxInsert, maxTrim))
	 {
            candidate->setQual(pattern, minPercentAgreement, minOverlap);

	    if (candidate->qual == MatchQual::QUALIFIED &&
                (!best || candidate->betterThan(*best)))
	    {
               delete best;
	       best = candidate;
	    }
	    else
               delete candidate;
	 }
	 else
            delete candidate;
      }

   return best;
}

//------------------------------------------------------------------------------------
// getPairMatches() finds the best qualified match for each pattern that matches the
// given read pair; the number of these matches is returned and if greater than zero,
// readSeq1, readSeq2, and the matches in pairMap must be de-allocated by the caller

int getPairMatches(const String& readStr1, const String& readStr2,
                   const MinimizerWindowLength w, const KmerRankTable *rankTable,
		   const Minimizer maxMinimizer, const PatternMap& pmap,
		   const BoolVector& inPmap, const PatternVector& pvector,
		   const double minPercentAgreement, const int minOverlap,
		   const int maxInsert, const int maxTrim,
		   Seq *& readSeq1, Seq *& readSeq2, PairMatchMap& pairMap)
{
   SingleMatchMap mmap1, mmap2;

   if (getSingleMatches(readStr1, w, rankTable, maxMinimizer, pmap, inPmap, pvector,
                        minPercentAgreement, minOverlap, readSeq1, mmap1) == 0)
      return 0; // there are no matches of read1

   const String revcomp = stringReverseComplement(readStr2);

   if (getSingleMatches(revcomp, w, rankTable, maxMinimizer, pmap, inPmap, pvector,
                        minPercentAgreement, minOverlap, readSeq2, mmap2,
			&mmap1) == 0)
   {
      // there are no matches of read2, so discard the matches of read1
      freeSingleMatchMap(mmap1);
      delete readSeq1;
      readSeq1 = nullptr;
      return 0;
   }

   // there are matches of read1 and read2

   for (SingleMatchMap::const_iterator mpos1 = mmap1.cbegin();
        mpos1 != mmap1.cend(); ++mpos1) // for each pattern that matches read1
   {
      const int pindex = mpos1->first;

      SingleMatchMap::const_iterator mpos2 = mmap2.find(pindex);
      if (mpos2 == mmap2.cend())
         continue; // there is no match of read2 to this pattern

      const Pattern& pattern = *pvector[pindex];

      PairMatch *best =
         bestPairMatch(pattern, minPercentAgreement, minOverlap, maxInsert, maxTrim,
                       mpos1->second, mpos2->second);
      if (best)
         pairMap.insert(std::make_pair(pindex, best));
   }

   freeSingleMatchMap(mmap1);
   freeSingleMatchMap(mmap2);

   const int numMatches = pairMap.size();
   if (numMatches == 0)
   {
      delete readSeq1; readSeq1 = nullptr;
      delete readSeq2; readSeq2 = nullptr;
   }

   return numMatches;
}

//------------------------------------------------------------------------------------
// selectBestPairMatch() identifies the best overall qualified match

void selectBestPairMatch(PairMatchMap& pairMap, PairMatchMap& bestMap)
{
   int pindex;
   PairMatch *best = nullptr;

   for (PairMatchMap::iterator ppos = pairMap.begin(); ppos != pairMap.end(); ++ppos)
   {
      PairMatch *candidate = ppos->second;

      if (!best || candidate->betterThan(*best))
      {
         pindex = ppos->first;
	 best   = candidate;
      }
   }

   if (best)
      bestMap.insert(std::make_pair(pindex, new PairMatch(*best)));
}

//------------------------------------------------------------------------------------
// freePairMatchMap() de-allocates the matches in the given map

void freePairMatchMap(PairMatchMap& pairMap)
{
   for (PairMatchMap::iterator ppos = pairMap.begin(); ppos != pairMap.end(); ++ppos)
      delete ppos->second;
}
