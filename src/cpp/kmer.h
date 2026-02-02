//------------------------------------------------------------------------------------
//
// kmer.h - module for representing and processing k-mers in 2-bit format
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef KMER_H
#define KMER_H

#include "util.h"

typedef uint8_t Base;

const Base BASE_A     = 0;   // 2-bit representation of adenine
const Base BASE_C     = 1;   // 2-bit representation of cytosine
const Base BASE_G     = 2;   // 2-bit representation of guanine
const Base BASE_T     = 3;   // 2-bit representation of thymine
const Base NUM_BASES  = 4;   // number of defined bases
const Base BASE_OTHER = 4;   // represents an undefined base, such as N

const char BASE_CHAR[NUM_BASES] = { 'A', 'C', 'G', 'T' };
const char UNDEFINED_BASE = 'N';

typedef uint32_t Kmer;       // holds a k-mer in 2-bit format
typedef uint8_t  KmerLength; // holds the length of a k-mer, i.e., the value of k

const KmerLength MAX_KMER_LENGTH = 15; // maximum supported length of a k-mer

inline Kmer numKmers(const KmerLength k)
{ return (k <= MAX_KMER_LENGTH ? 1 << (2 * k) : 0); }

//------------------------------------------------------------------------------------

inline char baseToChar(const Base base)
{ return (base < NUM_BASES ? BASE_CHAR[base] : UNDEFINED_BASE); }

inline Base charToBase(const char ch)
{
   switch (ch)
   {
      case 'A': case 'a': return BASE_A;
      case 'C': case 'c': return BASE_C;
      case 'G': case 'g': return BASE_G;
      case 'T': case 't': return BASE_T;
      default : return BASE_OTHER;
   }
}

//------------------------------------------------------------------------------------

inline char upperACGT(const char ch)
{
   switch (ch)
   {
      case 'a': case 'c': case 'g': case 't':
         return ch + 'A' - 'a';
      default:
         return ch;
   }
}

//------------------------------------------------------------------------------------

inline Base baseComplement(const Base base)
{ return (base < NUM_BASES ? 3 - base : BASE_OTHER); }

//------------------------------------------------------------------------------------

String kmerToString(KmerLength k, Kmer kmer);

Kmer stringToKmer(const String& s);

//------------------------------------------------------------------------------------

Kmer kmerReverseComplement(KmerLength k, Kmer kmer);

String stringReverseComplement(const String& s);

//------------------------------------------------------------------------------------

class KmerFinder // abstract class for finding k-mers in a sequence
{
public:
   KmerFinder(const char *sequence, int sequenceLen, KmerLength kmerLen);

   virtual ~KmerFinder() { }

   virtual void find(); // finds and reports k-mers in the given sequence

   // override this function to do whatever you want with each k-mer found;
   // return true to continue, false to discontinue the search
   virtual bool reportKmer(Kmer kmer, int startIndex) = 0;

   const char *seq;    // sequence to search
   const int seqlen;   // length of the sequence
   const KmerLength k; // length of k-mers to find
   const Kmer mask;    // mask used internally
};

//------------------------------------------------------------------------------------
#endif
