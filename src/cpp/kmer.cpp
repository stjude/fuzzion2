//------------------------------------------------------------------------------------
//
// kmer.cpp - module for representing and processing k-mers in 2-bit format
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "kmer.h"
#include <stdexcept>

//------------------------------------------------------------------------------------
// kmerToString() converts the given k-mer to a string

std::string kmerToString(KmerLength k, Kmer kmer)
{
   if (k > MAX_KMER_LENGTH)
      throw std::runtime_error("unsupported k-mer length");

   char buffer[MAX_KMER_LENGTH + 1];
   buffer[MAX_KMER_LENGTH] = 0; // null terminating byte

   for (int i = 1; i <= k; i++)
   {
      buffer[MAX_KMER_LENGTH - i] = baseToChar(kmer & 3);
      kmer >>= 2;
   }

   return &buffer[MAX_KMER_LENGTH - k];
}

//------------------------------------------------------------------------------------
// stringToKmer() converts the given string to a k-mer

Kmer stringToKmer(const std::string& s)
{
   if (s.length() > MAX_KMER_LENGTH)
      throw std::runtime_error("unsupported k-mer length");

   KmerLength k = s.length();

   Kmer kmer = 0;

   for (int i = 0; i < k; i++)
   {
      Base base = charToBase(s[i]);

      if (base < NUM_BASES)
         kmer = (kmer << 2) | base;
      else
         throw std::runtime_error("cannot convert " + s);
   }

   return kmer;
}

//------------------------------------------------------------------------------------
// kmerReverseComplement() returns the reverse complement of the given k-mer

Kmer kmerReverseComplement(KmerLength k, Kmer kmer)
{
   if (k > MAX_KMER_LENGTH)
      throw std::runtime_error("unsupported k-mer length");

   Kmer revcomp = 0;

   for (int i = 1; i <= k; i++)
   {
      revcomp = (revcomp << 2) | baseComplement(kmer & 3);
      kmer >>= 2;
   }

   return revcomp;
}

//------------------------------------------------------------------------------------
// charReverseComplement() returns the reverse complement of the given array of char
// in a newly allocated array of char; it is the caller's obligation to de-allocate it

char *charReverseComplement(const char *s, int slen)
{
   char *revcomp = new char[slen];

   for (int i = 0; i < slen; i++)
   {
      char ch    = s[slen - i - 1];
      Base base  = charToBase(ch);

      revcomp[i] = (base < NUM_BASES ? BASE_CHAR[3 - base] : ch);
   }

   return revcomp;
}

//------------------------------------------------------------------------------------
// stringReverseComplement() returns the reverse complement of the given string

std::string stringReverseComplement(const std::string& s)
{
   std::string revcomp = s;

   int slen = s.length();

   for (int i = 0; i < slen; i++)
   {
      char ch    = s[slen - i - 1];
      Base base  = charToBase(ch);

      revcomp[i] = (base < NUM_BASES ? BASE_CHAR[3 - base] : ch);
   }

   return revcomp;
}

//------------------------------------------------------------------------------------
// KmerFinder::KmerFinder() checks for a valid k-mer length and initializes data
// members

KmerFinder::KmerFinder(const char *sequence, int sequenceLen, KmerLength kmerLen)
   : seq(sequence), seqlen(sequenceLen), k(kmerLen)
{
   if (k < 1 || k > MAX_KMER_LENGTH)
      throw std::runtime_error("unsupported k-mer length");

   mask = numKmers(k) - 1;
}

//------------------------------------------------------------------------------------
// KmerFinder::find() searches for k-mers in the sequence; each k-mer that is found is
// reported along with its starting index in the sequence; the search terminates when
// the sequence has been exhausted or when report() returns false

void KmerFinder::find()
{
   Kmer kmer = 0;      // the k-mer under construction
   KmerLength len = 0; // length of the k-mer under construction

   for (int i = 0; i < seqlen; i++)
   {
      Base base = charToBase(seq[i]);

      if (base < NUM_BASES)
      {
         if (len < k)
	 {
            len++;
	    kmer |= base << 2 * (k - len);
	 }
	 else // len == k
            kmer = ((kmer << 2) & mask) | base;

	 if (len == k && !reportKmer(kmer, i - k + 1))
            break;
      }
      else // encountered an unrecognized base such as 'N'
         if (len > 0)
	 {
            kmer = 0;
	    len  = 0;
	 }
   }
}
