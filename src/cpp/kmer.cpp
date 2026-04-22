//------------------------------------------------------------------------------------
//
// kmer.cpp - module for representing and processing k-mers in 2-bit format
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "kmer.h"

//------------------------------------------------------------------------------------
// kmerToString() converts the given k-mer to a string

String kmerToString(const KmerLength k, Kmer kmer)
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

Kmer stringToKmer(const String& s)
{
   if (s.length() > MAX_KMER_LENGTH)
      throw std::runtime_error("unsupported k-mer length");

   const KmerLength k = s.length();
   Kmer kmer = 0;

   for (int i = 0; i < k; i++)
   {
      const Base base = charToBase(s[i]);

      if (base < NUM_BASES)
         kmer = (kmer << 2) | base;
      else
         throw std::runtime_error("cannot convert " + s);
   }

   return kmer;
}

//------------------------------------------------------------------------------------
// kmerReverseComplement() returns the reverse complement of the given k-mer

Kmer kmerReverseComplement(const KmerLength k, Kmer kmer)
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
// stringReverseComplement() returns the reverse complement of the given string;
// case (upper or lower) is retained

String stringReverseComplement(const String& s)
{
   const int A_to_T = 'T' - 'A';
   const int C_to_G = 'G' - 'C';

   String revcomp = s;
   const int slen = s.length();

   for (int i = 0; i < slen; i++)
   {
      char ch = s[slen - i - 1];

      switch (ch)
      {
         case 'A': case 'a': ch += A_to_T; break;
         case 'C': case 'c': ch += C_to_G; break;
         case 'G': case 'g': ch -= C_to_G; break;
         case 'T': case 't': ch -= A_to_T; break;
         default : break; // leave character unchanged
      }

      revcomp[i] = ch;
   }

   return revcomp;
}

//------------------------------------------------------------------------------------
// KmerFinder::KmerFinder() checks for a valid k-mer length and initializes data
// members

KmerFinder::KmerFinder(const char *sequence, const int sequenceLen,
                       const KmerLength kmerLen)
   : seq(sequence), seqlen(sequenceLen), k(kmerLen), mask(numKmers(kmerLen) - 1)
{
   if (k < 1 || k > MAX_KMER_LENGTH)
      throw std::runtime_error("unsupported k-mer length");
}

//------------------------------------------------------------------------------------
// KmerFinder::find() searches for k-mers in the sequence; each k-mer that is found is
// reported along with its starting index in the sequence; the search terminates when
// the sequence has been exhausted or when reportKmer() returns false; this function
// consumes a large fraction of the running time of fuzzion2 and so an effort has been
// made to optimize it

void KmerFinder::find()
{
   Kmer kmer = 0;      // the k-mer under construction
   KmerLength len = 0; // length of the k-mer under construction

   for (int i = 0; i < seqlen; i++)
      switch (seq[i])
      {
         case 'A': case 'a':
            if (len == k)
	    {
               if (!reportKmer(kmer = (kmer << 2) & mask, i - k + 1))
                  return;
	    }
	    else if (++len == k)
	    {
               if (!reportKmer(kmer, i - k + 1))
                  return;
	    }
	    break;

	 case 'C': case 'c':
            if (len == k)
	    {
               if (!reportKmer(kmer = ((kmer << 2) & mask) | BASE_C, i - k + 1))
                  return;
	    }
	    else if (++len == k)
	    {
               if (!reportKmer(kmer |= BASE_C, i - k + 1))
                  return;
	    }
	    else
               kmer |= BASE_C << ((k - len) << 1);
	    break;

	 case 'G': case 'g':
            if (len == k)
	    {
               if (!reportKmer(kmer = ((kmer << 2) & mask) | BASE_G, i - k + 1))
                  return;
	    }
	    else if (++len == k)
	    {
               if (!reportKmer(kmer |= BASE_G, i - k + 1))
                  return;
	    }
	    else
               kmer |= BASE_G << ((k - len) << 1);
	    break;

	 case 'T': case 't':
            if (len == k)
	    {
               if (!reportKmer(kmer = ((kmer << 2) & mask) | BASE_T, i - k + 1))
                  return;
	    }
	    else if (++len == k)
	    {
               if (!reportKmer(kmer |= BASE_T, i - k + 1))
                  return;
	    }
	    else
               kmer |= BASE_T << ((k - len) << 1);
	    break;

	 default: // some other base, such as 'N'
            if (len > 0)
	    {
               kmer = 0;
	       len  = 0;
	    }
      }
}
