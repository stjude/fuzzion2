//------------------------------------------------------------------------------------
//
// refgen.h - module for reading reference genome files in 2-bit format
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef REFGEN_H
#define REFGEN_H

#include "bin.h"

//------------------------------------------------------------------------------------

class RefGenSeq // represents a sequence of a reference
{
public:
   RefGenSeq(const int beginPos, const int endPos, const char *sequence)
      : begin(beginPos), end(endPos), seq(sequence) { }

   virtual ~RefGenSeq() { delete[] seq; }

   inline char getBase(const int position) const
   { return (position >= begin && position <= end ? seq[position - begin] : 'N'); }

   const int   begin; // first position of sequence, 1-based, inclusive
   const int   end;   // last  position of sequence, 1-based, inclusive
   const char *seq;   // sequence of (end - begin + 1) bases
};

//------------------------------------------------------------------------------------

class RefGenReader : public BinReader
{
public:
   RefGenReader(const int bufferSize=DEFAULT_BINARY_BUFFER_SIZE)
      : BinReader(bufferSize), numref(0), refName(), refOffset() { }

   virtual ~RefGenReader() { }

   void open(const String& refGenFilename);

   RefGenSeq *getRefGenSeq(const String& selectedRefName, int beginPos, int endPos);

   void close();

   int numref;                      // number of references
   StringVector refName;            // name of each reference
   std::vector<uint32_t> refOffset; // byte offset of each reference in 2-bit file

protected:
   uint32_t readUint32();
};

//------------------------------------------------------------------------------------
#endif
