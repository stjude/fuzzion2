//------------------------------------------------------------------------------------
//
// refgen.h - module for reading reference genome files in 2-bit format
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef REFGEN_H
#define REFGEN_H

#include "bin.h"
#include <vector>

//------------------------------------------------------------------------------------

class RefGenSeq // represents a sequence of a reference
{
public:
   RefGenSeq(int beginPos, int endPos, char *sequence)
      : begin(beginPos), end(endPos), seq(sequence) { }

   virtual ~RefGenSeq() { delete[] seq; }

   inline char getBase(int position) const
   { return (position >= begin && position <= end ? seq[position - begin] : 'N'); }

   int   begin; // first position of sequence, 1-based, inclusive
   int   end;   // last  position of sequence, 1-based, inclusive

   char *seq;   // sequence of (end - begin + 1) bases
};

//------------------------------------------------------------------------------------

class RefGenReader : public BinReader
{
public:
   RefGenReader(int bufferSize=DEFAULT_BINARY_BUFFER_SIZE)
      : BinReader(bufferSize), numref(0), refName(), refOffset() { }

   virtual ~RefGenReader() { }

   void open(const std::string& refGenFilename);

   RefGenSeq *getRefGenSeq(const std::string& selectedRefName,
                           int beginPos, int endPos);

   void close();

   int numref;                         // number of references
   std::vector<std::string> refName;   // name of each reference
   std::vector<uint32_t>    refOffset; // byte offset of each reference in 2-bit file

protected:
   uint32_t readUint32();
};

//------------------------------------------------------------------------------------
#endif
