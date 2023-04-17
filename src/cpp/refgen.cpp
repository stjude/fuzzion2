//------------------------------------------------------------------------------------
//
// refgen.cpp - module for reading reference genome files in 2-bit format
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "refgen.h"
#include "util.h"
#include <stdexcept>

const uint32_t TWO_BIT_SIGNATURE_NOSWAP = 0x1A412743;
const uint32_t TWO_BIT_SIGNATURE_SWAP   = 0x4327411A;

//------------------------------------------------------------------------------------
// RefGenReader::readUint32() reads and returns an unsigned 32-bit integer and throws
// an exception if end-of-file was reached

uint32_t RefGenReader::readUint32()
{
   uint32_t value;

   if (BinReader::readUint32(value))
      return value;

   throw std::runtime_error("invalid format in " + filename);
}

//------------------------------------------------------------------------------------
// RefGenReader::open() opens the named file containing a reference genome in 2-bit
// format; reference names and byte offsets are read from the file header

void RefGenReader::open(const std::string& refGenFilename)
{
   BinReader::open(refGenFilename);

   swap = false;

   uint32_t signature = readUint32();

   if (signature == TWO_BIT_SIGNATURE_SWAP)
      swap = true;
   else if (signature != TWO_BIT_SIGNATURE_NOSWAP)
      throw std::runtime_error(refGenFilename + " is not a 2-bit file");

   readUint32();          // skip version
   numref = readUint32(); // read number of references
   readUint32();          // skip reserved

   uint8_t namelen;
   char name[256];

   for (int i = 0; i < numref; i++)
      if (readUint8(namelen) && readBuffer(name, namelen))
      {
         name[namelen] = 0; // append null terminating byte
	 refName.push_back(name);
	 refOffset.push_back(readUint32());
      }
      else
         throw std::runtime_error("invalid format in " + refGenFilename);
}

//------------------------------------------------------------------------------------
// RefGenReader::getRefGenSeq() reads a sequence of the named reference and returns it
// in a newly allocated object; the caller is obligated to de-allocate it; the
// specified positions are 1-based and inclusive; to get all positions of the named
// reference, pass 1 for beginPos and a really large value such as 1000000000 (one
// billion) for endPos; in this case, the end position will be set to the reference
// length (i.e., the number of positions in the reference)

RefGenSeq *RefGenReader::getRefGenSeq(const std::string& selectedRefName,
                                      int beginPos, int endPos)
{
   int i = 0;
   while (i < numref && refName[i] != selectedRefName)
      i++;

   if (i == numref)
      throw std::runtime_error("unrecognized reference name \"" + selectedRefName +
                               "\"");

   seek(refOffset[i]);

   uint32_t reflen = readUint32();

   if (endPos > reflen)
      endPos = reflen;

   if (beginPos < 1 || beginPos > endPos)
      throw std::runtime_error("invalid position");

   char *sequence = new char[endPos - beginPos + 1];

   uint32_t nBlockCount = readUint32();

   IntVector nstart(nBlockCount); // vector of N-block first positions,
                                         // 1-based, inclusive

   for (uint32_t i = 0; i < nBlockCount; i++)
      nstart[i] = readUint32() + 1;      // make it 1-based

   IntVector nstop(nBlockCount);  // vector of N-block last positions,
                                         // 1-based, inclusive

   for (uint32_t i = 0; i < nBlockCount; i++)
      nstop[i] = nstart[i] + readUint32() - 1;

   uint32_t maskBlockCount = readUint32();

   uint32_t dnaOffset = refOffset[i] +
      2 * sizeof(uint32_t) * (nBlockCount + maskBlockCount + 2);

   seek(dnaOffset + ((beginPos - 1) >> 2)); // jump to first DNA byte

   // read the 2-bit sequence data and convert it to characters

   const char BASE[4] = { 'T', 'C', 'A', 'G' };

   uint8_t byte;
   bool readByte = true;

   for (int pos = beginPos; pos <= endPos; pos++)
   {
      if (readByte && !readUint8(byte))
         throw std::runtime_error("truncated 2-bit file " + filename);

      uint8_t shift = 2 * (3 - ((pos - 1) & 3));
      uint8_t index = (byte >> shift) & 3;

      sequence[pos - beginPos] = BASE[index];

      readByte = (shift == 0);
   }

   // now insert N's

   for (uint32_t i = 0; i < nBlockCount; i++)
      if (nstart[i] <= endPos && nstop[i] >= beginPos)
      {
         int start = std::max(nstart[i], beginPos);
	 int stop  = std::min(nstop[i],  endPos);

	 for (int pos = start; pos <= stop; pos++)
            sequence[pos - beginPos] = 'N';
      }

   return new RefGenSeq(beginPos, endPos, sequence);
}

//------------------------------------------------------------------------------------
// RefGenReader::close() closes the 2-bit file

void RefGenReader::close()
{
   numref = 0;
   refName.clear();
   refOffset.clear();

   BinReader::close();
}
