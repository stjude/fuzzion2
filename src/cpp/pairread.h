//------------------------------------------------------------------------------------
//
// pairread.h - module relating to paired reads
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef PAIRREAD_H
#define PAIRREAD_H

#include <string>
#include <vector>

//------------------------------------------------------------------------------------

class PairReader // abstract class for getting paired reads
{
public:
   PairReader() { }

   virtual ~PairReader() { }

   virtual void open() = 0;

   virtual bool getNextPair(std::string& name1, std::string& sequence1,
                            std::string& name2, std::string& sequence2) = 0;

   virtual void close() = 0;
};

typedef std::vector<PairReader *> PairReaderVector;

//------------------------------------------------------------------------------------

class InputReader : public PairReader // gets paired reads from a series of files
{
public:
   InputReader(const PairReaderVector& inReaderVector)
      : readerVector(inReaderVector), current(-1) { }

   virtual ~InputReader() { }

   void open();

   bool getNextPair(std::string& name1, std::string& sequence1,
                    std::string& name2, std::string& sequence2);

   void close();

   PairReaderVector readerVector;
   int current; // index into readerVector of open reader, or -1 if none
};

//------------------------------------------------------------------------------------

bool namesMatch(const std::string& name1, const std::string& name2);

#endif
