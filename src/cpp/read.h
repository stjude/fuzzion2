//------------------------------------------------------------------------------------
//
// read.h - module for reading FASTQ and Bam files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef READ_H
#define READ_H

#include "bamread.h"

//------------------------------------------------------------------------------------

class FileReader // abstract class for getting reads from a file
{
public:
   FileReader(const String& inFilename)
      : filename(inFilename), initialName(), initialStr(), initialConsumed(0),
        interleaved(false) { }

   virtual ~FileReader() { }

   virtual void open()  = 0;
   virtual bool readNext(String& name, String& str) = 0;
   virtual void close() = 0;

   void fileOpen();
   bool getNextSingle(String& name,  String& str);
   bool getNextPair  (String& name1, String& str1,
                      String& name2, String& str2);
   void fileClose() { close(); }

   String filename;                      // name of input file
   StringVector initialName, initialStr; // holds the initial reads
   int  initialConsumed;                 // #initial reads consumed
   bool interleaved;                     // true if consecutive reads are mates
};

typedef std::vector<FileReader *> FileVector;

//------------------------------------------------------------------------------------

class FastqFileReader : public FileReader // for reading a FASTQ file
{
public:
   static constexpr int LINEBUF_LEN = 1000000; // accommodate a long read

   FastqFileReader(const String& inFilename)
      : FileReader(inFilename), fileHandle(nullptr), isPipe(false),
        linebuf(new char[LINEBUF_LEN]) { }

   virtual ~FastqFileReader() { delete[] linebuf; }

   virtual void open() override;
   virtual bool readNext(String& name, String& str) override;
   virtual void close() override;

   std::FILE *fileHandle; // is non-null when the file is open
   bool  isPipe;          // is true when reading from a pipe
   char *linebuf;         // line buffer, holds a single line read from a FASTQ file
};

//------------------------------------------------------------------------------------

class BamFileReader : public FileReader // for reading a Bam file
{
public:
   BamFileReader(const String& inFilename)
      : FileReader(inFilename), bamReader(), bamRead() { }

   virtual ~BamFileReader() { }

   virtual void open() override { bamReader.open(filename); }
   virtual bool readNext(String& name, String& str) override;
   virtual void close() override { bamReader.close(); }

   BamReader bamReader; // used to read the Bam file
   BamRead   bamRead;   // buffer to hold a Bam read
};

//------------------------------------------------------------------------------------

class SingleReader // for getting single reads from one or more files
{
public:
   SingleReader(const FileVector& inFileVector)
      : fileVector(inFileVector), current(0) { }

   virtual ~SingleReader();

   bool getNext(String& name, String& str);
   void close();

   FileVector fileVector;
   int current; // index into fileVector of current file
};

//------------------------------------------------------------------------------------

class PairReader // for getting read pairs from one or more files
{
public:
   PairReader(const FileVector& inFileVector1, const FileVector& inFileVector2)
      : fileVector1(inFileVector1), fileVector2(inFileVector2), current(0) { }

   virtual ~PairReader();

   bool getNext(String& name1, String& str1,
                String& name2, String& str2);
   void close();

   FileVector fileVector1, fileVector2;
   int current; // index into fileVector1 and fileVector2 of current file
};

//------------------------------------------------------------------------------------

SingleReader *createSingleReader(const String& fastq1, const String& fastq2,
                                 const String& ifastq, const String& ubam,
                                 const StringVector& filename);

PairReader *createPairReader(const String& fastq1, const String& fastq2,
                             const String& ifastq, const String& ubam,
                             const StringVector& filename);

#endif
