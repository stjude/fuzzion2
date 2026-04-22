//------------------------------------------------------------------------------------
//
// read.cpp - module for reading FASTQ and Bam files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "read.h"

const String UNPAIRED_MESSAGE = "; use -single=1 for unpaired reads";
const int INITIAL_READS = 100; // #initial reads to get when opening a file

//------------------------------------------------------------------------------------
// namesMatch() returns true if the given read names match, indicating that the reads
// are mates of a pair

static bool namesMatch(const String& name1, const String& name2)
{
   const int len = name1.length();
   if (len != name2.length())
      return false;

   // see if one name ends with "1" and the other ends with "2"
   const int last = len - 1;

   if (last > 0 && (name1[last] == '1' && name2[last] == '2' ||
                    name1[last] == '2' && name2[last] == '1'))
      return (std::strncmp(name1.c_str(), name2.c_str(), last) == 0);

   return (name1 == name2);
}

//------------------------------------------------------------------------------------
// FileReader::fileOpen() opens the file and gets the initial reads; from these, it is
// determined whether consecutive reads are mates

void FileReader::fileOpen()
{
   open();

   String name, str;

   for (int i = 0; i < INITIAL_READS && readNext(name, str); i++)
   {
      initialName.push_back(name);
      initialStr .push_back(str);
   }

   const int numInitial = initialName.size();
   if (numInitial == 0)
      throw std::runtime_error("no reads in " + filename);

   if (INITIAL_READS < 2 ||
       numInitial < INITIAL_READS && numInitial % 2 == 1 /* odd #reads in file */)
      interleaved = false;
   else // check for matching read names
   {
      int i = 0;
      while (i + 1 < numInitial && namesMatch(initialName[i], initialName[i + 1]))
         i += 2;

      interleaved = (i + 1 >= numInitial);
   }
}

//------------------------------------------------------------------------------------
// FileReader::getNextSingle() gets the next read from the file and returns true if
// successful; false is returned if end-of-file is reached

bool FileReader::getNextSingle(String& name, String& str)
{
   const int numInitial = initialName.size();

   if (initialConsumed < numInitial)
   {
      // get one of the initial reads
      name = initialName[initialConsumed];
      str  = initialStr [initialConsumed++];
      return true;
   }

   return (numInitial == INITIAL_READS && readNext(name, str));
}

//------------------------------------------------------------------------------------
// FileReader::getNextPair() gets the next read pair from the file and returns true if
// successful; false is returned if end-of-file is reached

bool FileReader::getNextPair(String& name1, String& str1,
                             String& name2, String& str2)
{
   if (!getNextSingle(name1, str1))
      return false; // reached end-of-file

   if (!getNextSingle(name2, str2))
      throw std::runtime_error("odd #reads in " + filename + UNPAIRED_MESSAGE);

   // got two reads; we expect them to be mates of a pair
   if (!namesMatch(name1, name2))
      throw std::runtime_error("mismatched read names " + name1 + " and " + name2 +
                               " in " + filename + UNPAIRED_MESSAGE);

   return true;
}

//------------------------------------------------------------------------------------
// FastqFileReader::open() opens a FASTQ file; if the file is gzipped, it will be
// uncompressed

void FastqFileReader::open()
{
   if (fileHandle) // the file is already open
      close();

   const int len = filename.length();

   if (len > 3 && filename.substr(len - 3) == ".gz")  // the file is gzipped
   {
      const String command = "gunzip -c " + filename; // we will uncompress it
      fileHandle = popen(command.c_str(), "r");
      isPipe     = true;
   }
   else // the file is uncompressed
   {
      fileHandle = std::fopen(filename.c_str(), "r");
      isPipe     = false;
   }

   if (!fileHandle)
      throw std::runtime_error("unable to open " + filename);
}

//------------------------------------------------------------------------------------
// FastqFileReader::readNext() gets the next read from a FASTQ file and returns true
// if successful; false is returned if end-of-file is reached

bool FastqFileReader::readNext(String& name, String& str)
{
   if (!fileHandle)
      throw std::runtime_error("attempt to read from unopened FASTQ file " +
                               filename);

   if (!std::fgets(linebuf, LINEBUF_LEN, fileHandle))
      return false; // reached end-of-file

   if (linebuf[0] == '@')
   {
      unsigned char ch;
      int i = 1;
      while ((ch = linebuf[i]) && !std::isspace(ch))
         i++;

      linebuf[i] = 0;
      name = &linebuf[1];

      if (std::fgets(linebuf, LINEBUF_LEN, fileHandle))
      {
         i = 0;
	 while ((ch = linebuf[i]) && !std::isspace(ch))
            i++;

	 linebuf[i] = 0;
	 str = linebuf;

	 if (std::fgets(linebuf, LINEBUF_LEN, fileHandle) && linebuf[0] == '+' &&
             std::fgets(linebuf, LINEBUF_LEN, fileHandle))
            return true;
      }
   }

   throw std::runtime_error("unrecognized FASTQ file format in " + filename);
}

//------------------------------------------------------------------------------------
// FastqFileReader::close() closes a FASTQ file

void FastqFileReader::close()
{
   if (!fileHandle) // the file is not open
      return;

   if (isPipe)
   {
      pclose(fileHandle);
      isPipe = false;
   }
   else
      std::fclose(fileHandle);

   fileHandle = nullptr;
}

//------------------------------------------------------------------------------------
// BamFileReader::readNext() gets the next read from a Bam file and returns true if
// successful; false is returned if end-of-file is reached

bool BamFileReader::readNext(String& name, String& str)
{
   while (true)
      if (bamReader.getNext(bamRead))
      {
         if (!bamRead.isUnmapped() && bamRead.isSecondary())
            continue; // skip secondary alignment

	 name = bamRead.constName();

	 const String *s = bamRead.sequence();
	 str = *s;
	 delete s;

	 return true;
      }
      else
         return false; // reached end-of-file
}

//------------------------------------------------------------------------------------
// SingleReader::~SingleReader() de-allocates its FileReader objects

SingleReader::~SingleReader()
{
   const int numFiles = fileVector.size();

   for (int i = 0; i < numFiles; i++)
      delete fileVector[i];
}

//------------------------------------------------------------------------------------
// SingleReader::getNext() gets the next read and returns true if successful; false is
// returned when all files have been fully read

bool SingleReader::getNext(String& name, String& str)
{
   const int numFiles = fileVector.size();

   if (current < numFiles)
      if (fileVector[current]->getNextSingle(name, str))
         return true;
      else // reached the end of the current file
         if (++current < numFiles) // there is anothe file; get its first read
            return fileVector[current]->getNextSingle(name, str);

   return false; // reached the end of all files
}

//------------------------------------------------------------------------------------
// SingleReader::close() closes each of its files

void SingleReader::close()
{
   const int numFiles = fileVector.size();

   for (int i = 0; i < numFiles; i++)
      fileVector[i]->fileClose();
}

//------------------------------------------------------------------------------------
// PairReader::~PairReader() de-allocates its FileReader objects

PairReader::~PairReader()
{
   const int numFiles1 = fileVector1.size();
   const int numFiles2 = fileVector2.size();

   for (int i = 0; i < numFiles1; i++)
      delete fileVector1[i];

   for (int i = 0; i < numFiles2; i++)
      delete fileVector2[i];
}

//------------------------------------------------------------------------------------
// getPair() gets the next read pair from a single file (if fileReader2 is nullptr) or
// from a pair of files (if fileReader2 is non-null); true is returned if successful;
// false is returned at end-of-file

static bool getPair(FileReader *fileReader1, FileReader *fileReader2,
                    String& name1, String& str1,
		    String& name2, String& str2)
{
   if (!fileReader2)
      return fileReader1->getNextPair(name1, str1, name2, str2);

   const bool gotRead1 = fileReader1->getNextSingle(name1, str1);
   const bool gotRead2 = fileReader2->getNextSingle(name2, str2);

   if (!gotRead1 && !gotRead2)
      return false; // reached the end of both files

   if (gotRead1 && gotRead2)
      if (namesMatch(name1, name2))
         return true;
      else
         throw std::runtime_error("mismatched read names " +
                                  name1 + " in " + fileReader1->filename + " and " +
				  name2 + " in " + fileReader2->filename +
				  UNPAIRED_MESSAGE);

   // we reached the end of one file but not the other
   throw std::runtime_error(fileReader1->filename + " and " + fileReader2->filename +
                            " contain different #reads" + UNPAIRED_MESSAGE);
}

//------------------------------------------------------------------------------------
// PairReader::getNext() gets the next read pair and returns true if successful; false
// is returned when all files have been fully read

bool PairReader::getNext(String& name1, String& str1,
                         String& name2, String& str2)
{
   const int numFiles1 = fileVector1.size();
   const int numFiles2 = fileVector2.size();

   if (current < numFiles1 && current < numFiles2)
      if (getPair(fileVector1[current], fileVector2[current],
                  name1, str1, name2, str2))
         return true;
      else // reached the end of the current files
         if (++current < numFiles1 && current < numFiles2)
            return getPair(fileVector1[current], fileVector2[current],
                           name1, str1, name2, str2);

   return false; // reached the end of all files
}

//------------------------------------------------------------------------------------
// PairReader::close() closes each of its files

void PairReader::close()
{
   const int numFiles1 = fileVector1.size();
   const int numFiles2 = fileVector2.size();

   for (int i = 0; i < numFiles1; i++)
      fileVector1[i]->fileClose();

   for (int i = 0; i < numFiles2; i++)
      if (fileVector2[i])
         fileVector2[i]->fileClose();
}

//------------------------------------------------------------------------------------
// createSingleReader() returns a SingleReader that is set up to get single reads from
// the named files; it is the caller's obligation to de-allocate it

SingleReader *createSingleReader(const String& fastq1, const String& fastq2,
                                 const String& ifastq, const String& ubam,
				 const StringVector& filename)
{
   FileVector  fileVector;
   FileReader *fileReader;

   if (fastq1 != "")
   {
      fileReader = new FastqFileReader(fastq1);
      fileReader->fileOpen();
      fileVector.push_back(fileReader);
   }

   if (fastq2 != "")
   {
      fileReader = new FastqFileReader(fastq2);
      fileReader->fileOpen();
      fileVector.push_back(fileReader);
   }

   if (ifastq != "")
   {
      fileReader = new FastqFileReader(ifastq);
      fileReader->fileOpen();
      fileVector.push_back(fileReader);
   }

   if (ubam != "")
   {
      fileReader = new BamFileReader(ubam);
      fileReader->fileOpen();
      fileVector.push_back(fileReader);
   }

   const int numFiles = filename.size();

   for (int i = 0; i < numFiles; i++)
   {
      fileReader = new BamFileReader(filename[i]);

      try
      {
         fileReader->fileOpen();
      }
      catch (const std::runtime_error& error)
      {
         // must not be a Bam file
	 delete fileReader;
	 fileReader = nullptr;
      }

      if (!fileReader)
      {
         fileReader = new FastqFileReader(filename[i]);
	 fileReader->fileOpen();
      }

      fileVector.push_back(fileReader);
   }

   if (fileVector.size() == 0)
      throw std::runtime_error("no FASTQ or Bam file specified");

   return new SingleReader(fileVector);
}

//------------------------------------------------------------------------------------
// pairable() returns true if the initial reads from two files are mates

static bool pairable(const FileReader *fileReader1, const FileReader *fileReader2)
{
   const int numInitial = fileReader1->initialName.size();
   if (numInitial != fileReader2->initialName.size())
      return false; // the number of initial reads differs

   for (int i = 0; i < numInitial; i++)
      if (!namesMatch(fileReader1->initialName[i], fileReader2->initialName[i]))
         return false; // found read names that do not match

   return true;
}

//------------------------------------------------------------------------------------
// createPairReader() returns a PairReader that is set up to get read pairs from the
// named files; it is the caller's obligation to de-allocate it

PairReader *createPairReader(const String& fastq1, const String& fastq2,
                             const String& ifastq, const String& ubam,
                             const StringVector& filename)
{
   FileVector  fileVector1,  fileVector2;
   FileReader *fileReader1, *fileReader2;

   if (fastq1 != "" && fastq2 != "")
   {
      fileReader1 = new FastqFileReader(fastq1);
      fileReader2 = new FastqFileReader(fastq2);

      fileReader1->fileOpen();
      fileReader2->fileOpen();

      if (pairable(fileReader1, fileReader2))
      {
         fileVector1.push_back(fileReader1);
         fileVector2.push_back(fileReader2);
      }
      else
         throw std::runtime_error(fastq1 + " and " + fastq2 + " are incompatible");
   }
   else if (fastq1 != "" || fastq2 != "")
      throw std::runtime_error("unspecified FASTQ mate file" + UNPAIRED_MESSAGE);

   if (ifastq != "")
   {
      fileReader1 = new FastqFileReader(ifastq);
      fileReader1->fileOpen();

      if (fileReader1->interleaved)
      {
         fileVector1.push_back(fileReader1);
         fileVector2.push_back(nullptr);
      }
      else
         throw std::runtime_error("mates are not consecutive in " + ifastq +
                                  UNPAIRED_MESSAGE);
   }

   if (ubam != "")
   {
      fileReader1 = new BamFileReader(ubam);
      fileReader1->fileOpen();

      if (fileReader1->interleaved)
      {
         fileVector1.push_back(fileReader1);
         fileVector2.push_back(nullptr);
      }
      else
         throw std::runtime_error("mates are not consecutive in " + ubam +
                                  UNPAIRED_MESSAGE);
   }

   FileVector unpaired;
   const int numFiles = filename.size();

   for (int i = 0; i < numFiles; i++)
   {
      fileReader1 = new BamFileReader(filename[i]);

      try
      {
         fileReader1->fileOpen();
      }
      catch (const std::runtime_error& error)
      {
         // must not be a Bam file
	 delete fileReader1;
	 fileReader1 = nullptr;
      }

      if (fileReader1) // this is a Bam file
         if (fileReader1->interleaved)
         {
            fileVector1.push_back(fileReader1);
            fileVector2.push_back(nullptr);
         }
         else
            throw std::runtime_error("mates are not consecutive in " + filename[i] +
                                     UNPAIRED_MESSAGE);
      else // this is not a Bam file
      {
         fileReader1 = new FastqFileReader(filename[i]);
         fileReader1->fileOpen();

         if (fileReader1->interleaved)
         {
            fileVector1.push_back(fileReader1);
            fileVector2.push_back(nullptr);
         }
	 else
            unpaired.push_back(fileReader1);
      }

   }

   const int numUnpaired = unpaired.size();

   if (numUnpaired > 0)
   {
      BoolVector paired(numUnpaired, false);

      for (int i = 0; i < numUnpaired; i++)
         if (!paired[i])
	 {
            for (int j = i + 1; j < numUnpaired; j++)
               if (!paired[j] && pairable(unpaired[i], unpaired[j]))
	       {
                  paired[i] = true;
		  paired[j] = true;

		  fileVector1.push_back(unpaired[i]);
		  fileVector2.push_back(unpaired[j]);
		  break;
	       }

	    if (!paired[i])
               throw std::runtime_error("no FASTQ mate file for " +
                                        unpaired[i]->filename + UNPAIRED_MESSAGE);
	 }
   }

   if (fileVector1.size() == 0)
      throw std::runtime_error("no FASTQ or Bam file specified");

   return new PairReader(fileVector1, fileVector2);
}
