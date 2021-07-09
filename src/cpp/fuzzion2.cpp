//------------------------------------------------------------------------------------
//
// fuzzion2.cpp - this program is a "fuzzy fusion finder"; it does fuzzy matching of
//                read pairs to fusion patterns
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "fastq.h"
#include "match.h"
#include "ubam.h"
#include <iomanip>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <thread>

const std::string VERSION_NAME   = "fuzzion2 v1.0.1";
const std::string COPYRIGHT      =
   "copyright 2021 St. Jude Children's Research Hospital";

const int    THREAD_BATCH_SIZE   = 100000; // number of read pairs in a full batch

const double DEFAULT_MAX_RANK    = 95;  // default max rank percentile of minimizers
const double DEFAULT_MIN_BASES   = 90;  // default min percentile of matching bases
const int    DEFAULT_MAX_INSERT  = 500; // default max insert size in bases
const int    DEFAULT_MIN_MINS    = 3;   // default min number of matching minimizers
const int    DEFAULT_MIN_OVERLAP = 5;   // default min length of overlap in #bases
const int    DEFAULT_SHOW        = 1;   // default setting of -show option
const int    DEFAULT_THREADS     = 8;   // default number of threads
const int    DEFAULT_WINDOW_LEN  = 5;   // default length of windows in #bases

double maxRank    = DEFAULT_MAX_RANK;
double minBases   = DEFAULT_MIN_BASES;
int    maxInsert  = DEFAULT_MAX_INSERT;
int    minMins    = DEFAULT_MIN_MINS;
int    minOverlap = DEFAULT_MIN_OVERLAP;
int    show       = DEFAULT_SHOW;
int    numThreads = DEFAULT_THREADS;
int    w          = DEFAULT_WINDOW_LEN;

std::string    patternFilename = ""; // name of pattern input file
std::string    rankFilename    = ""; // name of k-mer rank table binary input file
std::string    fastqFilename1  = ""; // name of FASTQ input file containing Read 1
std::string    fastqFilename2  = ""; // name of FASTQ input file containing Read 2
std::string    ifastqFilename  = ""; // name of interleaved FASTQ input file
StringVector   ubamFilename;         // names of unaligned Bam input files

KmerRankTable *rankTable;            // holds the k-mer rank table
Minimizer      maxMinimizer;         // limit used to identify common minimizers

PatternVector *patternVector;        // holds the input patterns
PatternMap    *patternMap;           // index of pattern minimizers

PairReader    *pairReader;           // used to get read pairs
bool           endOfInput = false;   // set to true when done getting read pairs

std::mutex     inputMutex;           // for getting read pairs
std::mutex     outputMutex;          // for writing hits to std::cout

uint64_t       numReadPairs = 0;     // number of read pairs found in the input

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stdout

void showUsage(const char *progname)
{
   std::cout
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname << " OPTION ... [ubam_filename ...] > matches"
      << NEWLINE;

   std::cout
      << NEWLINE
      << "These options are required:" << NEWLINE
      << "  -pattern=filename  "
             << "name of pattern input file" << NEWLINE
      << "  -rank=filename     "
             << "name of binary  input file containing the k-mer rank table"
	     << NEWLINE;

   std::cout
      << NEWLINE
      << "The following are optional:" << NEWLINE
      << "  -fastq1=filename   "
             << "name of FASTQ Read 1 input file" << NEWLINE
      << "  -fastq2=filename   "
             << "name of FASTQ Read 2 input file" << NEWLINE
      << "  -ifastq=filename   "
             << "name of interleaved FASTQ input file (may be /dev/stdin)" << NEWLINE
      << "  -maxins=N          "
             << "maximum insert size in bases, default is "
	     << DEFAULT_MAX_INSERT << NEWLINE
      << "  -maxrank=N         "
             << "maximum rank percentile of minimizers, default is "
             << doubleToString(DEFAULT_MAX_RANK) << NEWLINE
      << "  -minbases=N        "
             << "minimum percentile of matching bases, default is "
	     << doubleToString(DEFAULT_MIN_BASES) << NEWLINE
      << "  -minmins=N         "
             << "minimum number of matching minimizers, default is "
	     << DEFAULT_MIN_MINS << NEWLINE
      << "  -minov=N           "
             << "minimum overlap in number of bases, default is "
             << DEFAULT_MIN_OVERLAP << NEWLINE
      << "  -show=N            "
             << "show best only (1) or all patterns (0) matching a read pair, "
             << "default is " << DEFAULT_SHOW << NEWLINE
      << "  -threads=N         "
             << "number of threads, default is "
	     << DEFAULT_THREADS << NEWLINE
      << "  -w=N               "
             << "window length in number of bases, default is "
             << DEFAULT_WINDOW_LEN << NEWLINE;
}

//------------------------------------------------------------------------------------
// parseArgs() parses the command-line arguments and returns true if all are valid

bool parseArgs(int argc, char *argv[])
{
   for (int i = 1; i < argc; i++)
   {
      std::string arg = argv[i];
      if (arg.length() == 0)
         continue;

      if (arg[0] != '-')
      {
         ubamFilename.push_back(arg);
	 continue;
      }

      // found an option
      StringVector opt;

      if (splitString(arg, opt, '=') != 2)
         return false; // incorrect option format

      if (doubleOpt(opt, "maxrank",  maxRank)         ||
          doubleOpt(opt, "minbases", minBases)        ||
	  intOpt   (opt, "maxins",   maxInsert)       ||
	  intOpt   (opt, "minmins",  minMins)         ||
          intOpt   (opt, "minov",    minOverlap)      ||
          intOpt   (opt, "show",     show)            ||
	  intOpt   (opt, "threads",  numThreads)      ||
          intOpt   (opt, "w",        w)               ||
          stringOpt(opt, "pattern",  patternFilename) ||
	  stringOpt(opt, "rank",     rankFilename)    ||
	  stringOpt(opt, "fastq1",   fastqFilename1)  ||
	  stringOpt(opt, "fastq2",   fastqFilename2)  ||
	  stringOpt(opt, "ifastq",   ifastqFilename))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   if (ubamFilename.size() > 0) // input is coming from unaligned Bam files
   {
      if (fastqFilename1 != "" || fastqFilename2 != "" || ifastqFilename != "")
         return false; // invalid specification of file names
   }
   else if (ifastqFilename != "") // input is coming from an interleaved FASTQ file
   {
      if (fastqFilename1 != "" || fastqFilename2 != "")
         return false; // invalid specification of file names
   }
   else // input is coming from a pair of FASTQ files
      if (fastqFilename1 == "" || fastqFilename2 == "")
         return false; // invalid specification of file names

   return (maxRank > 0.0 && maxRank <= 100.0 && minBases > 0.0 && minBases <= 100.0 &&
	   maxInsert >= 100 && minMins > 0 && minOverlap > 0 &&
	   (show == 0 || show == 1) && numThreads > 0 && numThreads <= 64 &&
	   w > 0 && w < 256 && patternFilename != "" && rankFilename != "");
}

//------------------------------------------------------------------------------------
// writePattern() writes a pattern to stdout

void writePattern(const Pattern& pattern, std::string sequence, int matchingBases,
                  int possible, int insertSize)
{
   std::cout << "pattern " << pattern.name
             << TAB << sequence
	     << TAB << matchingBases
	     << TAB << possible
	     << TAB << std::fixed << std::setprecision(1)
                    << 100.0 * matchingBases / possible
             << TAB << insertSize;

   int numAnnotations = pattern.annotation.size();

   for (int i = 0; i < numAnnotations; i++)
      std::cout << TAB << pattern.annotation[i];

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeRead() writes a read to stdout

void writeRead(const std::string& name, int leadingBlanks,
               const std::string& sequence, int matchingBases, int possible)
{
   std::cout << "read " << name
             << TAB << std::string(leadingBlanks, ' ') << sequence
	     << TAB << matchingBases
	     << TAB << possible
	     << TAB << std::fixed << std::setprecision(1)
                    << 100.0 * matchingBases / possible
             << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeMatch() writes a match to stdout

void writeMatch(const std::string& name1, const std::string& sequence1,
                const std::string& name2, const std::string& sequence2,
		const Match& match)
{
   const Pattern& pattern = (*patternVector)[match.c1.index];

   int offset1    = match.c1.offset;
   int offset2    = match.c2.offset;
   int offsetDiff = match.offsetDiff();

   int plen = pattern.sequence.length();

   int len1 = std::min((int)sequence1.length(), plen - offset1);
   int len2 = std::min((int)sequence2.length(), plen - offset2);

   int displayLen = std::max(len1, len2 + offsetDiff) + 2;

   int insertSize = offsetDiff + sequence2.length();

   int leadingBlanks = offsetDiff + (offset2 >= plen - pattern.rightBases ? 2 :
                                    (offset2 >= pattern.leftBases ? 1 : 0));

   writePattern(pattern, pattern.displaySequence.substr(offset1, displayLen),
                match.matchingBases(), sequence1.length() + sequence2.length(),
		insertSize);

   writeRead(name1, 0, sequence1, match.c1.matchingBases, sequence1.length());

   writeRead(name2, leadingBlanks, sequence2, match.c2.matchingBases,
             sequence2.length());
}

//------------------------------------------------------------------------------------
// processOrientation() matches the given read pair to the set of patterns

void processOrientation(const std::string& name1, const std::string& sequence1,
                        const std::string& name2, const std::string& sequence2)
{
   MatchVector matchVector;

   getMatches(sequence1, sequence2, patternVector, patternMap, w, rankTable,
              maxMinimizer, minBases, minMins, maxInsert, (show == 1), matchVector);

   int numMatches = matchVector.size();
   if (numMatches == 0)
      return; // no matches

   std::string revcomp = stringReverseComplement(sequence2);

   // determine the number of matches with valid overlaps
   int numValid = 0;

   while (numValid < numMatches &&
          validOverlaps(sequence1, revcomp, patternVector, minBases, minOverlap,
                        matchVector[numValid]))
      numValid++;

   if (numValid == 0)
      return; // no valid matches

   // write valid matches

   if (numThreads > 1)
      outputMutex.lock();

   for (int i = 0; i < numValid; i++)
      writeMatch(name1, sequence1, name2, revcomp, matchVector[i]);

   if (numThreads > 1)
      outputMutex.unlock();
}

//------------------------------------------------------------------------------------
// getBatch() gets the next batch of read pairs and returns the number of read pairs
// obtained; if less than a full batch, end-of-input was reached; if an exception was
// raised, its message is provided

int getBatch(StringVector& name1, StringVector& seq1,
             StringVector& name2, StringVector& seq2, std::string& message)
{
   int count = 0;

   if (!endOfInput)
   {
      if (numThreads > 1)
         inputMutex.lock();

      if (!endOfInput)
      {
         try
	 {
            while (count < THREAD_BATCH_SIZE &&
                   pairReader->getNextPair(name1[count], seq1[count],
                                           name2[count], seq2[count]))
               count++;

	    numReadPairs += count;
	 }
	 catch (const std::runtime_error& error)
	 {
            message = error.what();
	 }

	 if (message != "" || count < THREAD_BATCH_SIZE)
            endOfInput = true;
      }

      if (numThreads > 1)
         inputMutex.unlock();
   }

   return count;
}

//------------------------------------------------------------------------------------
// processBatch() processes the given batch of read pairs; if an exception is raised,
// its message is provided

void processBatch(int count,
                  const StringVector& name1, const StringVector& seq1,
                  const StringVector& name2, const StringVector& seq2,
		  std::string& message)
{
   try
   {
      for (int i = 0; i < count; i++)
      {
         processOrientation(name1[i], seq1[i], name2[i], seq2[i]);
         processOrientation(name2[i], seq2[i], name1[i], seq1[i]);
      }
   }
   catch (const std::runtime_error& error)
   {
      message = error.what();

      if (numThreads > 1)
         inputMutex.lock();

      endOfInput = true;

      if (numThreads > 1)
         inputMutex.unlock();
   }
}

//------------------------------------------------------------------------------------
// threadWork() contains the work that each thread performs, processing read pairs
// until end-of-input is detected; if an exception occurs, its message is provided

void threadWork(std::string *message)
{
   StringVector name1(THREAD_BATCH_SIZE, "");
   StringVector seq1 (THREAD_BATCH_SIZE, "");
   StringVector name2(THREAD_BATCH_SIZE, "");
   StringVector seq2 (THREAD_BATCH_SIZE, "");

   int count;

   // keep going until an exception is raised or end-of-input is reached
   do
   {
      count = getBatch(name1, seq1, name2, seq2, *message);

      if (*message == "" && count > 0)
         processBatch(count, name1, seq1, name2, seq2, *message);
   }
   while (*message == "" && count == THREAD_BATCH_SIZE);
}

//------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   if (!parseArgs(argc, argv))
   {
      showUsage(argv[0]);
      return 1;
   }

   try
   {
      rankTable = readRankTable(rankFilename);

      maxMinimizer = (maxRank / 100) * numKmers(rankTable->k);

      StringVector annotationHeading;

      patternVector = readPatterns(patternFilename, annotationHeading);

      if (patternVector->size() == 0)
         throw std::runtime_error("no patterns in " + patternFilename);

      // write output heading line
      std::cout << VERSION_NAME
                << TAB << "sequence"
		<< TAB << "matching bases"
		<< TAB << "possible"
		<< TAB << "% match"
		<< TAB << "insert size";

      int numAnnotations = annotationHeading.size();

      for (int i = 0; i < numAnnotations; i++)
         std::cout << TAB << annotationHeading[i];

      std::cout << NEWLINE;

      patternMap = createPatternMap(patternVector, w, rankTable, maxMinimizer);

      UbamPairReader             *ubam   = NULL;
      InterleavedFastqPairReader *ifastq = NULL;
      FastqPairReader            *fastq  = NULL;

      if (ubamFilename.size() > 0) // unaligned Bam files
      {
         ubam = new UbamPairReader();
	 ubam->open(ubamFilename);
	 pairReader = ubam;
      }
      else if (ifastqFilename != "") // interleaved FASTQ file
      {
         ifastq = new InterleavedFastqPairReader();
	 ifastq->open(ifastqFilename);
	 pairReader = ifastq;
      }
      else // pair of FASTQ files
      {
         fastq = new FastqPairReader();
	 fastq->open(fastqFilename1, fastqFilename2);
	 pairReader = fastq;
      }

      std::string **message = new std::string *[numThreads];
      std::thread **thread  = new std::thread *[numThreads];

      // start all threads
      for (int i = 0; i < numThreads; i++)
      {
         message[i] = new std::string("");
	 thread[i]  = new std::thread(threadWork, message[i]);
      }

      // wait for each thread to finish
      for (int i = 0; i < numThreads; i++)
      {
         thread[i]->join();
	 delete thread[i];
      }

      // all of the threads have finished; check for any exceptions
      for (int i = 0; i < numThreads; i++)
         if (*message[i] != "")
            throw std::runtime_error(*message[i]);

      if (ubam)   ubam->close();
      if (ifastq) ifastq->close();
      if (fastq)  fastq->close();

      delete pairReader;

      std::cout << "read-pairs " << numReadPairs << NEWLINE;
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
