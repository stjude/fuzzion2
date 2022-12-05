//------------------------------------------------------------------------------------
//
// fuzzion2.cpp - this program is a "fuzzy fusion finder"; it does fuzzy matching of
//                read pairs to fusion patterns
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2022 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "fastq.h"
#include "hit.h"
#include "match.h"
#include "ubam.h"
#include "version.h"
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <thread>

const std::string VERSION_NAME   = FUZZION2 + CURRENT_VERSION;

const int    THREAD_BATCH_SIZE   = 100000; // number of read pairs in a full batch

const double DEFAULT_MAX_RANK    = 99.9; // default max rank percentile of minimizers
const double DEFAULT_MIN_BASES   = 90.0; // default min percentile of matching bases
const int    DEFAULT_MAX_INSERT  = 500;  // default max insert size in bases
const int    DEFAULT_MAX_TRIM    = 5;    // default max bases, second ahead of first
const int    DEFAULT_MIN_MINS    = 1;    // default min number of matching minimizers
const int    DEFAULT_MIN_OVERLAP = 5;    // default min length of overlap in #bases
const int    DEFAULT_SHOW        = 1;    // default setting of -show option
const int    DEFAULT_SINGLE      = 0;    // default setting of -single option
const int    DEFAULT_THREADS     = 8;    // default number of threads
const int    DEFAULT_WINDOW_LEN  = 10;   // default length of windows in #bases

double maxRank    = DEFAULT_MAX_RANK;
double minBases   = DEFAULT_MIN_BASES;
int    maxInsert  = DEFAULT_MAX_INSERT;
int    maxTrim    = DEFAULT_MAX_TRIM;
int    minMins    = DEFAULT_MIN_MINS;
int    minOverlap = DEFAULT_MIN_OVERLAP;
int    show       = DEFAULT_SHOW;
int    single     = DEFAULT_SINGLE;
int    numThreads = DEFAULT_THREADS;
int    w          = DEFAULT_WINDOW_LEN;

std::string    patternFilename = ""; // name of pattern input file
std::string    rankFilename    = ""; // name of k-mer rank table binary input file

std::string    fastqFilename1  = ""; // name of FASTQ input file containing Read 1
std::string    fastqFilename2  = ""; // name of FASTQ input file containing Read 2
std::string    ifastqFilename  = ""; // name of interleaved FASTQ input file
std::string    ubamFilename    = ""; // name of unaligned Bam input file
StringVector   inputFilename;        // names of input files listed on command line

KmerRankTable *rankTable;            // holds the k-mer rank table
Minimizer      maxMinimizer;         // limit used to identify common minimizers

PatternVector *patternVector;        // holds the input patterns
PatternMap    *patternMap;           // index of pattern minimizers

PairReader    *pairReader;           // used to get read pairs from input files
bool           endOfInput = false;   // set to true when done getting read pairs

std::mutex     inputMutex;           // for getting read pairs
std::mutex     outputMutex;          // for writing hits to std::cout

uint64_t       numReadPairs = 0;     // number of read pairs found in the input

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname << " OPTION ... [filename ...] > hits"
      << NEWLINE;

   std::cerr
      << NEWLINE
      << "These options are required:" << NEWLINE
      << "  -pattern=filename   "
             << "name of pattern input file" << NEWLINE
      << "  -rank=filename      "
             << "name of binary  input file containing the k-mer rank table"
	     << NEWLINE;

   std::cerr
      << NEWLINE
      << "Specify -fastq1 and -fastq2, or -ifastq or -ubam, "
      << "or list filenames on command line" << NEWLINE
      << "  -fastq1=filename    "
             << "name of FASTQ Read 1 input file" << NEWLINE
      << "  -fastq2=filename    "
             << "name of FASTQ Read 2 input file" << NEWLINE
      << "  -ifastq=filename    "
             << "name of interleaved FASTQ input file (may be /dev/stdin)" << NEWLINE
      << "  -ubam=filename      "
             << "name of unaligned Bam input file" << NEWLINE;

   std::cerr
      << NEWLINE
      << "The following are optional:" << NEWLINE
      << "   N is a numeric value, e.g., -threads=4" << NEWLINE
      << "  -maxins=N     "
             << "maximum insert size in bases. . . . . . . . . . . . default "
	     << DEFAULT_MAX_INSERT << NEWLINE
      << "  -maxrank=N    "
             << "maximum rank percentile of minimizers . . . . . . . default "
             << doubleToString(DEFAULT_MAX_RANK) << NEWLINE
      << "  -maxtrim=N    "
             << "maximum bases second read aligned ahead of first. . default "
	     << DEFAULT_MAX_TRIM << NEWLINE
      << "  -minbases=N   "
             << "minimum percentile of matching bases. . . . . . . . default "
	     << doubleToString(DEFAULT_MIN_BASES) << NEWLINE
      << "  -minmins=N    "
             << "minimum number of matching minimizers . . . . . . . default "
	     << DEFAULT_MIN_MINS << NEWLINE
      << "  -minov=N      "
             << "minimum overlap in number of bases. . . . . . . . . default "
             << DEFAULT_MIN_OVERLAP << NEWLINE
      << "  -show=N       "
             << "show best only (1) or all patterns (0) that match . default "
	     << DEFAULT_SHOW << NEWLINE
      << "  -single=N     "
             << "show single-read (1) or just read-pair (0) matches. default "
	     << DEFAULT_SINGLE << NEWLINE
      << "  -threads=N    "
             << "number of threads . . . . . . . . . . . . . . . . . default "
	     << DEFAULT_THREADS << NEWLINE
      << "  -w=N          "
             << "window length in number of bases. . . . . . . . . . default "
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
         inputFilename.push_back(arg);
	 continue;
      }

      // found an option
      StringVector opt;

      if (splitString(arg, opt, '=') != 2)
         return false; // incorrect option format

      if (doubleOpt(opt, "maxrank",  maxRank)         ||
          doubleOpt(opt, "minbases", minBases)        ||
	  intOpt   (opt, "maxins",   maxInsert)       ||
	  intOpt   (opt, "maxtrim",  maxTrim)         ||
	  intOpt   (opt, "minmins",  minMins)         ||
          intOpt   (opt, "minov",    minOverlap)      ||
          intOpt   (opt, "show",     show)            ||
	  intOpt   (opt, "single",   single)          ||
	  intOpt   (opt, "threads",  numThreads)      ||
          intOpt   (opt, "w",        w)               ||
          stringOpt(opt, "pattern",  patternFilename) ||
	  stringOpt(opt, "rank",     rankFilename)    ||
	  stringOpt(opt, "fastq1",   fastqFilename1)  ||
	  stringOpt(opt, "fastq2",   fastqFilename2)  ||
	  stringOpt(opt, "ifastq",   ifastqFilename)  ||
	  stringOpt(opt, "ubam",     ubamFilename))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   if (inputFilename.size() > 0)  // list of file names on the command line
   {
      if (fastqFilename1 != "" || fastqFilename2 != "" || ifastqFilename != "" ||
          ubamFilename != "")
         return false;
   }
   else if (fastqFilename1 != "") // pair of FASTQ input files
   {
      if (fastqFilename2 == "" || ifastqFilename != "" || ubamFilename != "")
         return false;
   }
   else if (ifastqFilename != "") // interleaved FASTQ input file
   {
      if (fastqFilename2 != "" || ubamFilename != "")
         return false;
   }
   else if (ubamFilename != "")   // unaligned Bam input file
   {
      if (fastqFilename2 != "")
         return false;
   }
   else // missing input file names
      return false;

   return (maxRank > 0.0 && maxRank <= 100.0 && minBases > 0.0 && minBases <= 100.0 &&
	   maxInsert > 0 && maxTrim >= 0 && minMins > 0 && minOverlap > 0 &&
	   (show == 0 || show == 1) && (single == 0 || single == 1) &&
	   numThreads > 0 && numThreads <= 64 && w > 0 && w < 256 &&
	   patternFilename != "" && rankFilename != "");
}

//------------------------------------------------------------------------------------
// createInputReader() analyzes the given vector of input file names and returns a
// newly allocated InputReader object for reading the input files

InputReader *createInputReader(const StringVector& inputFilename)
{
   int numInputFiles = inputFilename.size();
   if (numInputFiles == 0)
      throw std::runtime_error("no input files");

   BoolVector   processed(numInputFiles, false);
   StringVector name1(numInputFiles, "");
   StringVector name2(numInputFiles, "");

   PairReaderVector readerVector;

   for (int i = 0; i < numInputFiles; i++)
      if (isUbamFile(inputFilename[i]))
      {
         readerVector.push_back(new UbamPairReader(inputFilename[i]));
	 processed[i] = true;
      }
      else
      {
         bool interleaved;

	 if (!isFastqFile(inputFilename[i], name1[i], name2[i], interleaved))
            throw std::runtime_error("unsupported file type " + inputFilename[i]);

	 if (interleaved)
	 {
            readerVector.push_back(new InterleavedFastqPairReader(inputFilename[i]));
	    processed[i] = true;
	 }
      }

   // now pair up the unprocessed FASTQ files

   for (int i = 0; i < numInputFiles; i++)
      if (!processed[i])
      {
         for (int j = i + 1; j < numInputFiles; j++)
            if (!processed[j] && namesMatch(name1[i], name1[j]) &&
                                 namesMatch(name2[i], name2[j]))
	    {
               readerVector.push_back(new FastqPairReader(inputFilename[i],
                                                          inputFilename[j]));
	       processed[i] = true;
	       processed[j] = true;
	       break;
	    }

	 if (!processed[i])
            throw std::runtime_error("unsupported FASTQ file " + inputFilename[i]);
      }

   return new InputReader(readerVector);
}

//------------------------------------------------------------------------------------
// writeMatch() writes a match to stdout showing the alignment of reads to a pattern;
// pattern delimiters (i.e., brackets and braces) are accounted for

void writeMatch(const std::string& name1, const std::string& sequence1,
                const std::string& name2, const std::string& sequence2,
		const Match& match)
{
   const Pattern& pattern = (*patternVector)[match.c1.index];

   int offset1 = match.c1.offset;
   int offset2 = match.c2.offset;

   int leftOffset = std::min(offset1, offset2);
   int fullLength = pattern.displaySequence.length();
   int displayLen = std::min(match.insertSize() + 2, fullLength - leftOffset);

   std::string displaySeq = pattern.displaySequence.substr(leftOffset, displayLen);

   int leading1 = 0, leading2 = 0; // number of leading blanks

   if (offset1 < offset2)      // first read aligns ahead of second read
      leading2 = offset2 - offset1 +
                 (offset2 >= pattern.sequence.length() - pattern.rightBases ? 2 :
		 (offset2 >= pattern.leftBases ? 1 : 0));
   else if (offset2 < offset1) // second read aligns ahead of first read
      leading1 = offset1 - offset2 +
                 (offset1 >= pattern.sequence.length() - pattern.rightBases ? 2 :
		 (offset1 >= pattern.leftBases ? 1 : 0));

   HitPattern *hitPattern = new HitPattern(pattern.name, displaySeq,
                                           pattern.annotation, match.matchingBases(),
					   match.possible(), match.insertSize());

   HitRead *hitRead1 = new HitRead(name1, leading1, sequence1,
                                   match.c1.matchingBases);

   HitRead *hitRead2 = new HitRead(name2, leading2, sequence2,
                                   match.c2.matchingBases);

   Hit hit(hitPattern, hitRead1, hitRead2);
   hit.write();
}

//------------------------------------------------------------------------------------
// processOrientation() matches the given read pair to the set of patterns

void processOrientation(const std::string& name1, const std::string& sequence1,
                        const std::string& name2, const std::string& sequence2)
{
   MatchVector matchVector;

   getMatches(sequence1, sequence2, patternVector, patternMap, w, rankTable,
              maxMinimizer, minBases, minMins, maxInsert, maxTrim, (show == 1),
	      (single == 1), matchVector);

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

      writeHitHeadingLine(CURRENT_VERSION, annotationHeading);

      patternMap = createPatternMap(patternVector, w, rankTable, maxMinimizer);

      if (fastqFilename1 != "")
         pairReader = new FastqPairReader(fastqFilename1, fastqFilename2);
      else if (ifastqFilename != "")
         pairReader = new InterleavedFastqPairReader(ifastqFilename);
      else if (ubamFilename != "")
         pairReader = new UbamPairReader(ubamFilename);
      else
         pairReader = createInputReader(inputFilename);

      pairReader->open();

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

      pairReader->close();

      writeReadPairLine(numReadPairs);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
