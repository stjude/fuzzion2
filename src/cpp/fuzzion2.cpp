//------------------------------------------------------------------------------------
//
// fuzzion2.cpp - this program is a "fuzzy fusion finder"; it does fuzzy matching of
//                reads to patterns
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "hit.h"
#include "read.h"
#include "version.h"
#include <mutex>
#include <thread>

const String VERSION_NAME = FUZZION2 + CURRENT_VERSION;
const String DOC_LINK     = "https://github.com/stjude/fuzzion2";

const int    READ_BATCH_SIZE        = 10000; // number of reads in a full batch
const int    MAX_MID_LEN            = 25;    // visualization parameter

// option defaults
const double DEFAULT_MAX_RANK       = 99.9;  // max rank percentile of minimizers
const double DEFAULT_MIN_BASES      = 90.0;  // min percentile of matching bases
const int    DEFAULT_MAX_INSERT     = 600;   // max insert size in bases
const int    DEFAULT_MAX_TRIM       = 5;     // max bases, second ahead of first
const int    DEFAULT_MIN_OVERLAP    = 7;     // min length of overlap in #bases
const int    DEFAULT_SHOW           = 1;     // setting of -show option
const int    DEFAULT_SINGLE         = 0;     // setting of -single option
const int    DEFAULT_THREADS        = 8;     // number of threads
const int    DEFAULT_WINDOW_LEN     = 10;    // length of windows in #bases

double maxRank    = DEFAULT_MAX_RANK;
double minBases   = DEFAULT_MIN_BASES;
int    maxInsert  = DEFAULT_MAX_INSERT;
int    maxTrim    = DEFAULT_MAX_TRIM;
int    minOverlap = DEFAULT_MIN_OVERLAP;
int    show       = DEFAULT_SHOW;
int    single     = DEFAULT_SINGLE;
int    numThreads = DEFAULT_THREADS;
int    w          = DEFAULT_WINDOW_LEN;

String logFilename     = "";  // name of output file for logging

String patternFilename = "";  // name of pattern input file
String rankFilename    = "";  // name of k-mer rank table binary input file

String fastqFilename1  = "";  // name of FASTQ input file containing Read 1
String fastqFilename2  = "";  // name of FASTQ input file containing Read 2
String ifastqFilename  = "";  // name of interleaved FASTQ input file
String ubamFilename    = "";  // name of unaligned Bam input file
StringVector inputFilename;   // names of FASTQ/Bam input files listed on command line

KmerRankTable *rankTable;     // holds the k-mer rank table
Minimizer      maxMinimizer;  // limit used to identify common minimizers

PatternVector *patternVector; // holds the input patterns
PatternMap    *patternMap;    // index of pattern minimizers
BoolVector    *inPatternMap;  // indicates which minimizers are in patternMap

std::mutex     inputMutex;    // enforces mutual exclusion when getting reads
std::mutex     outputMutex;   // enforces mutual exclusion when writing hits

SingleReader  *singleReader = nullptr;  // for getting single reads from input files
PairReader    *pairReader   = nullptr;  // for getting read pairs from input files
bool           endOfInput   = false;    // set to true when done getting reads
uint64_t       numReads     = 0;        // number of reads found in the input

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
      << "Specify -fastq1 and -fastq2, or -ifastq or -ubam," << NEWLINE
      << "or list the names of FASTQ and Bam files on the command line" << NEWLINE
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
      << "  -minov=N      "
             << "minimum overlap in number of bases. . . . . . . . . default "
             << DEFAULT_MIN_OVERLAP << NEWLINE
      << "  -show=N       "
             << "show best only (1) or all patterns (0) that match . default "
             << DEFAULT_SHOW << NEWLINE
      << "  -single=N     "
             << "show single-read (1) or read-pair (0) matches . . . default "
             << DEFAULT_SINGLE << NEWLINE
      << "  -threads=N    "
             << "number of threads . . . . . . . . . . . . . . . . . default "
             << DEFAULT_THREADS << NEWLINE
      << "  -w=N          "
             << "window length in number of bases. . . . . . . . . . default "
             << DEFAULT_WINDOW_LEN << NEWLINE;

   std::cerr
      << NEWLINE
      << "For more information, see documentation at " << DOC_LINK << NEWLINE;
}

//------------------------------------------------------------------------------------
// parseArgs() parses the command-line arguments and returns true if all are valid

bool parseArgs(const int argc, const char *argv[])
{
   for (int i = 1; i < argc; i++)
   {
      const String arg = argv[i];
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
          intOpt   (opt, "minov",    minOverlap)      ||
          intOpt   (opt, "show",     show)            ||
          intOpt   (opt, "single",   single)          ||
          intOpt   (opt, "threads",  numThreads)      ||
          intOpt   (opt, "w",        w)               ||
          stringOpt(opt, "log",      logFilename)     || // log option for debugging
          stringOpt(opt, "pattern",  patternFilename) ||
          stringOpt(opt, "rank",     rankFilename)    ||
          stringOpt(opt, "fastq1",   fastqFilename1)  ||
          stringOpt(opt, "fastq2",   fastqFilename2)  ||
          stringOpt(opt, "ifastq",   ifastqFilename)  ||
          stringOpt(opt, "ubam",     ubamFilename))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   return (maxRank > 0.0 && maxRank <= 100.0 &&
           minBases >= 50.0 && minBases <= 100.0 &&
           maxInsert > 0 && maxTrim >= 0 && minOverlap >= TRIMER_LEN &&
           (show == 0 || show == 1) && (single == 0 || single == 1) &&
           numThreads > 0 && numThreads <= 64 && w > 0 && w < 256 &&
           patternFilename != "" && rankFilename != "");
}

//------------------------------------------------------------------------------------
// writeHits() writes to stdout each hit in the given vector; the hits are then
// de-allocated

void writeHits(HitVector& hitVector)
{
   const int numHits = hitVector.size();
   if (numHits == 0)
      return;

   if (numThreads > 1)
      outputMutex.lock();

   for (int i = 0; i < numHits; i++)
      hitVector[i]->write(std::cout);

   if (numThreads > 1)
      outputMutex.unlock();

   // cleanup
   for (int i = 0; i < numHits; i++)
      delete hitVector[i];
}

//------------------------------------------------------------------------------------
// processSingleRead() finds matches of the given single read and writes those that
// are qualified as hits

void processSingleRead(const String& readName, const String& readStr)
{
   Seq *readSeq;
   SingleMatchMap mmap;

   if (getSingleMatches(readStr, w, rankTable, maxMinimizer, *patternMap,
                        *inPatternMap, *patternVector, minBases, minOverlap, readSeq,
                        mmap) == 0)
      return; // no matches

   SingleMatchMap bestMap;
   selectBestSingleMatches((show == 1), mmap, bestMap);

   HitVector hitVector;

   for (SingleMatchMap::const_iterator mpos = bestMap.cbegin();
        mpos != bestMap.cend(); ++mpos)
   {
      const Pattern& pattern = *(*patternVector)[mpos->first];
      const SingleMatchVector& mvector = mpos->second;

      Hit *hit = createHitFromSingleMatch(pattern, MAX_MID_LEN, minBases, minOverlap,
                                          readName, *readSeq, *mvector[0]);
      hitVector.push_back(hit);
   }

   writeHits(hitVector);

   // cleanup
   delete readSeq;
   freeSingleMatchMap(mmap);
   freeSingleMatchMap(bestMap);
}

//------------------------------------------------------------------------------------
// processReadPair() finds matches of the given read pair and writes those that are
// qualified as hits

void processReadPair(const String& readName1, const String& readStr1,
                     const String& readName2, const String& readStr2)
{
   Seq *readSeq1, *readSeq2;
   PairMatchMap pairMap;

   if (getPairMatches(readStr1, readStr2, w, rankTable, maxMinimizer, *patternMap,
                      *inPatternMap, *patternVector, minBases, minOverlap, maxInsert,
                      maxTrim, readSeq1, readSeq2, pairMap) == 0)
      return; // no matches

   PairMatchMap bestMap;

   if (show == 1)
      selectBestPairMatch(pairMap, bestMap);

   const PairMatchMap& pmap = (show == 1 ? bestMap : pairMap);

   HitVector hitVector;

   for (PairMatchMap::const_iterator ppos = pmap.cbegin();
        ppos != pmap.cend(); ++ppos)
   {
      const Pattern& pattern = *(*patternVector)[ppos->first];
      const PairMatch *pairMatch = ppos->second;

      Hit *hit = createHitFromPairMatch(pattern, MAX_MID_LEN, minBases, minOverlap,
                                        readName1, *readSeq1, readName2, *readSeq2,
                                        *pairMatch);
      hitVector.push_back(hit);
   }

   writeHits(hitVector);

   // cleanup
   delete readSeq1;
   delete readSeq2;
   freePairMatchMap(pairMap);
   freePairMatchMap(bestMap);
}

//------------------------------------------------------------------------------------
// processBatch() processes the given batch of reads; if an exception is raised, its
// message is provided

void processBatch(const int count, String& message,
                  String *readName1, String *readStr1,
                  String *readName2=nullptr, String *readStr2=nullptr)
{
   try
   {  // process both orientations

      if (single == 1)
         for (int i = 0; i < count; i++)
         {
            processSingleRead(readName1[i], readStr1[i]);

            const String revcomp = stringReverseComplement(readStr1[i]);
            processSingleRead(readName1[i], revcomp);
         }
      else
         for (int i = 0; i < count; i++)
         {
            processReadPair(readName1[i], readStr1[i], readName2[i], readStr2[i]);
            processReadPair(readName2[i], readStr2[i], readName1[i], readStr1[i]);
         }
   }
   catch (const Error& error)
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
// getBatch() gets the next batch of reads and returns the number of reads or read
// pairs obtained; if less than a full batch, end-of-input was reached; if an
// exception was raised, its message is provided

int getBatch(const int batchSize, String& message,
             String *readName1, String *readStr1,
             String *readName2=nullptr, String *readStr2=nullptr)
{
   if (endOfInput)
      return 0;

   int count = 0;

   if (numThreads > 1)
      inputMutex.lock();

   if (!endOfInput)
   {
      try
      {
         if (single == 1)
            while (count < batchSize &&
                   singleReader->getNext(readName1[count], readStr1[count]))
               count++;
         else
            while (count < batchSize &&
                   pairReader->getNext(readName1[count], readStr1[count],
                                       readName2[count], readStr2[count]))
               count++;
      }
      catch (const Error& error) { message = error.what(); }

      if (single == 1)
         numReads += count;
      else
         numReads += 2 * count;

      if (message != "" || count < batchSize)
         endOfInput = true;
   }

   if (numThreads > 1)
      inputMutex.unlock();

   return count;
}

//------------------------------------------------------------------------------------
// threadWork() contains the work that each thread performs, processing reads until
// end-of-input is detected; if an exception occurs, its message is provided

void threadWork(String *message)
{
   int count;

   if (single == 1)
   {
      const int batchSize = READ_BATCH_SIZE;

      String *readName = new String[batchSize];
      String *readStr  = new String[batchSize];

      do
      {
         count = getBatch(batchSize, *message, readName, readStr);

         if (*message == "" && count > 0)
            processBatch(count, *message, readName, readStr);
      }
      while (*message == "" && count == batchSize);

      delete[] readName;
      delete[] readStr;
   }
   else
   {
      const int batchSize = READ_BATCH_SIZE / 2;

      String *readName1 = new String[batchSize];
      String *readStr1  = new String[batchSize];
      String *readName2 = new String[batchSize];
      String *readStr2  = new String[batchSize];

      do
      {
         count = getBatch(batchSize, *message, readName1, readStr1, readName2,
                          readStr2);

         if (*message == "" && count > 0)
            processBatch(count, *message, readName1, readStr1, readName2, readStr2);
      }
      while (*message == "" && count == batchSize);

      delete[] readName1;
      delete[] readStr1;
      delete[] readName2;
      delete[] readStr2;
   }
}

//------------------------------------------------------------------------------------

int main(const int argc, const char *argv[])
{
   if (!parseArgs(argc, argv))
   {
      showUsage(argv[0]);
      return 1;
   }

   try
   {
      if (logFilename != "")
         logOpen(logFilename);

      // open all FASTQ and Bam input files to verify they exist and are as expected

      if (single == 1)
         singleReader =
            createSingleReader(fastqFilename1, fastqFilename2, ifastqFilename,
                               ubamFilename, inputFilename);
      else
         pairReader =
            createPairReader(fastqFilename1, fastqFilename2, ifastqFilename,
                             ubamFilename, inputFilename);

      // now read the pattern set and k-mer rank table, and create the pattern index
      StringVector annotationHeading;
      patternVector = readPatterns(patternFilename, annotationHeading);
      rankTable     = readRankTable(rankFilename);
      maxMinimizer  = (maxRank / 100) * numKmers(rankTable->k) - 1;
      patternMap    = createPatternMap(*patternVector, w, rankTable, maxMinimizer,
                                       inPatternMap);

      writeHitHeadingLine(std::cout, CURRENT_VERSION, annotationHeading);

      String **message     = new String *[numThreads];
      std::thread **thread = new std::thread *[numThreads];

      // start all threads
      for (int i = 0; i < numThreads; i++)
      {
         message[i] = new String("");
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
            throw Error(*message[i]);

      writeReadCountLine(std::cout, numReads);

      if (single == 1)
         singleReader->close();
      else
         pairReader->close();

      if (logFilename != "")
         logClose();
   }
   catch (const Error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
