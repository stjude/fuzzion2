//------------------------------------------------------------------------------------
//
// fuzzion2.cpp - this program is a "fuzzy fusion finder"; it does fuzzy matching of
//                read-pairs to fusion patterns
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "fastq.h"
#include "pattern.h"
#include "ubam.h"
#include <cmath>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <thread>

const std::string VERSION_ID     = "fuzzion2 v1.0.0, copyright 2021 "
                                   "St. Jude Children's Research Hospital";

const int    THREAD_BATCH_SIZE   = 100000; // number of read-pairs in a full batch

const double DEFAULT_MAX_RANK    = 95;  // default max rank percentile of minimizers
const double DEFAULT_MIN_BASES   = 90;  // default min percentile of matching bases
const int    DEFAULT_MAX_INSERT  = 500; // default max insert size in bases
const int    DEFAULT_MIN_MINS    = 3;   // default min number of matching minimizers
const int    DEFAULT_MIN_OVERLAP = 5;   // default min length of overlap in #bases
const int    DEFAULT_THREADS     = 8;   // default number of threads
const int    DEFAULT_WINDOW_LEN  = 5;   // default length of windows in #bases

double maxRank    = DEFAULT_MAX_RANK;
double minBases   = DEFAULT_MIN_BASES;
int    maxInsert  = DEFAULT_MAX_INSERT;
int    minMins    = DEFAULT_MIN_MINS;
int    minOverlap = DEFAULT_MIN_OVERLAP;
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

PairReader    *pairReader;           // used to get read-pairs
bool           endOfInput = false;   // set to true when done getting read-pairs

std::mutex     inputMutex;           // for getting read-pairs
std::mutex     outputMutex;          // for writing hits to std::cout

uint64_t       numReadPairs = 0;     // number of read-pairs found in the input

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stdout

void showUsage(const char *progname)
{
   std::cout
      << VERSION_ID << NEWLINE << NEWLINE
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
	   numThreads > 0 && numThreads <= 64 && w > 0 && w < 256 &&
	   patternFilename != "" && rankFilename != "");
}

//------------------------------------------------------------------------------------
// computeMinMatches() returns the minimum number of matching bases for the given
// sequence length

int computeMinMatches(int seqlen)
{
   return std::ceil((minBases / 100) * seqlen);
}

//------------------------------------------------------------------------------------
// hasCandidates() returns true if each of the given sequences has at least one
// candidate match in the pattern set; if true is returned, the candidate matches are
// provided in the vectors and it is the caller's obligation to de-allocate these
// vectors; if false is returned, the caller must not access or de-allocate the
// vectors

bool hasCandidates(const std::string& sequence1,  const std::string& sequence2,
                   CandidateVector *& candidate1, CandidateVector *& candidate2)
{
   candidate1 = findCandidates(sequence1, patternMap, patternVector,
		               w, rankTable, maxMinimizer, minMins,
			       computeMinMatches(sequence1.length()));

   int numCandidates1 = candidate1->size();
   if (numCandidates1 == 0)
   {
      delete candidate1;
      return false;
   }

   std::string revcomp = stringReverseComplement(sequence2);

   PatternVector eligiblePattern(patternVector->size(), NULL);

   for (int i = 0; i < numCandidates1; i++)
   {
      int index = (*candidate1)[i].location.index;
      eligiblePattern[index] = (*patternVector)[index];
   }

   candidate2 = findCandidates(revcomp, patternMap, &eligiblePattern,
		               w, rankTable, maxMinimizer, minMins,
			       computeMinMatches(sequence2.length()));

   if (candidate2->size() > 0)
      return true;

   delete candidate1;
   delete candidate2;
   return false;
}

//------------------------------------------------------------------------------------
// bestMatch() returns the best matching pattern for the given sequences or NULL if
// none; the match offsets within the pattern sequence are provided

const Pattern *bestMatch(const std::string& sequence1, const std::string& sequence2,
                         int& bestOffset1, int& bestOffset2)
{
   CandidateVector *candidate1, *candidate2;

   if (!hasCandidates(sequence1, sequence2, candidate1, candidate2))
      return NULL; // no match

   int numCandidates1 = candidate1->size();
   int numCandidates2 = candidate2->size();

   int bestIndex = -1, maxMatches;

   for (int i = 0; i < numCandidates1; i++)
   {
      const Location& loc1 = (*candidate1)[i].location;

      for (int j = 0; j < numCandidates2; j++)
      {
         const Location& loc2 = (*candidate2)[j].location;

	 if (loc1.index != loc2.index)
            continue; // different pattern

	 if (loc1.offset > loc2.offset)
            continue; // incompatible offsets

	 if (loc2.offset - loc1.offset + sequence2.length() > maxInsert)
            continue; // alignment is longer than the maximum insert size

	 int matches = (*candidate1)[i].matches + (*candidate2)[j].matches;

	 if (bestIndex >= 0 && (matches < maxMatches || matches == maxMatches &&
             loc2.offset - loc1.offset >= bestOffset2 - bestOffset1))
            continue; // not the best match

	 // save the best match
	 bestIndex   = loc1.index;
	 bestOffset1 = loc1.offset;
	 bestOffset2 = loc2.offset;
	 maxMatches  = matches;
      }
   }

   delete candidate1;
   delete candidate2;

   return (bestIndex >= 0 ? (*patternVector)[bestIndex] : NULL);
}

//------------------------------------------------------------------------------------
// overlapLeft() returns true if the given sequence aligned at the specified offset
// overlaps the left side of the pattern

bool overlapLeft(const std::string& sequence, const Pattern *pattern, int offset)
{
   int maxPatternLen = pattern->sequence.length() - offset;
   int patternLen    = std::min((int)sequence.length(), maxPatternLen);
   int overlapBases  = std::min(patternLen, pattern->leftBases - offset);

   if (overlapBases < minOverlap)
      return false; // insufficient number of overlapping bases

   if (overlapBases == patternLen)
      return true;  // full overlap; we have already verified agreement of bases

   // verify agreement of bases for partial overlap
   int lcslen = lengthOfLCS(sequence, 0, overlapBases, pattern->sequence, offset,
                            overlapBases);

   return (lcslen >= computeMinMatches(overlapBases));
}

//------------------------------------------------------------------------------------
// overlapRight() returns true if the given sequence aligned at the specified offset
// overlaps the right side of the pattern

bool overlapRight(const std::string& sequence, const Pattern *pattern, int offset)
{
   int maxPatternLen = pattern->sequence.length() - offset;
   int patternLen    = std::min((int)sequence.length(), maxPatternLen);
   int overlapBases  = patternLen + std::min(0, pattern->rightBases - maxPatternLen);

   if (overlapBases < minOverlap)
      return false; // insufficient number of overlapping bases

   if (overlapBases == patternLen)
      return true;  // full overlap; we have already verified agreement of bases

   // verify agreement of bases for partial overlap
   int lcslen = lengthOfLCS(sequence, sequence.length() - overlapBases, overlapBases,
                            pattern->sequence, offset + patternLen - overlapBases,
			    overlapBases);

   return (lcslen >= computeMinMatches(overlapBases));
}

//------------------------------------------------------------------------------------
// processOrientation() matches the given paired reads to the set of patterns

void processOrientation(const std::string& name1, const std::string& sequence1,
                        const std::string& name2, const std::string& sequence2)
{
   int offset1, offset2;

   const Pattern *pattern = bestMatch(sequence1, sequence2, offset1, offset2);
   if (!pattern)
      return;

   std::string revcomp = stringReverseComplement(sequence2);

   bool foundHit = false;

   if (pattern->hasBraces)
   {
      if (overlapLeft (sequence1, pattern, offset1) &&
          overlapRight(sequence1, pattern, offset1) ||
	  overlapLeft (revcomp,   pattern, offset2) &&
	  overlapRight(revcomp,   pattern, offset2))
         foundHit = true;
   }
   else // brackets
      if (overlapLeft (sequence1, pattern, offset1) &&
          overlapRight(revcomp,   pattern, offset2))
         foundHit = true;

   if (!foundHit)
      return;

   int insertSize = offset2 - offset1 + sequence2.length();

   int len     = pattern->sequence.length();
   int len1    = std::min((int)sequence1.length(), len - offset1);
   int len2    = std::min((int)sequence2.length(), len - offset2);
   int dispLen = std::max(len1, len2 + offset2 - offset1) + 2;

   int leadingBlanks = offset2 - offset1;

   if (offset2 >= pattern->sequence.length() - pattern->rightBases)
      leadingBlanks += 2; // account for two delimiter characters
   else if (offset2 >= pattern->leftBases)
      leadingBlanks++;    // account for one delimiter character

   // write the match
   if (numThreads > 1)
      outputMutex.lock();

   std::cout << "pattern " << pattern->name << " " << insertSize
             << TAB << pattern->displaySequence.substr(offset1, dispLen)
	     << NEWLINE;

   std::cout << name1
             << TAB << sequence1
	     << NEWLINE;

   std::cout << name2
             << TAB << std::string(leadingBlanks, ' ') << revcomp
	     << NEWLINE;

   if (numThreads > 1)
      outputMutex.unlock();
}

//------------------------------------------------------------------------------------
// getBatch() gets the next batch of read-pairs and returns the number of read-pairs
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
// processBatch() processes the given batch of read-pairs; if an exception is raised,
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
// threadWork() contains the work that each thread performs, processing read-pairs
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

      patternVector = readPatterns(patternFilename);

      if (patternVector->size() == 0)
         throw std::runtime_error("no patterns in " + patternFilename);

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

      std::cout << numReadPairs
                << TAB << "read pairs"
		<< NEWLINE;
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
