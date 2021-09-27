//------------------------------------------------------------------------------------
//
// kmerank.cpp - this program reads a reference genome file in 2-bit format, counts
//               the k-mers in the reference genome, and produces a k-mer rank table
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "rank.h"
#include "util.h"
#include <iostream>
#include <stdexcept>

const std::string VERSION_ID  = "kmerank v1.0.2, copyright 2021 "
                                "St. Jude Children's Research Hospital";

const int DEFAULT_KMER_LENGTH = 15;

int k = DEFAULT_KMER_LENGTH;     // k-mer length

std::string refGenFilename = ""; // name of reference genome input file
std::string binaryFilename = ""; // name of binary output file
std::string textFilename   = ""; // name of text output file

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_ID << NEWLINE << NEWLINE
      << "Usage: " << progname << " OPTION ..." << NEWLINE;

   std::cerr
      << NEWLINE
      << "These options are required:" << NEWLINE
      << "  -ref=filename   "
             << "name of reference genome input file in 2-bit format" << NEWLINE
      << "  -bin=filename   "
             << "name of binary output file" << NEWLINE;

   std::cerr
      << NEWLINE
      << "The following are optional:" << NEWLINE
      << "  -k=N            "
             << "k-mer length, default is " << DEFAULT_KMER_LENGTH
	     << ", maximum is " << static_cast<int>(MAX_KMER_LENGTH) << NEWLINE
      << "  -txt=filename   "
             << "name of text output file, default is none" << NEWLINE;
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
         return false; // not an option

      StringVector opt;

      if (splitString(arg, opt, '=') != 2)
         return false; // incorrect option format

      if (intOpt(opt, "k", k) ||
          stringOpt(opt, "ref", refGenFilename) ||
	  stringOpt(opt, "bin", binaryFilename) ||
	  stringOpt(opt, "txt", textFilename))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   if (k < 1 || k > MAX_KMER_LENGTH || refGenFilename == "" || binaryFilename == "")
      return false; // missing or invalid option

   return true;
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
      KmerRankTable *table = createRankTable(k, refGenFilename);

      table->writeBinary(binaryFilename);

      if (textFilename != "")
         table->writeText(textFilename);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
