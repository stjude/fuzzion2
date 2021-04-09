//------------------------------------------------------------------------------------
//
// pairread.h - module relating to paired reads
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef PAIRREAD_H
#define PAIRREAD_H

#include <string>

//------------------------------------------------------------------------------------

class PairReader // abstract class for getting paired reads
{
public:
   PairReader() { }

   virtual ~PairReader() { }

   virtual bool getNextPair(std::string& name1, std::string& sequence1,
                            std::string& name2, std::string& sequence2) = 0;
};

bool namesMatch(const std::string& name1, const std::string& name2);

//------------------------------------------------------------------------------------
#endif
