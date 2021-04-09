//------------------------------------------------------------------------------------
//
// pairread.cpp - module relating to paired reads
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "pairread.h"
#include <cstring>

//------------------------------------------------------------------------------------
// namesMatch() returns true if the given read names match

bool namesMatch(const std::string& name1, const std::string& name2)
{
   int len = name1.length();
   if (len != name2.length())
      return false; // names have unequal lengths

   // see if one name ends with "1" and the other name ends with "2"
   int last = len - 1;

   if (last > 0 && (name1[last] == '1' && name2[last] == '2' ||
                    name1[last] == '2' && name2[last] == '1'))
      return (std::strncmp(name1.c_str(), name2.c_str(), last) == 0);

   return (name1 == name2);
}
