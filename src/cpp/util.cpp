//------------------------------------------------------------------------------------
//
// util.cpp - module defining some utilities
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "util.h"
#include <iomanip>
#include <sstream>

//------------------------------------------------------------------------------------
// splitString() splits the given string delimited by the specified delimiter and
// saves the component strings as elements of a vector; the number of component
// strings is returned and is equal to the number of delimiters in the string plus one

int splitString(const std::string& s, StringVector& v, char delimiter)
{
   v.clear();

   const char *str = s.c_str();

   int  startpos = 0;
   char c;

   do
   {
      int i = startpos;
      c = str[i];

      while (c && c != delimiter)
         c = str[++i];

      if (i == startpos)
         v.push_back("");
      else
         v.push_back(s.substr(startpos, i - startpos));

      startpos = i + 1;
   }
   while (c);

   return v.size();
}

//------------------------------------------------------------------------------------
// stringToNonnegInt() converts a string to a nonnegative integer; -1 is returned if
// the conversion cannot be performed

int stringToNonnegInt(const std::string& s)
{
   const char *str = s.c_str();

   int  value = 0;
   int  i     = 0;
   char c     = str[0];

   while (c >= '0' && c <= '9')
   {
      value = 10 * value + c - '0';
      c = str[++i];
   }

   return (i > 0 && !c ? value : -1);
}

//------------------------------------------------------------------------------------
// stringToNonnegDouble() converts a string to a nonnegative double; -1 is returned if
// the conversion cannot be performed

double stringToNonnegDouble(const std::string& s)
{
   double value = -1.0;

   std::istringstream stream(s);
   stream >> value;

   return (value >= 0.0 ? value : -1.0);
}

//------------------------------------------------------------------------------------
// stringOpt() returns true if opt contains the named option and the option value is
// provided as a string

bool stringOpt(const StringVector& opt, std::string optname, std::string& optvalue)
{
   if (opt[0] == "-" + optname)
   {
      optvalue = opt[1];
      return true;
   }
   else
      return false;
}

//------------------------------------------------------------------------------------
// intOpt() returns true if opt contains the named option and the option value is
// provided as an integer

bool intOpt(const StringVector& opt, std::string optname, int& optvalue)
{
   if (opt[0] == "-" + optname)
   {
      optvalue = stringToNonnegInt(opt[1]);
      return true;
   }
   else
      return false;
}

//------------------------------------------------------------------------------------
// doubleOpt() returns true if opt contains the named option and the option value is
// provided as a double

bool doubleOpt(const StringVector& opt, std::string optname, double& optvalue)
{
   if (opt[0] == "-" + optname)
   {
      optvalue = stringToNonnegDouble(opt[1]);
      return true;
   }
   else
      return false;
}

//------------------------------------------------------------------------------------
// intToString() returns a string representation of the given integer

std::string intToString(int i)
{
   char buffer[100];

   std::sprintf(buffer, "%d", i);

   return buffer;
}

//------------------------------------------------------------------------------------
// doubleToString() returns a string representation of the given double

std::string doubleToString(double d)
{
   char buffer[100];

   std::sprintf(buffer, "%.1f", d);

   return buffer;
}

//------------------------------------------------------------------------------------
// intToStringLeadingZeros() returns a string representation of the given integer of
// the specified width, with leading zeros

std::string intToStringLeadingZeros(int i, int width)
{
   std::ostringstream stream;
   stream << std::setfill('0') << std::setw(width) << i;

   return stream.str();
}
