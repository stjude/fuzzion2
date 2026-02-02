//------------------------------------------------------------------------------------
//
// util.cpp - module defining some utilities
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "util.h"
#include <fstream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>

std::ofstream *logfile = nullptr; // is non-null when the log file is open
std::mutex     logMutex; // enforces mutual exclusion when writing to log file

//------------------------------------------------------------------------------------
// getline() reads the next line from the input stream and returns true, or returns
// false if end-of-file was encountered

bool getline(std::istream& stream, String& line)
{
   if (!std::getline(stream, line))
      return false;

   // check for Windows line ending
   const int len = line.length();
   if (len > 0 && line[len - 1] == CRETURN)
      line = line.substr(0, len - 1); // remove trailing carriage return

   return true;
}

//------------------------------------------------------------------------------------
// splitString() splits the given string delimited by the specified delimiter and
// saves the component strings as elements of a vector; the number of component
// strings is returned and is equal to the number of delimiters in the string plus one

int splitString(const String& s, StringVector& v, const char delimiter)
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
// stringToInt() converts a string to a nonnegative integer; -1 is returned if the
// conversion cannot be performed

int stringToInt(const String& s)
{
   const char *str = s.c_str();

   const uint64_t MAXINT = std::numeric_limits<int>::max();
   uint64_t value  = 0;

   int  i = 0;
   char c = str[0];

   while (c >= '0' && c <= '9')
   {
      value = 10 * value + c - '0';

      if (value > MAXINT) // integer too large
         return -1;

      c = str[++i];
   }

   return (i > 0 && !c ? value : -1);
}

//------------------------------------------------------------------------------------
// stringToDouble() converts a string to a nonnegative double; -1.0 is returned if the
// conversion cannot be performed

double stringToDouble(const String& s)
{
   double value = -1.0;

   std::istringstream stream(s);
   stream >> value;

   return (value >= 0.0 ? value : -1.0);
}

//------------------------------------------------------------------------------------
// stringOpt() returns true if opt contains the named option and the option value is
// provided as a string

bool stringOpt(const StringVector& opt, const String& optname, String& optvalue)
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

bool intOpt(const StringVector& opt, const String& optname, int& optvalue)
{
   if (opt[0] == "-" + optname)
   {
      optvalue = stringToInt(opt[1]);
      return true;
   }
   else
      return false;
}

//------------------------------------------------------------------------------------
// doubleOpt() returns true if opt contains the named option and the option value is
// provided as a double

bool doubleOpt(const StringVector& opt, const String& optname, double& optvalue)
{
   if (opt[0] == "-" + optname)
   {
      optvalue = stringToDouble(opt[1]);
      return true;
   }
   else
      return false;
}

//------------------------------------------------------------------------------------
// intToString() returns a string representation of the given integer

String intToString(const int i)
{
   char buffer[100];

   std::snprintf(buffer, sizeof(buffer), "%d", i);

   return buffer;
}

//------------------------------------------------------------------------------------
// doubleToString() returns a string representation of the given double

String doubleToString(const double d)
{
   char buffer[100];

   std::snprintf(buffer, sizeof(buffer), "%.1f", d);

   return buffer;
}

//------------------------------------------------------------------------------------
// hasPrefix() returns true if the given string has the given prefix

bool hasPrefix(const String& s, const String& prefix)
{
   const int prefixLen = prefix.length();

   if (s.length() < prefixLen)
      return false;

   return (s.substr(0, prefixLen) == prefix);
}

//------------------------------------------------------------------------------------
// logOpen() opens the named log file

void logOpen(const String& filename)
{
   if (logfile)
      logClose();

   logfile = new std::ofstream(filename.c_str());
   if (!logfile->is_open())
      throw std::runtime_error("unable to open " + filename);
}

//------------------------------------------------------------------------------------
// logging() returns true if the log file is open

bool logging()
{
   return (logfile ? true : false);
}

//------------------------------------------------------------------------------------
// logWrite() writes a message to the log file

void logWrite(const String& message)
{
   if (!logfile) // log file is not open
      return;

   logMutex.lock();
   *logfile << message;
   logMutex.unlock();
}

//------------------------------------------------------------------------------------
// logClose() closes the log file

void logClose()
{
   if (logfile)
   {
      logfile->close();
      delete logfile;
      logfile = nullptr;
   }
}
