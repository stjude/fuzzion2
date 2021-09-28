//------------------------------------------------------------------------------------
//
// util.h - module defining some utilities
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>

inline int roundit(double d) { return static_cast<int>(d + 0.5); }

typedef std::vector<std::string> StringVector;

const char TAB     = '\t';
const char NEWLINE = '\n';
const char CRETURN = '\r';

bool getline(std::istream& stream, std::string& line);
int  splitString(const std::string& s, StringVector& v, char delimiter=TAB);

int    stringToNonnegInt(const std::string& s);
double stringToNonnegDouble(const std::string& s);

bool stringOpt(const StringVector& opt, std::string optname, std::string& optvalue);
bool intOpt   (const StringVector& opt, std::string optname, int& optvalue);
bool doubleOpt(const StringVector& opt, std::string optname, double& optvalue);

std::string intToString(int i);
std::string doubleToString(double d);

std::string intToStringLeadingZeros(int i, int width);

bool hasPrefix(const std::string& s, const std::string& prefix);

#endif
