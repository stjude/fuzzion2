//------------------------------------------------------------------------------------
//
// util.h - module defining utilities
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

typedef std::string String;
typedef std::vector<String> StringVector;
typedef std::vector<bool>   BoolVector;
typedef std::vector<int>    IntVector;
typedef std::runtime_error  Error;

const char TAB     = '\t';
const char NEWLINE = '\n';
const char CRETURN = '\r';

bool getline(std::istream& stream, String& line);
int  splitString(const String& s, StringVector& v, char delimiter=TAB);

int    stringToInt(const String& s);
double stringToDouble(const String& s);

bool stringOpt(const StringVector& opt, const String& optname, String& optvalue);
bool intOpt   (const StringVector& opt, const String& optname, int& optvalue);
bool doubleOpt(const StringVector& opt, const String& optname, double& optvalue);

String intToString(int i);
String doubleToString(double d);

bool hasPrefix(const String& s, const String& prefix);

void logOpen(const String& filename);
bool logging();
void logWrite(const String& message);
void logClose();

#endif
