//------------------------------------------------------------------------------------
//
// fuzzion2html.cpp - this program reads fuzzion2 hits from stdin and writes them in
//                    HTML format to stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2021 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "util.h"
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

const std::string VERSION_ID     = "fuzzion2html v1.0.2, copyright 2021 "
                                   "St. Jude Children's Research Hospital";

const std::string FUZZION2       = "fuzzion2 ";
const std::string PATTERN        = "pattern ";
const std::string READ           = "read ";
const std::string READ_PAIRS     = "read-pairs ";

const std::string SEQUENCE       = "sequence";       const int SEQUENCE_COL = 1;
const std::string MBASES         = "matching bases"; const int MBASES_COL   = 2;
const std::string POSSIBLE       = "possible";       const int POSSIBLE_COL = 3;
const std::string PERCENT        = "% match";        const int PERCENT_COL  = 4;
const std::string ISIZE          = "insert size";    const int ISIZE_COL    = 5;

const int MIN_HEADING_COLS       = 6;
const int READ_COLS              = 5;

const int SECTION_INDENT         = 2;
const int ANNOTATION_INDENT      = 11;

const std::string BLANK          = "&nbsp;"; // non-breaking space character
const std::string DELETE         = "-";      // represents a deleted character

const std::string DELIM_COLOR    = "color:darkred";
const std::string MISMATCH_COLOR = "background-color:cyan";
const std::string INSERT_COLOR   = "background-color:yellow";
const std::string DELETE_COLOR   = "background-color:lime";

std::string title = ""; // optional title

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_ID << NEWLINE << NEWLINE
      << "Usage: " << progname << " OPTION ... < fuzzion2_hits > html" << NEWLINE;

   std::cerr
      << NEWLINE
      << "The following is optional:" << NEWLINE
      << "  -title=string   "
             << "string to include in the title of the HTML page" << NEWLINE;
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

      if (stringOpt(opt, "title", title))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   return true;
}

//------------------------------------------------------------------------------------

class MatchPattern
{
public:
   MatchPattern(const std::string& inPercent, const std::string& inSequence,
                const std::string& inInsertSize);

   virtual ~MatchPattern() { }

   std::string highlightSequence(int maxPatternLen) const;

   std::string percent, sequence, insertSize;
   int delim1, delim2; // index into sequence of brackets or braces
};

//------------------------------------------------------------------------------------

class MatchRead
{
public:
   MatchRead(const std::string& inPercent, const std::string& inSequence,
             const std::string& inLength,  const std::string& inName)
      : percent(inPercent), sequence(inSequence), length(inLength),
	name(inName) { }

   virtual ~MatchRead() { }

   std::string highlightSequence(const MatchPattern& pattern) const;

   std::string percent, sequence, length, name;
};

//------------------------------------------------------------------------------------

class Match
{
public:
   Match(const MatchPattern& inPattern, const MatchRead& inRead1,
         const MatchRead& inRead2)
      : pattern(inPattern), read1(inRead1), read2(inRead2) { }

   virtual ~Match() { }

   void write(int number, int maxPatternLen) const;

   MatchPattern pattern;
   MatchRead    read1, read2;
};

typedef std::vector<Match *> MatchVector;

struct MatchCompare
{
   bool operator()(Match* const& a, Match* const& b) const
   {
      // order by ascending delim1, then by ascending read name
      return (a->pattern.delim1 == b->pattern.delim1 ?
              (a->read1.name < b->read1.name) :
	      (a->pattern.delim1 < b->pattern.delim1));
   }
};

//------------------------------------------------------------------------------------

class Pattern
{
public:
   Pattern(const std::string& inName, const StringVector& inAnnotation)
      : name(inName), annotation(inAnnotation), match() { }

   virtual ~Pattern();

   void write(const StringVector& annotationHeading) const;

   std::string  name;
   StringVector annotation;
   MatchVector  match;
};

typedef std::map<std::string, Pattern> PatternMap; // key is pattern name

//------------------------------------------------------------------------------------
// openTag() returns a string that "opens" an HTML tag, e.g., "<b>"

std::string openTag(std::string tag, std::string style="")
{
   if (style == "")
      return "<" + tag + ">";
   else
      return "<" + tag + " style=\"" + style + "\">";
}

//------------------------------------------------------------------------------------
// closeTag() returns a string that "closes" an HTML tag, e.g., "</b>"

std::string closeTag(std::string tag)
{
   return "</" + tag + ">";
}

//------------------------------------------------------------------------------------
// wrap() returns a string that wraps opening and closing tags around a string

std::string wrap(std::string content, std::string tag, std::string style="")
{
   return openTag(tag, style) + content + closeTag(tag);
}

//------------------------------------------------------------------------------------
// writeHtmlBegin() writes to stdout the opening lines of the HTML document

void writeHtmlBegin(const std::string& version, uint64_t numReadPairs,
                    uint64_t numMatches, int numPatterns)
{
   const std::string TITLE =
      version + " results" + (title != "" ? " : " : "") + title;

   std::cout << openTag("!DOCTYPE html") << NEWLINE
             << openTag("html")   << NEWLINE
	     << openTag("head")   << NEWLINE
	     << wrap(TITLE, "title") << NEWLINE
	     << openTag("style")  << NEWLINE
	     << "table { color:black; background-color:ghostwhite; "
                "font-family:'Lucida Console', monospace; }" << NEWLINE
	     << closeTag("style") << NEWLINE
	     << closeTag("head")  << NEWLINE
	     << openTag("body")   << NEWLINE
	     << openTag("main", "font-family:arial") << NEWLINE
	     << wrap(TITLE, "h2") << NEWLINE
	     << openTag("p") << numReadPairs << " read pairs processed; ";

   if (numMatches == 0)
      std::cout << "no matches";
   else if (numMatches == 1)
      std::cout << "1 read pair matches 1 pattern";
   else
      std::cout << numMatches << " read pairs match " << numPatterns
                << (numPatterns == 1 ? " pattern" : " patterns");

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeHtmlEnd() writes to stdout the closing lines of the HTML document

void writeHtmlEnd()
{
   std::cout << closeTag("main") << NEWLINE
             << closeTag("body") << NEWLINE
	     << closeTag("html") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writePatternBegin() writes to stdout the opening lines of a pattern that has at
// least one match

void writePatternBegin(const std::string& patternName, int numMatches)
{
   std::cout << openTag("p") << NEWLINE
             << openTag("details") << NEWLINE
	     << openTag("summary")
	     << openTag("b")
	     << "pattern"
	     << "<a id=\"" << patternName << "\">" << BLANK << closeTag("a")
	     << patternName << " has " << numMatches
	     << " matching read " << (numMatches == 1 ? "pair" : "pairs")
	     << closeTag("b")
	     << closeTag("summary") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writePatternEnd() writes to stdout the closing line of a pattern that has at least
// one match

void writePatternEnd()
{
   std::cout << closeTag("details") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeAnnotationSectionBegin() writes to stdout the opening lines of an annotation
// section

void writeAnnotationSectionBegin()
{
   std::cout << openTag("details", "color:darkred") << NEWLINE
             << openTag("summary");

   for (int i = 1; i <= SECTION_INDENT; i++)
      std::cout << BLANK;

   std::cout << "pattern annotations" << closeTag("summary") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeAnnotationSectionEnd() writes to stdout the closing line of an annotation
// section

void writeAnnotationSectionEnd()
{
   std::cout << closeTag("details") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeAnnotation() writes an annotation to stdout

void writeAnnotation(const std::string& name, const std::string& value)
{
   for (int i = 1; i <= ANNOTATION_INDENT; i++)
      std::cout << BLANK;

   std::cout << wrap(name, "i") << " : " << value << openTag("br") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeMatchSectionBegin() writes to stdout the opening lines of a match section

void writeMatchSectionBegin()
{
   std::cout << openTag("details", "color:darkblue") << NEWLINE
             << openTag("summary");

   for (int i = 1; i <= SECTION_INDENT; i++)
      std::cout << BLANK;

   std::cout << "matching read pairs" << closeTag("summary") << NEWLINE
             << openTag("table") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeMatchSectionEnd() writes to stdout the closing lines of a match section

void writeMatchSectionEnd()
{
   std::cout << closeTag("table")   << NEWLINE
             << closeTag("details") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeMatchRow() writes to stdout one row of a match; if number is positive, it is a
// pattern row; otherwise, it is a read row

void writeMatchRow(int number, const std::string& matchPercent,
                   const std::string& sequence, const std::string& length,
		   std::string name="")
{
   std::cout << openTag("tr") << NEWLINE;

   if (number > 0)
      std::cout << openTag("td", "text-align:right")
	        << openTag("b") << number << closeTag("b")
		<< closeTag("td") << NEWLINE;
   else
      std::cout << openTag("td") << closeTag("td") << NEWLINE;

   std::cout << openTag("td", "text-align:right")
	     << openTag("nobr")
	     << BLANK << (matchPercent == "0.0" ? "N/A" : matchPercent) + BLANK
	     << closeTag("nobr")
             << closeTag("td") << NEWLINE;

   std::cout << openTag("td")
             << openTag("nobr") << sequence << closeTag("nobr")
	     << closeTag("td") << NEWLINE;

   std::cout << openTag("td", "text-align:right")
	     << openTag("nobr")
             << BLANK << (number > 0 ? "isize=" : "length=") << length
	     << closeTag("nobr")
	     << closeTag("td") << NEWLINE;

   if (name == "")
      std::cout << openTag("td") << closeTag("td") << NEWLINE;
   else
      std::cout << openTag("td")
	        << openTag("nobr")
	        << BLANK << name
	        << closeTag("nobr")
		<< closeTag("td") << NEWLINE;

   std::cout << closeTag("tr") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeBlankRow() writes to stdout a blank row

void writeBlankRow()
{
   const int NUM_COLS = 5;

   std::cout << openTag("tr") << NEWLINE;

   for (int i = 1; i <= NUM_COLS; i++)
      std::cout << openTag("td") << BLANK << closeTag("td") << NEWLINE;

   std::cout << closeTag("tr") << NEWLINE;
}

//------------------------------------------------------------------------------------
// MatchPattern::MatchPattern() looks for delimiter characters in the given sequence
// and saves their offsets

MatchPattern::MatchPattern(const std::string& inPercent,
                           const std::string& inSequence,
			   const std::string& inInsertSize)
   : percent(inPercent), sequence(inSequence), insertSize(inInsertSize)
{
   if ((delim1 = sequence.find(']'))             != std::string::npos &&
       (delim2 = sequence.find('[', delim1 + 1)) != std::string::npos)
      return; // this is a gene-fusion pattern

   if ((delim1 = sequence.find('}'))             != std::string::npos &&
       (delim2 = sequence.find('{', delim1 + 1)) != std::string::npos)
      return; // this is an ITD pattern

   throw std::runtime_error("missing delimiter in " + sequence);
}

//------------------------------------------------------------------------------------
// MatchPattern::highlightSequence() returns a highlighted representation of the
// pattern sequence

std::string MatchPattern::highlightSequence(int maxPatternLen) const
{
   int numTrailingBlanks = maxPatternLen - sequence.length();
   std::string trailingBlanks = "";

   for (int i = 1; i <= numTrailingBlanks; i++)
      trailingBlanks += BLANK;

   return wrap(sequence.substr(0, delim1), "span", "background-color:#ffe0b0") +
          wrap(wrap(sequence.substr(delim1, 1), "b"), "span", DELIM_COLOR) +
	  sequence.substr(delim1 + 1, delim2 - delim1 - 1) +
          wrap(wrap(sequence.substr(delim2, 1), "b"), "span", DELIM_COLOR) +
          wrap(sequence.substr(delim2 + 1) + trailingBlanks, "span",
               "background-color:#ffecbc");
}

//------------------------------------------------------------------------------------
// isDelimiter() returns true if the given character is a bracket or brace

bool isDelimiter(char ch)
{
   return (ch == '[' || ch == ']' || ch == '{' || ch == '}');
}

//------------------------------------------------------------------------------------
// highlightRead() finds a longest common subsequence of two strings and uses it to
// produce a highlighted version of the second string; the first string is a pattern
// sequence that may contain delimiter characters (brackets or braces)

std::string highlightRead(const std::string& strA, const std::string& strB)
{
   int lenA = strA.length();
   int lenB = strB.length();

   const char *a = strA.c_str();
   const char *b = strB.c_str();

   int **c = new int *[lenA + 1];

   for (int i = 0; i <= lenA; i++)
      c[i] = new int[lenB + 1]();

   for (int i = 0; i < lenA; i++)
      for (int j = 0; j < lenB; j++)
         if (a[i] == b[j])
            c[i + 1][j + 1] = c[i][j] + 1;
         else
            c[i + 1][j + 1] = std::max(c[i + 1][j], c[i][j + 1]);

   std::string highlight = "";

   int i = lenA, j = lenB;

   while (i > 0 || j > 0)
      if (i > 0 && j > 0 && a[i - 1] == b[j - 1]) // match
      {
         highlight = b[j - 1] + highlight;
	 i--;
	 j--;
      }
      else if (i > 0 && j > 0 && !isDelimiter(a[i - 1]) &&
               c[i - 1][j - 1] == c[i][j - 1] &&
	       c[i - 1][j - 1] == c[i - 1][j]) // mismatch
      {
         highlight =
            wrap(std::string(1, b[j - 1]), "span", MISMATCH_COLOR) + highlight;
	 i--;
	 j--;
      }
      else if (j > 0 && (i == 0 || c[i][j - 1] >= c[i - 1][j])) // insertion
      {
         if (i == lenA) // don't highlight if on last row
            highlight = b[j - 1] + highlight;
	 else
            highlight =
               wrap(std::string(1, b[j - 1]), "span", INSERT_COLOR) + highlight;
	 j--;
      }
      else // i > 0
      {
         if (isDelimiter(a[i - 1]))
            highlight =
               wrap(wrap(std::string(1, a[i - 1]), "b"), "span", DELIM_COLOR) +
	       highlight;
	 else // deletion
            highlight = wrap(DELETE, "span", DELETE_COLOR) + highlight;
	 i--;
      }

   for (int i = 0; i <= lenA; i++)
      delete c[i];

   delete c;

   return highlight;
}

//------------------------------------------------------------------------------------
// MatchRead::highlightSequence() returns a highlighted representation of a read
// sequence that matches the given pattern

std::string MatchRead::highlightSequence(const MatchPattern& pattern) const
{
   if (percent == "0.0") // unmatched mate
      return wrap(sequence, "span", "background-color:#c8c8cf");

   int patternLen = pattern.sequence.length();
   int readLen    = sequence.length();

   std::string leadingBlanks = "";
   int leadingBlankCount     = 0;

   for (int i = 0; i < readLen && sequence[i] == ' '; i++)
   {
      leadingBlanks += BLANK;
      leadingBlankCount++;
   }

   int beginOffset = leadingBlankCount;
   int endOffset   = readLen - 1;

   if (pattern.delim1 < beginOffset)
   {
      if (pattern.delim2 >= beginOffset && pattern.delim2 < readLen)
         endOffset++; // account for one delimiter
   }
   else if (pattern.delim2 <= readLen)
      endOffset += 2; // account for two delimiters
   else if (pattern.delim1 <  readLen)
      endOffset++;    // account for one delimiter

   endOffset = std::min(endOffset, patternLen - 1);

   if (beginOffset > endOffset)
      throw std::runtime_error("truncated pattern sequence: " + pattern.sequence);

   std::string patternSeq =
      pattern.sequence.substr(beginOffset, endOffset - beginOffset + 1);

   std::string readSeq = sequence.substr(beginOffset);

   return leadingBlanks + highlightRead(patternSeq, readSeq);
}

//------------------------------------------------------------------------------------
// Match::write() writes to stdout four rows representing a match: one pattern row
// followed by two read rows representing the matching read pair and one blank row
// for spacing

void Match::write(int number, int maxPatternLen) const
{
   std::string highlightPattern = pattern.highlightSequence(maxPatternLen);
   std::string highlightRead1   = read1.highlightSequence(pattern);
   std::string highlightRead2   = read2.highlightSequence(pattern);

   writeMatchRow(number, pattern.percent, highlightPattern, pattern.insertSize);
   writeMatchRow(0, read1.percent, highlightRead1, read1.length, read1.name);
   writeMatchRow(0, read2.percent, highlightRead2, read2.length, read2.name);
   writeBlankRow();
}

//------------------------------------------------------------------------------------
// Pattern::~Pattern() de-allocates each match of the pattern

Pattern::~Pattern()
{
   int numMatches = match.size();

   for (int i = 0; i < numMatches; i++)
      delete match[i];
}

//------------------------------------------------------------------------------------
// Pattern::write() writes to stdout all information about a pattern that has at least
// one match

void Pattern::write(const StringVector& annotationHeading) const
{
   int numAnnotations = annotation.size();
   int numMatches     = match.size();

   writePatternBegin(name, numMatches);

   if (numAnnotations > 0)
   {
      writeAnnotationSectionBegin();

      for (int i = 0; i < numAnnotations; i++)
         if (annotation[i] != "")
            writeAnnotation(annotationHeading[i], annotation[i]);

      writeAnnotationSectionEnd();
   }

   writeMatchSectionBegin();

   int maxPatternLen = 0;

   for (int i = 0; i < numMatches; i++)
      maxPatternLen = std::max(maxPatternLen,
		               static_cast<int>(match[i]->pattern.sequence.length()));

   for (int i = 0; i < numMatches; i++)
      match[i]->write(i + 1, maxPatternLen);

   writeMatchSectionEnd();

   writePatternEnd();
}

//------------------------------------------------------------------------------------
// readMatches() reads the hits from stdin and stores them in the map; the number of
// read pairs is summed and returned to the caller

uint64_t readMatches(PatternMap& pmap, std::string& version,
                     StringVector& annotationHeading)
{
   uint64_t numReadPairs = 0;
   std::string headingLine, line1, line2, line3;
   StringVector headingCol;

   if (!getline(std::cin, headingLine))
      throw std::runtime_error("no input");

   int numCols = splitString(headingLine, headingCol);

   if (numCols < MIN_HEADING_COLS ||
       !hasPrefix(headingCol[0], FUZZION2)  ||
       headingCol[SEQUENCE_COL] != SEQUENCE ||
       headingCol[MBASES_COL]   != MBASES   ||
       headingCol[POSSIBLE_COL] != POSSIBLE ||
       headingCol[PERCENT_COL]  != PERCENT  ||
       headingCol[ISIZE_COL]    != ISIZE)
      throw std::runtime_error("unexpected heading line");

   version = headingCol[0];

   for (int i = MIN_HEADING_COLS; i < numCols; i++)
      annotationHeading.push_back(headingCol[i]);

   while (getline(std::cin, line1))
   {
      if (hasPrefix(line1, FUZZION2)) // found another heading line
         if (line1 == headingLine)
            continue; // ignore it
         else
            throw std::runtime_error("inconsistent heading lines");

      StringVector col1, part1, col2, part2, col3, part3;

      if (hasPrefix(line1, READ_PAIRS))
         if (splitString(line1, col1) == 1 && splitString(col1[0], part1, ' ') == 2)
	 {
            std::istringstream stream(part1[1]);
	    uint64_t count = 0;
	    stream >> count;
	    numReadPairs += count;
	    continue;
	 }
	 else
            throw std::runtime_error("unexpected line: " + line1);

      if (!hasPrefix(line1, PATTERN) ||
          splitString(line1, col1) != numCols ||
          splitString(col1[0], part1, ' ') != 2 || part1[1] == "" ||
          !getline(std::cin, line2) || !hasPrefix(line2, READ) ||
	  splitString(line2, col2) != READ_COLS ||
	  splitString(col2[0], part2, ' ') != 2 || part2[1] == "" ||
          !getline(std::cin, line3) || !hasPrefix(line3, READ) ||
	  splitString(line3, col3) != READ_COLS ||
	  splitString(col3[0], part3, ' ') != 2 || part3[1] == "")
         throw std::runtime_error("invalid match format in: " + line1);

      const std::string& patternName = part1[1];

      PatternMap::iterator ppos = pmap.find(patternName);

      if (ppos == pmap.end())
      {
         StringVector annotation;

	 for (int i = MIN_HEADING_COLS; i < numCols; i++)
            annotation.push_back(col1[i]);

	 pmap.insert(std::make_pair(patternName, Pattern(patternName, annotation)));

	 ppos = pmap.find(patternName);
      }

      Pattern& pattern = ppos->second;

      if (pattern.match.size() >= std::numeric_limits<int>::max())
         continue; // too many matches; don't show them all

      MatchPattern mpattern(col1[PERCENT_COL], col1[SEQUENCE_COL], col1[ISIZE_COL]);

      MatchRead read1(col2[PERCENT_COL], col2[SEQUENCE_COL], col2[POSSIBLE_COL],
                      part2[1]);
      MatchRead read2(col3[PERCENT_COL], col3[SEQUENCE_COL], col3[POSSIBLE_COL],
                      part3[1]);

      Match *match = new Match(mpattern, read1, read2);

      pattern.match.push_back(match);
   }

   return numReadPairs;
}

//------------------------------------------------------------------------------------
// writeHtml() writes an HTML document to stdout

void writeHtml(PatternMap& pmap, const std::string& version,
               const StringVector& annotationHeading, uint64_t numReadPairs)
{
   // first, count the total number of matching read pairs and sort the matches
   uint64_t numMatches = 0;

   for (PatternMap::iterator ppos = pmap.begin(); ppos != pmap.end(); ++ppos)
   {
      Pattern& pattern = ppos->second;
      numMatches += pattern.match.size();
      std::sort(pattern.match.begin(), pattern.match.end(), MatchCompare());
   }

   // now write the HTML
   writeHtmlBegin(version, numReadPairs, numMatches, pmap.size());

   for (PatternMap::iterator ppos = pmap.begin(); ppos != pmap.end(); ++ppos)
   {
      const Pattern& pattern = ppos->second;
      pattern.write(annotationHeading);
   }

   writeHtmlEnd();
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
      PatternMap   pmap;
      std::string  version;
      StringVector annotationHeading;

      uint64_t numReadPairs = readMatches(pmap, version, annotationHeading);

      writeHtml(pmap, version, annotationHeading, numReadPairs);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
