//------------------------------------------------------------------------------------
//
// fuzzion2html.cpp - this program reads fuzzion2 hits from stdin and writes them in
//                    HTML format to stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "group.h"
#include "version.h"
#include <iostream>
#include <stdexcept>

const std::string VERSION_NAME    = "fuzzion2html " + CURRENT_VERSION;

const int SECTION_INDENT          = 2;
const int ANNOTATION_INDENT       = 11;

const std::string BLANK           = "&nbsp;";  // non-breaking space character
const std::string HYPHEN          = "&#8209;"; // non-breaking hyphen
const std::string DELETE          = "-";       // represents a deleted character
const std::string NA              = "N/A";     // not applicable

const std::string ALIGN_CENTER    = "text-align:center";
const std::string ALIGN_RIGHT     = "text-align:right";

const std::string NAME_COLOR      = "color:darkgreen";
const std::string DELIM_COLOR     = "color:darkred";
const std::string MISMATCH_COLOR  = "background-color:cyan";
const std::string INSERT_COLOR    = "background-color:yellow";
const std::string DELETE_COLOR    = "background-color:lime";

const std::string STRONG_MATCH    = "strong";
const std::string WEAK_MATCH      = "weak";
const std::string DUPLICATE_MATCH = "dup";

int minStrong = DEFAULT_MIN_STRONG; // minimum overlap for a strong match

std::string title = "";             // optional title
std::string groupColList = "";      // comma-separated list of group column headings

//------------------------------------------------------------------------------------
// showUsage() writes the program's usage to stderr

void showUsage(const char *progname)
{
   std::cerr
      << VERSION_NAME << ", " << COPYRIGHT << NEWLINE << NEWLINE
      << "Usage: " << progname << " OPTION ... < fuzzion2_hits > html"
      << NEWLINE;

   std::cerr
      << NEWLINE
      << "The following are optional:" << NEWLINE
      << "  -group=string   "
             << "comma-separated list of column headings, default is no grouping"
             << NEWLINE
      << "  -strong=N       "
             << "minimum overlap of a strong match in #bases, default is "
	     << DEFAULT_MIN_STRONG << NEWLINE
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

      if (intOpt   (opt, "strong", minStrong)    ||
          stringOpt(opt, "group",  groupColList) ||
          stringOpt(opt, "title",  title))
         continue;  // this option has been recognized

      return false; // unrecognized option
   }

   return (minStrong > 0);
}

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

void writeHtmlBegin(const std::string& fuzzion2Version, uint64_t numReadPairs,
                    int numMatches, int numMatched, bool grouping)
{
   std::string fullTitle = FUZZION2 + fuzzion2Version + " results";

   if (title != "")
      fullTitle += " : " + title;

   std::cout << openTag("!DOCTYPE html") << NEWLINE
             << openTag("html")   << NEWLINE
	     << openTag("head")   << NEWLINE
	     << wrap(fullTitle, "title") << NEWLINE
	     << openTag("style")  << NEWLINE
	     << "table { color:black; background-color:ghostwhite; "
                "font-family:'Lucida Console', monospace; }" << NEWLINE
	     << closeTag("style") << NEWLINE
	     << closeTag("head")  << NEWLINE
	     << openTag("body")   << NEWLINE
	     << openTag("main", "font-family:arial") << NEWLINE
	     << wrap(fullTitle, "h2") << NEWLINE
	     << openTag("p") << numReadPairs << " read pairs processed; ";

   if (numMatches == 0)
      std::cout << "no matches";
   else
   {
      if (numMatches == 1)
         std::cout << "1 read pair matches 1 ";
      else
         std::cout << numMatches << " read pairs match " << numMatched << " ";

      if (grouping)
         std::cout << (numMatched == 1 ? "pattern group" : "pattern groups");
      else
         std::cout << (numMatched == 1 ? "pattern" : "patterns");
   }

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
// writeSummaryBegin() writes to stdout the opening lines of a pattern or group

void writeSummaryBegin(const Summary *summary, bool grouping)
{
   std::cout << openTag("p") << NEWLINE
             << openTag("details") << NEWLINE
	     << openTag("summary")
	     << openTag("b")
	     << (grouping ? "group" : "pattern")
	     << "<a id=\"" << summary->name << "\">" << BLANK << closeTag("a")
	     << wrap(summary->name, "span", NAME_COLOR)
	     << " has " << summary->readPairs << " matching read "
	     << (summary->readPairs == 1 ? "pair" : "pairs")
	     << closeTag("b") << " ("
	     << summary->distinct()   << " distinct, "
	     << summary->weak         << " weak, "
	     << summary->strongNospan << " strong-, "
	     << summary->strongSpan   << " strong+)"
	     << closeTag("summary") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeSummaryEnd() writes to stdout the closing line of a pattern or group

void writeSummaryEnd()
{
   std::cout << closeTag("details") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeAnnotationSectionBegin() writes to stdout the opening lines of an annotation
// section

void writeAnnotationSectionBegin(bool grouping)
{
   std::cout << openTag("details", "color:darkred") << NEWLINE
             << openTag("summary");

   for (int i = 1; i <= SECTION_INDENT; i++)
      std::cout << BLANK;

   std::cout << (grouping ? "group" : "pattern") << " annotations"
             << closeTag("summary") << NEWLINE;
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
   if (value == "")
      return;

   for (int i = 1; i <= ANNOTATION_INDENT; i++)
      std::cout << BLANK;

   if (name != "")
      std::cout << wrap(name, "i") << " : ";
  
   std::cout << value << openTag("br") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeAllAnnotations() writes the given annotations to stdout

void writeAllAnnotations(const StringVector& annotationHeading,
                         const StringVector& annotation, bool grouping)
{
   int numAnnotationHeadings = annotationHeading.size();
   int numAnnotations        = annotation.size();

   bool hasValue = false;

   for (int i = 0; i < numAnnotations && !hasValue; i++)
      if (annotation[i] != "")
         hasValue = true;

   if (!hasValue)
      return; // there are no annotations

   writeAnnotationSectionBegin(grouping);

   const std::string NO_HEADING = "";

   for (int i = 0; i < numAnnotations; i++)
   {
      const std::string& heading =
         (i < numAnnotationHeadings ? annotationHeading[i] : NO_HEADING);

      writeAnnotation(heading, annotation[i]);
   }

   writeAnnotationSectionEnd();
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
// writeBlankColumn() writes to stdout a blank column

void writeBlankColumn()
{
   std::cout << openTag("td") << BLANK << closeTag("td") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeBlankRow() writes to stdout a blank row

void writeBlankRow()
{
   const int NUM_COLS = 6;

   std::cout << openTag("tr") << NEWLINE;

   for (int i = 1; i <= NUM_COLS; i++)
      writeBlankColumn();

   std::cout << closeTag("tr") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeMatchRow() writes to stdout one row of a match; if number is positive, it is a
// pattern row; otherwise, it is a read row

void writeMatchRow(std::string matchLabel, int number, int qualifier,
                   double percentMatch, const std::string& sequence, int length,
		   std::string name)
{
   std::cout << openTag("tr") << NEWLINE;

   // column 1 of 6
   std::cout << openTag("td", ALIGN_CENTER) << matchLabel
             << closeTag("td") << NEWLINE;

   // column 2 of 6
   if (number <= 0)
      writeBlankColumn();
   else
      std::cout << openTag("td", ALIGN_CENTER)
	        << BLANK << openTag("b")
		<< number << (qualifier > 0 ? HYPHEN + intToString(qualifier) : "")
	       	<< closeTag("b")
		<< closeTag("td") << NEWLINE;

   // column 3 of 6
   std::cout << openTag("td", ALIGN_RIGHT)
	     << BLANK << (percentMatch == 0.0 ? NA : doubleToString(percentMatch))
             << closeTag("td") << NEWLINE;

   // column 4 of 6
   std::cout << openTag("td")
             << openTag("nobr") << BLANK << sequence << closeTag("nobr")
	     << closeTag("td") << NEWLINE;

   // column 5 of 6
   std::cout << openTag("td", ALIGN_RIGHT)
	     << openTag("nobr")
             << BLANK << (number > 0 ? "isize=" : "length=") << length
	     << closeTag("nobr")
	     << closeTag("td") << NEWLINE;

   // column 6 of 6
   std::cout << openTag("td")
	     << openTag("nobr") << BLANK << name << closeTag("nobr")
             << closeTag("td") << NEWLINE;

   std::cout << closeTag("tr") << NEWLINE;
}

//------------------------------------------------------------------------------------
// highlightPatternSequence() returns a highlighted representation of a pattern
// sequence

std::string highlightPatternSequence(const std::string& sequence,
                                     int delim1, int delim2, int maxPatternLen)
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

inline bool isDelimiter(char ch)
{
   return (ch == '[' || ch == ']' || ch == '{' || ch == '}');
}

//------------------------------------------------------------------------------------
// highlight() finds a longest common subsequence of two strings and uses it to
// produce a highlighted version of the second string; the first string is a pattern
// sequence that may contain delimiter characters (brackets or braces); the second
// string is a read sequence

std::string highlight(const std::string& strA, const std::string& strB)
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

   std::string highlightString = "";

   int i = lenA, j = lenB;

   while (i > 0 || j > 0)
      if (i > 0 && j > 0 && a[i - 1] == b[j - 1]) // match
      {
         highlightString = b[j - 1] + highlightString;
	 i--;
	 j--;
      }
      else if (i > 0 && j > 0 && !isDelimiter(a[i - 1]) &&
               c[i - 1][j - 1] == c[i][j - 1] &&
	       c[i - 1][j - 1] == c[i - 1][j]) // mismatch
      {
         highlightString =
            wrap(std::string(1, b[j - 1]), "span", MISMATCH_COLOR) + highlightString;
	 i--;
	 j--;
      }
      else if (j > 0 && (i == 0 || c[i][j - 1] >= c[i - 1][j])) // insertion
      {
         if (i == lenA) // don't highlight if on last row
            highlightString = b[j - 1] + highlightString;
	 else
            highlightString =
               wrap(std::string(1, b[j - 1]), "span", INSERT_COLOR) + highlightString;
	 j--;
      }
      else // i > 0
      {
         if (isDelimiter(a[i - 1]))
            highlightString =
               wrap(wrap(std::string(1, a[i - 1]), "b"), "span", DELIM_COLOR) +
	       highlightString;
	 else // deletion
            highlightString = wrap(DELETE, "span", DELETE_COLOR) + highlightString;
	 i--;
      }

   for (int i = 0; i <= lenA; i++)
      delete c[i];

   delete c;

   return highlightString;
}

//------------------------------------------------------------------------------------
// highlightReadSequence() returns a highlighted representation of a read sequence
// that matches the given pattern sequence

std::string highlightReadSequence(int numLeadingBlanks,
                                  const std::string& readSequence,
				  double percentMatch,
				  const std::string& patternSequence,
				  int delim1, int delim2)
{
   std::string leadingBlanks = "";

   for (int i = 1; i <= numLeadingBlanks; i++)
      leadingBlanks += BLANK;

   if (percentMatch == 0.0) // unmatched mate
      return wrap(leadingBlanks + readSequence, "span", "background-color:#c8c8cf");

   int readLen     = readSequence.length();
   int patternLen  = patternSequence.length();

   int beginOffset = numLeadingBlanks;           // inclusive
   int endOffset   = numLeadingBlanks + readLen; // exclusive

   if (delim1 < beginOffset)
   {
      if (delim2 >= beginOffset && delim2 < endOffset)
         endOffset++; // account for one delimiter
   }
   else if (delim2 <= endOffset)
      endOffset += 2; // account for two delimiters
   else if (delim1 <  endOffset)
      endOffset++;    // account for one delimiter

   endOffset = std::min(endOffset, patternLen);

   if (beginOffset >= endOffset)
      throw std::runtime_error("truncated pattern sequence: " + patternSequence);

   std::string patternSubstr =
      patternSequence.substr(beginOffset, endOffset - beginOffset);

   return (leadingBlanks + highlight(patternSubstr, readSequence));
}

//------------------------------------------------------------------------------------
// writeMatch() writes to stdout four rows representing a match: one pattern row
// followed by two read rows representing the matching read pair and one blank row
// for spacing

void writeMatch(int number, int qualifier, const Hit& hit, int maxPatternLen)
{
   std::string highlightPattern =
      highlightPatternSequence(hit.pattern->displaySequence, hit.pattern->leftBases,
                               hit.pattern->delim2, maxPatternLen);

   std::string highlightRead1 =
      highlightReadSequence(hit.read1->leadingBlanks,  hit.read1->sequence,
                            hit.read1->percentMatch(), hit.pattern->displaySequence,
			    hit.pattern->leftBases, hit.pattern->delim2);

   std::string highlightRead2 =
      highlightReadSequence(hit.read2->leadingBlanks,  hit.read2->sequence,
                            hit.read2->percentMatch(), hit.pattern->displaySequence,
			    hit.pattern->leftBases, hit.pattern->delim2);

   writeMatchRow(hit.label(minStrong), number, qualifier, hit.pattern->percentMatch(),
                 highlightPattern, hit.pattern->insertSize, hit.pattern->name);

   writeMatchRow((hit.read1->isSpanning ? "+" : "-"), 0, 0, hit.read1->percentMatch(),
		 highlightRead1, hit.read1->sequence.length(), hit.read1->name);

   writeMatchRow((hit.read2->isSpanning ? "+" : "-"), 0, 0, hit.read2->percentMatch(),
		 highlightRead2, hit.read2->sequence.length(), hit.read2->name);

   writeBlankRow();
}

//------------------------------------------------------------------------------------
// writePattern() writes to stdout the annotations and matches of a pattern; the
// matches are given in hitVector from begin (inclusive) to end (exclusive)

void writePattern(const StringVector& annotationHeading, const HitVector& hitVector,
                  int begin, int end)
{
   Summary *summary = summarizeHits(hitVector, begin, end, minStrong);
   writeSummaryBegin(summary, false);
   delete summary;

   writeAllAnnotations(annotationHeading, hitVector[begin]->pattern->annotation,
                       false);

   writeMatchSectionBegin();

   int maxPatternLen = maxDisplayLength(hitVector, begin, end);

   for (int i = begin; i < end; i++)
   {
      const Hit *hit = hitVector[i];
      writeMatch(i - begin + 1, 0, *hit, maxPatternLen);
   }

   writeMatchSectionEnd();

   writeSummaryEnd();
}

//------------------------------------------------------------------------------------
// writeGroup() writes to stdout the annotations and matches of a group

void writeGroup(const StringVector& annotationHeading, Group& group)
{
   Summary *summary = group.summarize(minStrong);
   writeSummaryBegin(summary, true);
   delete summary;

   writeAllAnnotations(annotationHeading, group.annotation, true);

   writeMatchSectionBegin();

   int maxPatternLen = group.maxGroupDisplayLength();

   int number = 1;

   for (ReadMap::iterator rpos = group.rmap.begin(); rpos != group.rmap.end(); ++rpos)
   {
      const HitVector& hitVector = rpos->second;
      int numHits = hitVector.size();

      for (int i = 0; i < numHits; i++)
      {
         const Hit *hit = hitVector[i];
         writeMatch(number, (numHits > 1 ? i + 1 : 0), *hit, maxPatternLen);
      }

      number++;
   }

   writeMatchSectionEnd();

   writeSummaryEnd();
}

//------------------------------------------------------------------------------------
// writeAllPatterns() writes to stdout an HTML document displaying matches of patterns

void writeAllPatterns(const std::string& fuzzion2Version,
                      const StringVector& annotationHeading,
                      const HitVector& hitVector, uint64_t numReadPairs)
{
   IntVector index;
   getPatternIndices(hitVector, index);
   int numPatterns = index.size();

   writeHtmlBegin(fuzzion2Version, numReadPairs, hitVector.size(), numPatterns,
                  false);

   for (int i = 0; i < numPatterns; i++)
   {
      int begin = index[i];
      int end   = (i + 1 < numPatterns ? index[i + 1] : hitVector.size());

      writePattern(annotationHeading, hitVector, begin, end);
   }

   writeHtmlEnd();
}

//------------------------------------------------------------------------------------
// writeAllGroups() writes to stdout an HTML document displaying matches of groups

void writeAllGroups(const std::string& fuzzion2Version, uint64_t numReadPairs,
                    GroupManager& groupManager)
{
   writeHtmlBegin(fuzzion2Version, numReadPairs, groupManager.readPairCount(),
                  groupManager.gmap.size(), true);

   GroupMap& gmap = groupManager.gmap;

   for (GroupMap::iterator gpos = gmap.begin(); gpos != gmap.end(); ++gpos)
      writeGroup(groupManager.annotationHeading, gpos->second);

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
      std::string  fuzzion2Version;
      StringVector annotationHeading;
      HitVector    hitVector;

      uint64_t numReadPairs =
         readHits(std::cin, fuzzion2Version, annotationHeading, hitVector);

      if (groupColList == "")
         writeAllPatterns(fuzzion2Version, annotationHeading, hitVector,
                          numReadPairs);
      else
      {
         GroupManager groupManager(groupColList, annotationHeading, hitVector);
	 writeAllGroups(fuzzion2Version, numReadPairs, groupManager);
      }
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
