//------------------------------------------------------------------------------------
//
// fuzzion2html.cpp - this program reads fuzzion2 hits from stdin and writes them in
//                    HTML format to stdout
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "group.h"
#include "version.h"

const String VERSION_NAME   = "fuzzion2html " + CURRENT_VERSION;

const int SECTION_INDENT    = 2;
const int ANNOTATION_INDENT = 11;

const String BLANK          = "&nbsp;";  // non-breaking space character
const String ELLIPSIS       = "&mldr;";  // an ellipsis ...
const String HYPHEN         = "&#8209;"; // non-breaking hyphen

const String ALIGN_CENTER   = "text-align:center";
const String ALIGN_RIGHT    = "text-align:right";

const String NAME_COLOR     = "color:darkgreen";
const String DELIM_COLOR    = "color:darkred";
const String LEFT_COLOR     = "background-color:#ffe0b0";
const String MIDDLE_COLOR   = "background-color:#f8ecd8";
const String RIGHT_COLOR    = "background-color:#ffecbc";
const String MISMATCH_COLOR = "background-color:cyan";
const String INSERT_COLOR   = "background-color:yellow";
const String DELETE_COLOR   = "background-color:lime";

int minStrong = DEFAULT_MIN_STRONG; // minimum overlap for a strong match

String title = "";         // optional title
String groupColList = "";  // comma-separated list of group column headings

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

bool parseArgs(const int argc, const char *argv[])
{
   for (int i = 1; i < argc; i++)
   {
      const String arg = argv[i];
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
         continue;     // this option has been recognized

      return false;    // unrecognized option
   }

   return (minStrong > 0);
}

//------------------------------------------------------------------------------------
// openTag() returns a string that "opens" an HTML tag, e.g., "<b>"

String openTag(const String& tag, const String& style="")
{
   if (style == "")
      return "<" + tag + ">";
   else
      return "<" + tag + " style=\"" + style + "\">";
}

//------------------------------------------------------------------------------------
// closeTag() returns a string that "closes" an HTML tag, e.g., "</b>"

String closeTag(const String& tag)
{
   return "</" + tag + ">";
}

//------------------------------------------------------------------------------------
// wrap() returns a string that wraps opening and closing tags around a string

String wrap(const String& content, const String& tag, const String& style="")
{
   return openTag(tag, style) + content + closeTag(tag);
}

//------------------------------------------------------------------------------------
// colorString() returns a string that applies a color to the given string

String colorString(const String& content, const String& color)
{
   return wrap(content, "span", color);
}

//------------------------------------------------------------------------------------
// colorChar() returns a string that applies a color to the given character

String colorChar(const char ch, const String& color)
{
   const char content[2] = { ch, 0 };
   return colorString(content, color);
}

//------------------------------------------------------------------------------------
// writeHtmlBegin() writes the opening lines of the HTML document

void writeHtmlBegin(const String& fuzzion2Version, const uint64_t numReads,
                    const int numHits, const int numMatched, const bool grouping)
{
   String fullTitle = FUZZION2 + fuzzion2Version + " results";

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
             << openTag("p") << numReads << " reads processed; ";

   if (numHits == 0)
      std::cout << "no matches";
   else
   {
      if (numHits == 1)
         std::cout << "1 match";
      else
         std::cout << numHits << " matches";

      std::cout << " to " << numMatched << " ";

      if (grouping)
         std::cout << (numMatched == 1 ? "pattern group" : "pattern groups");
      else
         std::cout << (numMatched == 1 ? "pattern" : "patterns");
   }

   std::cout << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeHtmlEnd() writes the closing lines of the HTML document

void writeHtmlEnd()
{
   std::cout << closeTag("main") << NEWLINE
             << closeTag("body") << NEWLINE
             << closeTag("html") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeSummaryBegin() writes the opening lines of a pattern or group

void writeSummaryBegin(const Summary *summary, const bool grouping)
{
   std::cout << openTag("p") << NEWLINE
             << openTag("details") << NEWLINE
             << openTag("summary")
             << openTag("b")
             << (grouping ? "group" : "pattern")
             << "<a id=\"" << summary->name << "\">" << BLANK << closeTag("a")
             << colorString(summary->name, NAME_COLOR)
             << " has " << summary->numMatches
             << (summary->numMatches == 1 ? " match" : " matches")
             << closeTag("b") << " ("
             << summary->distinct()   << " distinct, "
             << summary->weak         << " " << HIT_WEAK          << ", "
             << summary->strongNospan << " " << HIT_STRONG_NOSPAN << ", "
             << summary->strongSpan   << " " << HIT_STRONG_SPAN   << ")"
             << closeTag("summary") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeSummaryEnd() writes the closing line of a pattern or group

void writeSummaryEnd()
{
   std::cout << closeTag("details") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeAnnotationSectionBegin() writes the opening lines of an annotation section

void writeAnnotationSectionBegin(const bool grouping)
{
   std::cout << openTag("details", "color:darkred") << NEWLINE
             << openTag("summary");

   for (int i = 1; i <= SECTION_INDENT; i++)
      std::cout << BLANK;

   std::cout << (grouping ? "group" : "pattern") << " annotations"
             << closeTag("summary") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeAnnotationSectionEnd() writes the closing line of an annotation section

void writeAnnotationSectionEnd()
{
   std::cout << closeTag("details") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeAnnotation() writes an annotation

void writeAnnotation(const String& name, const String& value)
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
// writeAllAnnotations() writes the given annotations

void writeAllAnnotations(const StringVector& annotationHeading,
                         const StringVector& annotation, const bool grouping)
{
   const int numAnnotationHeadings = annotationHeading.size();
   const int numAnnotations        = annotation.size();

   bool hasValue = false;

   for (int i = 0; i < numAnnotations && !hasValue; i++)
      if (annotation[i] != "")
         hasValue = true;

   if (!hasValue)
      return; // there are no annotations

   writeAnnotationSectionBegin(grouping);

   const String NO_HEADING = "";

   for (int i = 0; i < numAnnotations; i++)
   {
      const String& heading =
         (i < numAnnotationHeadings ? annotationHeading[i] : NO_HEADING);

      writeAnnotation(heading, annotation[i]);
   }

   writeAnnotationSectionEnd();
}

//------------------------------------------------------------------------------------
// writeMatchSectionBegin() writes the opening lines of a match section

void writeMatchSectionBegin()
{
   std::cout << openTag("details", "color:darkblue") << NEWLINE
             << openTag("summary");

   for (int i = 1; i <= SECTION_INDENT; i++)
      std::cout << BLANK;

   std::cout << "matches" << closeTag("summary") << NEWLINE
             << openTag("table") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeMatchSectionEnd() writes the closing lines of a match section

void writeMatchSectionEnd()
{
   std::cout << closeTag("table")   << NEWLINE
             << closeTag("details") << NEWLINE;
}

//------------------------------------------------------------------------------------
// isLeftDelimiter() returns true if the given character is a left pattern delimiter

bool isLeftDelimiter(const char ch)
{
   return (ch == ']' || ch == '}' || ch == '<' || ch == '(');
}

//------------------------------------------------------------------------------------
// isRightDelimiter() returns true if the given character is a right pattern delimiter

bool isRightDelimiter(const char ch)
{
   return (ch == '[' || ch == '{' || ch == '>' || ch == ')');
}

//------------------------------------------------------------------------------------
// isInclusion() returns true if the given character identifies an inclusion pattern

bool isInclusion(const char ch)
{
   return (ch == '<' || ch == '>');
}

//------------------------------------------------------------------------------------
// isExclusion() returns true if the given character identifies an exclusion pattern

bool isExclusion(const char ch)
{
   return (ch == '(' || ch == ')');
}

//------------------------------------------------------------------------------------
// highlightDelimiter() returns a string that highlights the given delimiter character

String highlightDelimiter(const char ch)
{
   const char delimiter[2] = { ch, 0 };
   return colorString(wrap(delimiter, "b"), DELIM_COLOR);
}

//------------------------------------------------------------------------------------
// adjust() replaces blanks with HTML non-breaking blanks and periods with HTML
// ellipses

String adjust(const String& content)
{
   const int len = content.length();
   String adjustedContent = "";

   for (int i = 0; i < len; i++)
      if (content[i] == ' ')
         adjustedContent += BLANK;
      else if (content[i] == ELLIPSIS_MARK)
         adjustedContent += ELLIPSIS;
      else
         adjustedContent += content[i];

   return adjustedContent;
}

//------------------------------------------------------------------------------------
// highlightPattern() returns a highlighted pattern visualization padded with trailing
// blanks

String highlightPattern(const String& patternVis, const int maxVisLen,
                        int& delim1, int& delim2, bool& inclusion, bool& exclusion)
{
   const int len = patternVis.length();
   const int numTrailing = maxVisLen - len;
   String trailing = ""; // trailing blanks

   for (int i = 1; i <= numTrailing; i++)
      trailing += BLANK;

   delim1    = -1;       // index of left  delimiter or -1 if not found
   delim2    = -1;       // index of right delimiter or -1 if not found
   inclusion = false;    // will be set to true if this is an inclusion pattern
   exclusion = false;    // will be set to true if this is an exclusion pattern

   for (int i = 0; i < len; i++)
   {
      const char ch = patternVis[i];

      if (isLeftDelimiter(ch))
      {
         delim1    = i;
         inclusion = isInclusion(ch);
         exclusion = isExclusion(ch);
      }
      else if (isRightDelimiter(ch))
      {
         delim2    = i;
         inclusion = isInclusion(ch);
         exclusion = isExclusion(ch);
         break;
      }
   }

   if (delim1 >= 0 && delim2 >= delim1 + 1)
   {
      const int midlen = delim2 - delim1 - 1;

      return 
      (colorString(adjust(patternVis.substr(0, delim1)), LEFT_COLOR) +
       highlightDelimiter(patternVis[delim1]) +
       colorString(adjust(patternVis.substr(delim1 + 1, midlen)), MIDDLE_COLOR) +
       highlightDelimiter(patternVis[delim2]) +
       colorString(adjust(patternVis.substr(delim2 + 1)) + trailing, RIGHT_COLOR));
   }
   else // unsided pattern
      return colorString(adjust(patternVis) + trailing, LEFT_COLOR);
}

//------------------------------------------------------------------------------------
// highlightRead() returns a highlighted read visualization

String highlightRead(const String& patternVis, const String& readVis,
                     const int delim1, const int delim2,
                     const bool inclusion, const bool exclusion)
{
   const int len = patternVis.length();
   if (len != readVis.length())
      throw Error("hit has unexpected visualization length");

   const bool extra = (inclusion || exclusion);

   String highlight = "";

   for (int i = 0; i < len; i++)
   {
      const char p = patternVis[i];
      const char r = readVis[i];

      if (r == ' ')
         highlight += BLANK;
      else if (r == ELLIPSIS_MARK)
         highlight += ELLIPSIS;
      else if (p == ' ')
         highlight += r;
      else if (i == delim1 || i == delim2)
         highlight += highlightDelimiter(r);
      else if (p == SPACER)
         highlight += colorChar(r, INSERT_COLOR);
      else if (r == SPACER)
         highlight += colorChar(r, DELETE_COLOR);
      else if ((!extra || i < delim1 || i > delim2) && upperACGT(p) != upperACGT(r))
         highlight += colorChar(r, MISMATCH_COLOR);
      else
         highlight += r;
   }

   return highlight;
}

//------------------------------------------------------------------------------------
// writeBlankColumn() writes a blank column

void writeBlankColumn()
{
   std::cout << openTag("td") << BLANK << closeTag("td") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeBlankRow() writes a blank row

void writeBlankRow()
{
   const int NUM_COLS = 6;

   std::cout << openTag("tr") << NEWLINE;

   for (int i = 1; i <= NUM_COLS; i++)
      writeBlankColumn();

   std::cout << closeTag("tr") << NEWLINE;
}

//------------------------------------------------------------------------------------
// writeMatchRow() writes one row of a match; if number is positive, it is a pattern
// row; otherwise, it is a read row

void writeMatchRow(const String& matchLabel, const int number, const int qualifier,
                   const int matches, const int possible, const String& vis,
                   const int length, const String& name, bool showAgreement=true)
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
   if (showAgreement)
      std::cout << openTag("td", ALIGN_RIGHT)
                << BLANK << doubleToString(100.0 * matches / possible)
                << closeTag("td") << NEWLINE;
   else
      writeBlankColumn();

   // column 4 of 6
   std::cout << openTag("td")
             << openTag("nobr") << BLANK << vis << closeTag("nobr")
             << closeTag("td") << NEWLINE;

   // column 5 of 6
   if (length == 0)
      writeBlankColumn();
   else
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
// writeMatch() writes 3 or 4 rows representing a match: 1 pattern row followed by 1
// or 2 read rows representing the matching read(s) and 1 blank row for spacing

void writeMatch(const int number, const int qualifier, const Hit& hit,
                const int maxVisLen)
{
   int  delim1, delim2;
   bool inclusion, exclusion;

   String vis = highlightPattern(hit.patternVis, maxVisLen, delim1, delim2, inclusion,
                                 exclusion);

   writeMatchRow(hit.label(minStrong), number, qualifier, hit.matches, hit.possible,
                 vis, hit.insertSize, hit.patternName);

   vis = highlightRead(hit.patternVis, hit.read1->vis, delim1, delim2, inclusion,
                       exclusion);

   writeMatchRow(BLANK, 0, 0, hit.read1->matches, hit.read1->possible, vis,
                 hit.read1->len, hit.read1->name, (hit.read2 ? true : false));

   if (hit.read2)
   {
      vis = highlightRead(hit.patternVis, hit.read2->vis, delim1, delim2, inclusion,
                          exclusion);

      writeMatchRow(BLANK, 0, 0, hit.read2->matches, hit.read2->possible, vis,
                    hit.read2->len, hit.read2->name);
   }

   writeBlankRow();
}

//------------------------------------------------------------------------------------
// writePattern() writes the annotations and matches of a pattern; the matches are
// given in hitVector from begin (inclusive) to end (exclusive)

void writePattern(const StringVector& annotationHeading, const HitVector& hitVector,
                  const int begin, const int end)
{
   const Summary *summary = summarizeHits(hitVector, begin, end, minStrong);
   writeSummaryBegin(summary, false);
   delete summary;

   writeAllAnnotations(annotationHeading, hitVector[begin]->annotation, false);

   writeMatchSectionBegin();

   const int maxVisLen = maxVisLength(hitVector, begin, end);

   for (int i = begin; i < end; i++)
   {
      const Hit *hit = hitVector[i];
      writeMatch(i - begin + 1, 0, *hit, maxVisLen);
   }

   writeMatchSectionEnd();

   writeSummaryEnd();
}

//------------------------------------------------------------------------------------
// writeGroup() writes the annotations and matches of a group

void writeGroup(const StringVector& annotationHeading, const Group& group)
{
   const Summary *summary = group.summarize(minStrong);
   writeSummaryBegin(summary, true);
   delete summary;

   writeAllAnnotations(annotationHeading, group.annotation, true);

   writeMatchSectionBegin();

   const int maxVisLen = group.maxGroupVisLength();

   int number = 1;
   const ReadMap& rmap = group.rmap;

   for (ReadMap::const_iterator rpos = rmap.cbegin(); rpos != rmap.cend(); ++rpos)
   {
      const HitVector& hitVector = rpos->second;
      const int numHits = hitVector.size();

      for (int i = 0; i < numHits; i++)
      {
         const Hit *hit = hitVector[i];
         writeMatch(number, (numHits > 1 ? i + 1 : 0), *hit, maxVisLen);
      }

      number++;
   }

   writeMatchSectionEnd();

   writeSummaryEnd();
}

//------------------------------------------------------------------------------------
// writeAllPatterns() writes an HTML document displaying matches of patterns

void writeAllPatterns(const String& fuzzion2Version,
                      const StringVector& annotationHeading,
                      const HitVector& hitVector, const uint64_t numReads)
{
   IntVector index;
   getPatternIndices(hitVector, index);
   const int numPatterns = index.size();

   writeHtmlBegin(fuzzion2Version, numReads, hitVector.size(), numPatterns, false);

   for (int i = 0; i < numPatterns; i++)
   {
      const int begin = index[i];
      const int end   = (i + 1 < numPatterns ? index[i + 1] : hitVector.size());

      writePattern(annotationHeading, hitVector, begin, end);
   }

   writeHtmlEnd();
}

//------------------------------------------------------------------------------------
// writeAllGroups() writes an HTML document displaying matches of groups

void writeAllGroups(const String& fuzzion2Version, const GroupManager& groupManager,
                    const uint64_t numReads)
{
   writeHtmlBegin(fuzzion2Version, numReads, groupManager.readCount(),
                  groupManager.gmap.size(), true);

   const GroupMap& gmap = groupManager.gmap;

   for (GroupMap::const_iterator gpos = gmap.cbegin(); gpos != gmap.cend(); ++gpos)
      writeGroup(groupManager.annotationHeading, gpos->second);

   writeHtmlEnd();
}

//------------------------------------------------------------------------------------

int main(const int argc, const char *argv[])
{
   if (!parseArgs(argc, argv))
   {
      showUsage(argv[0]);
      return 1;
   }

   try
   {
      String       fuzzion2Version;
      StringVector annotationHeading;
      HitVector    hitVector;
      uint64_t     numReads;

      readHits(std::cin, fuzzion2Version, annotationHeading, hitVector, numReads);

      if (groupColList == "")
         writeAllPatterns(fuzzion2Version, annotationHeading, hitVector, numReads);
      else
      {
         GroupManager groupManager(groupColList, annotationHeading, hitVector);
         writeAllGroups(fuzzion2Version, groupManager, numReads);
      }
   }
   catch (const Error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
