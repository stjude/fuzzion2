//------------------------------------------------------------------------------------
//
// group.cpp - module with support for pattern groups
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "group.h"

//------------------------------------------------------------------------------------
// Group::addHit() saves the given hit

void Group::addHit(Hit *hit)
{
   const String& readName1 = hit->read1->name;

   ReadMap::iterator rpos = rmap.find(readName1);
   if (rpos == rmap.end())
      rpos = rmap.insert(std::make_pair(readName1, HitVector())).first;

   HitVector& hitVector = rpos->second;
   hitVector.push_back(hit);
}

//------------------------------------------------------------------------------------
// Group::summarize() returns a pointer to a newly-allocated Summary of the group

Summary *Group::summarize(const int minStrong, const String& sampleID) const
{
   const int numMatches = rmap.size();
   int weak = 0, strongNospan = 0, strongSpan = 0;

   for (ReadMap::const_iterator rpos = rmap.cbegin(); rpos != rmap.cend(); ++rpos)
   {
      const HitVector& hitVector = rpos->second;
      const int numHits = hitVector.size();

      // find the "best" hit, assuming weak < strongNospan < strongSpan

      bool isWeak = false, isStrongNospan = false, isStrongSpan = false;

      for (int i = 0; i < numHits; i++)
      {
         const String label = hitVector[i]->label(minStrong);

	 if (label == HIT_STRONG_SPAN)
	 {
            isStrongSpan = true;
	    break;
	 }

	 if (label == HIT_STRONG_NOSPAN)
            isStrongNospan = true;
	 else if (label == HIT_WEAK)
            isWeak = true;
      }

      if (isStrongSpan)
         strongSpan++;
      else if (isStrongNospan)
         strongNospan++;
      else if (isWeak)
         weak++;
   }

   return new Summary(sampleID, numMatches, weak, strongNospan, strongSpan, name,
                      annotation);
}

//------------------------------------------------------------------------------------
// Group::maxGroupVisLength() returns the maximum length of the visualization sequence
// for the hits of the group

int Group::maxGroupVisLength() const
{
   int overallMaxlen = 0;

   for (ReadMap::const_iterator rpos = rmap.cbegin(); rpos != rmap.cend(); ++rpos)
   {
      const HitVector& hitVector = rpos->second;

      const int maxlen = maxVisLength(hitVector, 0, hitVector.size());
      if (maxlen > overallMaxlen)
         overallMaxlen = maxlen;
   }

   return overallMaxlen;
}

//------------------------------------------------------------------------------------
// GroupManager::GroupManager() parses a comma-separated list of column headings
// identifying the group key column first, followed by zero or more group annotation
// columns, and initializes the group map to contain the given hits

GroupManager::GroupManager(const String& groupColList,
                           const StringVector& patternAnnotationHeading,
			   const HitVector& hitVector)
   : annotationHeading(), gmap()
{
   StringVector groupCol;
   const int numGroupCols = splitString(groupColList, groupCol, ',');

   for (int i = 0; i < numGroupCols; i++)
      if (groupCol[i].find_first_not_of(' ') == String::npos)
         throw std::runtime_error("invalid group column list");

   int keyIndex;
   IntVector annotationIndex;

   const int numGroupAnnotations = numGroupCols - 1;
   int numPatternAnnotations = patternAnnotationHeading.size();

   for (int i = 0; i < numGroupCols; i++)
   {
      int j = 0;
      while (j < numPatternAnnotations && groupCol[i] != patternAnnotationHeading[j])
         j++;

      if (j < numPatternAnnotations)
         if (i == 0)
            keyIndex = j;
         else
	 {
            annotationIndex.push_back(j);
	    annotationHeading.push_back(groupCol[i]);
	 }
      else
         throw std::runtime_error("missing group column " + groupCol[i]);
   }

   const int numHits = hitVector.size();

   for (int i = 0; i < numHits; i++)
   {
      Hit *hit = hitVector[i];

      const StringVector& patternAnnotation = hit->annotation;
      numPatternAnnotations = patternAnnotation.size();

      if (keyIndex >= numPatternAnnotations)
         continue; // no group name

      String groupName = patternAnnotation[keyIndex];

      const std::size_t first = groupName.find_first_not_of(' ');
      if (first == String::npos)
         continue; // no group name

      const std::size_t last = groupName.find_last_not_of(' ');
      groupName = groupName.substr(first, last - first + 1); // trim blanks from name

      GroupMap::iterator gpos = gmap.find(groupName);

      if (gpos == gmap.end())
      {
         StringVector groupAnnotation;

	 for (int j = 0; j < numGroupAnnotations; j++)
	 {
            const int k = annotationIndex[j];

	    if (k < numPatternAnnotations)
               groupAnnotation.push_back(patternAnnotation[k]);
	    else
               groupAnnotation.push_back(" ");
	 }

	 gpos = gmap.insert(
             std::make_pair(groupName, Group(groupName, groupAnnotation))).first;
      }

      Group& group = gpos->second;
      group.addHit(hit);
   }
}

//------------------------------------------------------------------------------------
// GroupManager::readCount() returns the total number of reads in groups

int GroupManager::readCount() const
{
   int count = 0;

   for (GroupMap::const_iterator gpos = gmap.cbegin(); gpos != gmap.cend(); ++gpos)
   {
      const Group& group = gpos->second;
      count += group.rmap.size();
   }

   return count;
}
