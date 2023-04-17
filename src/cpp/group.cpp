//------------------------------------------------------------------------------------
//
// group.cpp - module with support for pattern groups
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "group.h"
#include <stdexcept>

//------------------------------------------------------------------------------------
// Group::addHit() saves the given hit

void Group::addHit(Hit *hit)
{
   const std::string& read1Name = hit->read1->name;

   ReadMap::iterator rpos = rmap.find(read1Name);
   if (rpos == rmap.end())
   {
      rmap.insert(std::make_pair(read1Name, HitVector()));
      rpos = rmap.find(read1Name);
   }

   HitVector& hitVector = rpos->second;
   hitVector.push_back(hit);
}

//------------------------------------------------------------------------------------
// Group::summarize() returns a pointer to a newly-allocated Summary of the group

Summary *Group::summarize(int minStrong, std::string sampleID)
{
   int readPairs = rmap.size();
   int weak = 0, strongNospan = 0, strongSpan = 0;

   for (ReadMap::iterator rpos = rmap.begin(); rpos != rmap.end(); ++rpos)
   {
      const HitVector& hitVector = rpos->second;
      int numHits = hitVector.size();

      // find the "best" hit, assuming weak < strongNospan < strongSpan

      bool isWeak = false, isStrongNospan = false, isStrongSpan = false;

      for (int i = 0; i < numHits; i++)
      {
         std::string label = hitVector[i]->label(minStrong);

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

   return new Summary(sampleID, readPairs, weak, strongNospan, strongSpan, name,
                      annotation);
}

//------------------------------------------------------------------------------------
// Group::maxGroupDisplayLength() returns the maximum length of the pattern display
// sequence for the hits of the group

int Group::maxGroupDisplayLength()
{
   int overallMaxLength = 0;

   for (ReadMap::iterator rpos = rmap.begin(); rpos != rmap.end(); ++rpos)
   {
      const HitVector& hitVector = rpos->second;

      int maxLength = maxDisplayLength(hitVector, 0, hitVector.size());
      if (maxLength > overallMaxLength)
         overallMaxLength = maxLength;
   }

   return overallMaxLength;
}

//------------------------------------------------------------------------------------
// GroupManager::GroupManager() parses a comma-separated list of column headings
// identifying the group key column first, followed by zero or more group annotation
// columns, and initializes the group map to contain the given hits

GroupManager::GroupManager(const std::string& groupColList,
                           const StringVector& patternAnnotationHeading,
			   const HitVector& hitVector)
   : annotationHeading(), gmap()
{
   StringVector groupCol;
   int numGroupCols = splitString(groupColList, groupCol, ',');

   for (int i = 0; i < numGroupCols; i++)
      if (groupCol[i].find_first_not_of(' ') == std::string::npos)
         throw std::runtime_error("invalid group column list");

   int keyIndex;
   IntVector annotationIndex;

   int numGroupAnnotations   = numGroupCols - 1;
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

   int numHits = hitVector.size();

   for (int i = 0; i < numHits; i++)
   {
      Hit *hit = hitVector[i];

      const StringVector& patternAnnotation = hit->pattern->annotation;
      numPatternAnnotations = patternAnnotation.size();

      if (keyIndex >= numPatternAnnotations)
         continue; // no group name

      std::string groupName = patternAnnotation[keyIndex];

      int first = groupName.find_first_not_of(' ');
      if (first == std::string::npos)
         continue; // no group name

      int last  = groupName.find_last_not_of(' ');
      groupName = groupName.substr(first, last - first + 1); // trim blanks from name

      GroupMap::iterator gpos = gmap.find(groupName);

      if (gpos == gmap.end())
      {
         StringVector groupAnnotation;

	 for (int j = 0; j < numGroupAnnotations; j++)
	 {
            int k = annotationIndex[j];

	    if (k < numPatternAnnotations)
               groupAnnotation.push_back(patternAnnotation[k]);
	    else
               groupAnnotation.push_back(" ");
	 }

	 gmap.insert(std::make_pair(groupName, Group(groupName, groupAnnotation)));

	 gpos = gmap.find(groupName);
      }

      Group& group = gpos->second;
      group.addHit(hit);
   }
}

//------------------------------------------------------------------------------------
// GroupManager::readPairCount() returns the number of read pairs in the group

int GroupManager::readPairCount()
{
   int count = 0;

   for (GroupMap::iterator gpos = gmap.begin(); gpos != gmap.end(); ++gpos)
   {
      const Group& group = gpos->second;
      count += group.rmap.size();
   }

   return count;
}
