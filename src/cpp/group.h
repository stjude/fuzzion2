//------------------------------------------------------------------------------------
//
// group.h - module with support for pattern groups
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2023 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef GROUP_H
#define GROUP_H

#include "summary.h"
#include <map>

typedef std::map<std::string, HitVector> ReadMap; // key is read1 name

//------------------------------------------------------------------------------------

class Group // represents a pattern group
{
public:
   Group(const std::string& inName, const StringVector& inAnnotation)
      : name(inName), annotation(inAnnotation), rmap() { }

   virtual ~Group() { }

   void addHit(Hit *hit);

   Summary *summarize(int minStrong, std::string sampleID="");

   int maxGroupDisplayLength();

   std::string name;        // name of the pattern group
   StringVector annotation; // group annotations
   ReadMap rmap;            // hits of patterns in the group, by read1 name
};

typedef std::map<std::string, Group> GroupMap; // key is group name

//------------------------------------------------------------------------------------

class GroupManager
{
public:
   GroupManager(const std::string& groupColList,
                const StringVector& patternAnnotationHeading,
		const HitVector& hitVector);

   virtual ~GroupManager() { }

   int readPairCount();

   StringVector annotationHeading; // group annotation headings
   GroupMap gmap;                  // groups with hits, by group name
};

#endif
