//------------------------------------------------------------------------------------
//
// group.h - module with support for pattern groups
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2026 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef GROUP_H
#define GROUP_H

#include "summary.h"
#include <map>

typedef std::map<String, HitVector> ReadMap; // key is read1 name

//------------------------------------------------------------------------------------

class Group // represents a pattern group
{
public:
   Group(const String& inName, const StringVector& inAnnotation)
      : name(inName), annotation(inAnnotation), rmap() { }

   virtual ~Group() { } // the hits in rmap are not de-allocated

   void addHit(Hit *hit);

   Summary *summarize(int minStrong, const String& sampleID="") const;

   int maxGroupVisLength() const;

   const String name;             // name of the pattern group
   const StringVector annotation; // group annotations
   ReadMap rmap;                  // hits of patterns in the group, by read1 name
};

typedef std::map<String, Group> GroupMap; // key is group name

//------------------------------------------------------------------------------------

class GroupManager
{
public:
   GroupManager(const String& groupColList,
                const StringVector& patternAnnotationHeading,
		const HitVector& hitVector);

   virtual ~GroupManager() { }

   int readCount() const;

   StringVector annotationHeading; // group annotation headings
   GroupMap gmap;                  // groups with hits, by group name
};

#endif
