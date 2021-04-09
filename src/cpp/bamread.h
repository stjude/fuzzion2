//------------------------------------------------------------------------------------
//
// bamread.h - module providing functionality for reading Bam files
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef BAMREAD_H
#define BAMREAD_H

#include <string>

//------------------------------------------------------------------------------------

class BamRead // a single read from a Bam file
{
public:
   BamRead();

   virtual ~BamRead();

   std::string *name()         const;
   const char  *constName()    const;

   int  bitFlag()              const;
   bool isPaired()             const;
   bool isProperPair()         const;
   bool isUnmapped()           const;
   bool isMateUnmapped()       const;
   bool isReverseStrand()      const;
   bool isMateReverseStrand()  const;
   bool isRead1()              const;
   bool isRead2()              const;
   bool isSecondary()          const;
   bool isFailedQC()           const;
   bool isDuplicate()          const;
   bool isSupplementary()      const;

   int  refID()                const;
   int  position()             const; // 1-based, inclusive
   int  endPosition()          const; // 1-based, inclusive

   int  mapQuality()           const;

   int  numCigarOps()          const;
   void cigarOp(int index, char& op, int& len) const;

   int  mateRefID()            const;
   int  matePosition()         const; // 1-based, inclusive

   int  insertSize()           const;

   int  length()               const;

   std::string *sequence()     const;
   std::string *qualities()    const;

protected:
   void *aptr; // opaque pointer to internal data structure

   friend class BamReader;
};

//------------------------------------------------------------------------------------

class BamReader // used to get BamReads from a Bam file
{
public:
   BamReader();

   virtual ~BamReader();

   void open(const std::string& filename);

   int  getNumRef() const;

   int  getRefID(const std::string& refName);
   int  getRefIDAlt(const std::string& refName);

   int  getRefLen(int refID) const;
   std::string getRefName(int refID) const;

   void jump(int refID, int startPosition=1);

   bool getNext(BamRead& read);

   void close();

protected:
   void *rptr; // opaque pointer to internal data structure
};

//------------------------------------------------------------------------------------
#endif
