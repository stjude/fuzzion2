//------------------------------------------------------------------------------------
//
// bamread.cpp - module providing functionality for reading Bam files; htslib is used
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "bamread.h"
#include "htslib/sam.h"
#include <cstring>
#include <map>
#include <stdexcept>

#define bam1_ptr ((bam1_t *)aptr)
#define data_ptr ((InternalData *)rptr)

typedef std::map<std::string, int> RefMap; // refName, refID

//------------------------------------------------------------------------------------

class InternalData
{
public:
   InternalData()
      : filename(""), file(NULL), header(NULL), refmap(), index(NULL), iter(NULL) { }

   virtual ~InternalData() { reinitialize(); }

   void reinitialize();

   std::string  filename;
   htsFile     *file;
   bam_hdr_t   *header;
   RefMap       refmap;
   hts_idx_t   *index;
   hts_itr_t   *iter;
};

//------------------------------------------------------------------------------------
// InternalData::reinitialize() discards any existing data and reinitializes the
// fields

void InternalData::reinitialize()
{
   if (iter)   { hts_itr_destroy(iter);   iter   = NULL; }
   if (index)  { hts_idx_destroy(index);  index  = NULL; }
   if (header) { bam_hdr_destroy(header); header = NULL; }
   if (file)   { hts_close(file);         file   = NULL; }

   refmap.clear();
   filename = "";
}

//------------------------------------------------------------------------------------
// BamRead::BamRead() allocates the data structure used by htslib to represent a read

BamRead::BamRead()
   : aptr(bam_init1())
{
}

//------------------------------------------------------------------------------------
// BamRead::~BamRead() de-allocates the data structure used by htslib to represent a
// read

BamRead::~BamRead()
{
   bam_destroy1(bam1_ptr);
}

//------------------------------------------------------------------------------------
// BamRead::name() returns the name associated with this read; the caller is obligated
// to de-allocate it

std::string *BamRead::name() const
{
   return new std::string(bam_get_qname(bam1_ptr));
}

//------------------------------------------------------------------------------------
// BamRead::constName() returns a const pointer to a C-style string containing the
// name associated with this read; this is more efficient than BamRead::name()
// because the string is not copied; the caller must not attempt to de-allocate this
// C-style string

const char *BamRead::constName() const
{
   return bam_get_qname(bam1_ptr);
}

//------------------------------------------------------------------------------------
// BamRead::bitFlag() returns the bit flag associated with this read; call one of the
// Boolean functions to test a single bit of this flag

int BamRead::bitFlag() const
{
   return bam1_ptr->core.flag;
}

//------------------------------------------------------------------------------------
// BamRead::isPaired() returns true if this read is one mate of a read pair

bool BamRead::isPaired() const
{
   return (bam1_ptr->core.flag & BAM_FPAIRED ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isProperPair() returns true if this read is one mate of a "proper" read
// pair

bool BamRead::isProperPair() const
{
   return (bam1_ptr->core.flag & BAM_FPROPER_PAIR ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isUnmapped() returns true if this read is unmapped

bool BamRead::isUnmapped() const
{
   return (bam1_ptr->core.flag & BAM_FUNMAP ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isMateUnmapped() returns true if this read has an unmapped mate

bool BamRead::isMateUnmapped() const
{
   return (bam1_ptr->core.flag & BAM_FMUNMAP ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isReverseStrand() returns true if this read is mapped to the reverse
// strand

bool BamRead::isReverseStrand() const
{
   return (bam1_ptr->core.flag & BAM_FREVERSE ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isMateReverseStrand() returns true if this read has a mate mapped to the
// reverse strand

bool BamRead::isMateReverseStrand() const
{
   return (bam1_ptr->core.flag & BAM_FMREVERSE ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isRead1() returns true if this read is the first mate of a pair

bool BamRead::isRead1() const
{
   return (bam1_ptr->core.flag & BAM_FREAD1 ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isRead2() returns true if this read is the second mate of a pair

bool BamRead::isRead2() const
{
   return (bam1_ptr->core.flag & BAM_FREAD2 ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isSecondary() returns true if this is a secondary alignment, i.e., it
// is one of multiple mappings of the read sequence and it is not the primary
// alignment

bool BamRead::isSecondary() const
{
   return (bam1_ptr->core.flag & BAM_FSECONDARY ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isFailedQC() returns true if this read failed quality control

bool BamRead::isFailedQC() const
{
   return (bam1_ptr->core.flag & BAM_FQCFAIL ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isDuplicate() returns true if this read has been marked as a duplicate

bool BamRead::isDuplicate() const
{
   return (bam1_ptr->core.flag & BAM_FDUP ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::isSupplementary() returns true if this is a supplementary (chimeric)
// alignment

bool BamRead::isSupplementary() const
{
   return (bam1_ptr->core.flag & BAM_FSUPPLEMENTARY ? true : false);
}

//------------------------------------------------------------------------------------
// BamRead::refID() returns the ID of the reference sequence to which this read has
// been mapped

int BamRead::refID() const
{
   return bam1_ptr->core.tid;
}

//------------------------------------------------------------------------------------
// BamRead::position() returns the 1-based, inclusive leftmost position to which this
// read has been mapped

int BamRead::position() const
{
   return bam1_ptr->core.pos + 1;
}

//------------------------------------------------------------------------------------
// BamRead::endPosition() returns the 1-based, inclusive rightmost position to which
// this read has been mapped

int BamRead::endPosition() const
{
   return bam_endpos(bam1_ptr);
}

//------------------------------------------------------------------------------------
// BamRead::mapQuality() returns the mapping quality score for this read

int BamRead::mapQuality() const
{
   return bam1_ptr->core.qual;
}

//------------------------------------------------------------------------------------
// BamRead::numCigarOps() returns the number of CIGAR operatiosn for this read

int BamRead::numCigarOps() const
{
   return bam1_ptr->core.n_cigar;
}

//------------------------------------------------------------------------------------
// BamRead::cigarOp() returns the operation character (one of MIDNSHP=XB) and length
// for the CIGAR operation at the specified index

void BamRead::cigarOp(int index, char& op, int& len) const
{
   if (index >= 0 && index < numCigarOps())
   {
      const uint32_t *cigar = bam_get_cigar(bam1_ptr);
      op  = bam_cigar_opchr(cigar[index]);
      len = bam_cigar_oplen(cigar[index]);
   }
   else // the index is out of bounds
   {
      op  = '?';
      len = 0;
   }
}

//------------------------------------------------------------------------------------
// BamRead::mateRefID() returns the ID of the reference sequence to which the mate is
// mapped

int BamRead::mateRefID() const
{
   return bam1_ptr->core.mtid;
}

//------------------------------------------------------------------------------------
// BamRead::matePosition() returns the 1-based, inclusive leftmost position to which
// the mate is mapped

int BamRead::matePosition() const
{
   return bam1_ptr->core.mpos + 1;
}

//------------------------------------------------------------------------------------
// BamRead::insertSize() returns the insert size for this read pair

int BamRead::insertSize() const
{
   return bam1_ptr->core.isize;
}

//------------------------------------------------------------------------------------
// BamRead::length() returns the length of the read sequence

int BamRead::length() const
{
   return bam1_ptr->core.l_qseq;
}

//------------------------------------------------------------------------------------
// BamRead::sequence() returns the read sequence; the caller is obligated to
// de-allocate it

std::string *BamRead::sequence() const
{
   const uint8_t *qseq = bam_get_seq(bam1_ptr);

   int readlen = length();
   std::string *seq = new std::string(readlen, ' ');

   for (int i = 0; i < readlen; i++)
      (*seq)[i] = seq_nt16_str[bam_seqi(qseq, i)];

   return seq;
}

//------------------------------------------------------------------------------------
// BamRead::qualities() returns a string of ASCII quality scores; the caller is
// obligated to de-allocate it

std::string *BamRead::qualities() const
{
   const uint8_t *phred = bam_get_qual(bam1_ptr);

   int readlen = length();
   std::string *ascii = new std::string(readlen, ' ');

   for (int i = 0; i < readlen; i++)
      (*ascii)[i] = phred[i] + 33;

   return ascii;
}

//------------------------------------------------------------------------------------
// BamReader::BamReader() allocates data structures used by htslib for reading a Bam
// file

BamReader::BamReader()
   : rptr(new InternalData())
{
}

//------------------------------------------------------------------------------------
// BamReader::~BamReader() de-allocates the data structures used by htslib for reading
// a Bam file

BamReader::~BamReader()
{
   delete data_ptr;
}

//------------------------------------------------------------------------------------
// BamReader::open() opens the Bam file with the specified name

void BamReader::open(const std::string& filename)
{
   if (data_ptr->filename != "")
      throw std::runtime_error("attempt to open " + filename + " when " +
                               data_ptr->filename + " is already open");

   data_ptr->filename = filename;

   if (!(data_ptr->file = hts_open(filename.c_str(), "r")))
      throw std::runtime_error("unable to open " + filename);

   if (!(data_ptr->header = sam_hdr_read(data_ptr->file)))
      throw std::runtime_error("unable to read header in " + filename);

   int numref = data_ptr->header->n_targets;

   for (int refID = 0; refID < numref; refID++)
   {
      std::string refName = data_ptr->header->target_name[refID];

      if (data_ptr->refmap.find(refName) == data_ptr->refmap.end())
         data_ptr->refmap.insert(std::make_pair(refName, refID));
   }
}

//------------------------------------------------------------------------------------
// BamReader::getNumRef() returns the number of reference sequences in the Bam file
// header or zero if there is no file open

int BamReader::getNumRef() const
{
   return (data_ptr->header ? data_ptr->header->n_targets : 0);
}

//------------------------------------------------------------------------------------
// BamReader::getRefID() returns the ID of the reference sequence having the given
// name; (-1) is returned if there is no reference sequence having the given name

int BamReader::getRefID(const std::string& refName)
{
   RefMap::iterator refpos = data_ptr->refmap.find(refName);

   return (refpos == data_ptr->refmap.end() ? -1 : refpos->second);
}

//------------------------------------------------------------------------------------
// BamReader::getRefIDAlt() returns the ID of the reference sequence having the given
// name; (-1) is returned if there is no reference sequence having the given name;
// whereas BamReader::getRefID() requires the given name to match verbatim a
// reference sequence name in the Bam file header, this function also checks to see if
// the given name can be matched by removing or adding a "chr" prefix

int BamReader::getRefIDAlt(const std::string& refName)
{
   int refID = getRefID(refName);
   if (refID >= 0)
      return refID;

   // the given name does not appear verbatim in the Bam file header

   if (refName.length() > 3 && std::tolower(refName[0]) == 'c' &&
       std::tolower(refName[1]) == 'h' && std::tolower(refName[2]) == 'r')
   {
      // try it without the "chr" prefix
      std::string revisedRefName = refName.substr(3);
      return getRefID(revisedRefName);
   }

   // try it with the "chr" prefix
   std::string revisedRefName = "chr" + refName;
   return getRefID(revisedRefName);
}

//------------------------------------------------------------------------------------
// BamReader::getRefLen() returns the length of the identified reference sequence or
// (-1) if there is no file open or the reference ID is out of bounds

int BamReader::getRefLen(int refID) const
{
   if (data_ptr->header && refID >= 0 && refID < data_ptr->header->n_targets)
      return data_ptr->header->target_len[refID];
   else
      return -1;
}

//------------------------------------------------------------------------------------
// BamReader::getRefName() returns the name of the identified reference sequence or
// "UNKNOWN" if there is no file open or the reference ID is out of bounds

std::string BamReader::getRefName(int refID) const
{
   if (data_ptr->header && refID >= 0 && refID < data_ptr->header->n_targets)
      return data_ptr->header->target_name[refID];
   else
      return "UNKNOWN";
}

//------------------------------------------------------------------------------------
// BamReader::jump() does a seek operation to the specified location in the Bam file;
// the start position is assumed to be 1-based

void BamReader::jump(int refID, int startPosition)
{
   if (!data_ptr->file)
      throw std::runtime_error("attempt to jump in Bam file that is not open");

   if (refID < 0 || refID >= data_ptr->header->n_targets || startPosition < 1)
      throw std::runtime_error("invalid argument passed to BamReader::jump()");

   if (!data_ptr->index &&
       (!(data_ptr->index = sam_index_load(data_ptr->file,
                                           data_ptr->filename.c_str()))))
      throw std::runtime_error("unable to read index file for " + data_ptr->filename);

   if (data_ptr->iter)
      hts_itr_destroy(data_ptr->iter);

   if (!(data_ptr->iter = sam_itr_queryi(data_ptr->index, refID,
                                         startPosition - 1, // 0-based position
					 getRefLen(refID))))
      throw std::runtime_error("disk seek failed on reference sequence " +
                               getRefName(refID) + " in " + data_ptr->filename);
}

//------------------------------------------------------------------------------------
// BamReader::getNext() reads the next read from the Bam file; true is returned if
// successful; false is returned when there are no more reads (if a jump was
// performed, this will be when the end of the reference sequence is reached;
// otherwise, it will be when the end of the Bam file is reached)

bool BamReader::getNext(BamRead& read)
{
   if (!data_ptr->file)
      throw std::runtime_error("attempt to read from Bam file that is not open");

   int status;

   if (data_ptr->iter)
      status = sam_itr_next(data_ptr->file, data_ptr->iter, (bam1_t *)read.aptr);
   else
      status = sam_read1(data_ptr->file,  data_ptr->header, (bam1_t *)read.aptr);

   if (status >= 0)
      return true;

   if (status == -1)
      return false;

   // status < -1
   throw std::runtime_error("error reading " + data_ptr->filename);
}

//------------------------------------------------------------------------------------
// BamReader::close() closes the Bam file and reinitializes the data fields

void BamReader::close()
{
   data_ptr->reinitialize();
}
