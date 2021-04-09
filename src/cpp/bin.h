//------------------------------------------------------------------------------------
//
// bin.h - module for binary file I/O
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef BIN_H
#define BIN_H

#include <string>

const int DEFAULT_BINARY_BUFFER_SIZE = 16 * 1024 * 1024; // 16 MB
const int BINARY_FILE_NOT_OPEN = -1;

//------------------------------------------------------------------------------------

class BinReader // for reading a binary file
{
public:
   BinReader(int bufferSize=DEFAULT_BINARY_BUFFER_SIZE);

   virtual ~BinReader() { delete[] buf; }

   void open(const std::string& binFilename);
   bool isOpen() const { return (fd != BINARY_FILE_NOT_OPEN); }

   void seek(uint64_t byteOffset);

   bool readBuffer(void *buffer, int numBytes);
   bool skipBytes(int numBytes);

   bool readString(char *string, int maxlen);

   bool readUint8 (uint8_t&  value);
   bool readUint16(uint16_t& value);
   bool readUint32(uint32_t& value);
   bool readUint64(uint64_t& value);

   void close();

   std::string  filename; // name of binary input file
   int          fd;       // file descriptor or (-1) if no file is open
   uint8_t     *buf;      // internal read buffer
   int          size;     // number of bytes in buf
   int          len;      // number of bytes used in buf
   int          index;    // index into buf of next byte to read
   bool         swap;     // set to true to swap byte ordering of multi-byte integers

protected:
   bool fill();
};

//------------------------------------------------------------------------------------

class BinWriter // for writing a binary file
{
public:
   BinWriter(int bufferSize=DEFAULT_BINARY_BUFFER_SIZE);

   virtual ~BinWriter() { delete[] buf; }

   void open(const std::string& binFilename, bool newFile=true);
   bool isOpen() const { return (fd != BINARY_FILE_NOT_OPEN); }

   void writeBuffer(const void *buffer, int numBytes);

   void writeString(const char *string);

   void writeUint8 (uint8_t  value);
   void writeUint16(uint16_t value);
   void writeUint32(uint32_t value);
   void writeUint64(uint64_t value);

   uint64_t bytesWritten() const { return flushed + index; }

   void close();

   std::string  filename; // name of binary output file
   int          fd;       // file descriptor or (-1) if no file is open
   uint8_t     *buf;      // internal write buffer
   int          size;     // number of bytes in buf
   int          index;    // index into buf of next byte to write
   uint64_t     flushed;  // number of bytes flushed

protected:
   void flush();
};

//------------------------------------------------------------------------------------

void swapBytes(void *buffer, int numBytes);

//------------------------------------------------------------------------------------
#endif
