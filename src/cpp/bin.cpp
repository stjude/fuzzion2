//------------------------------------------------------------------------------------
//
// bin.cpp - module for binary file I/O
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2020 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "bin.h"
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

//------------------------------------------------------------------------------------
// BinReader::BinReader() allocates the internal read buffer

BinReader::BinReader(int bufferSize)
   : filename(""), fd(BINARY_FILE_NOT_OPEN), size(bufferSize), len(0), index(0),
     swap(false)
{
   buf = new uint8_t[size];
}

//------------------------------------------------------------------------------------
// BinReader::open() opens the named binary file for reading

void BinReader::open(const std::string& binFilename)
{
   if (isOpen())
      close();

   filename = binFilename;

   fd = ::open(filename.c_str(), O_RDONLY);

   if (!isOpen())
      throw std::runtime_error("unable to open " + filename);

   len   = 0;
   index = 0;
}

//------------------------------------------------------------------------------------
// BinReader::seek() seeks to the specified byte offset

void BinReader::seek(uint64_t byteOffset)
{
   if (!isOpen())
      throw std::runtime_error("attempt to seek in unopened binary file");

   if (lseek(fd, byteOffset, SEEK_SET) == -1)
      throw std::runtime_error("seek error in " + filename);

   len   = 0;
   index = 0;
}

//------------------------------------------------------------------------------------
// BinReader::fill() performs a disk read, filling the internal read buffer with the
// bytes read; true is returned if at least one byte was read; false is returned to
// indicate that end-of-file has been reached (no bytes were read)

bool BinReader::fill()
{
   if (!isOpen())
      throw std::runtime_error("attempt to read from unopened binary file");

   ssize_t bytes = read(fd, buf, size);

   if (bytes == -1)
      throw std::runtime_error("error reading from " + filename);

   len   = bytes;
   index = 0;

   return (bytes > 0);
}

//------------------------------------------------------------------------------------
// BinReader::readBuffer() copies the specified number of bytes into the given buffer;
// true is returned if the requested number of bytes were read and copied to the
// buffer; false is returned if end-of-file was reached

bool BinReader::readBuffer(void *buffer, int numBytes)
{
   uint8_t *destination = static_cast<uint8_t *>(buffer);

   while (numBytes > 0)
   {
      if (index >= len && !fill())
         return false;

      int bytesToCopy = std::min(len - index, numBytes);

      std::memcpy(destination, &buf[index], bytesToCopy);

      destination += bytesToCopy;
      index       += bytesToCopy;
      numBytes    -= bytesToCopy;
   }

   return true;
}

//------------------------------------------------------------------------------------
// BinReader::skipBytes() skips the specified number of bytes; true is returned if the
// requested number of bytes were skipped; false is returned if end-of-file was
// reached

bool BinReader::skipBytes(int numBytes)
{
   while (numBytes > 0)
   {
      if (index >= len && !fill())
         return false;

      int bytesToSkip = std::min(len - index, numBytes);

      index    += bytesToSkip;
      numBytes -= bytesToSkip;
   }

   return true;
}

//------------------------------------------------------------------------------------
// BinReader::readString() copies a string, up to and including its null terminating
// byte, but not more than maxlen characters; true is returned if successful; false
// is returned if end-of-file was reached

bool BinReader::readString(char *string, int maxlen)
{
   for (int i = 0; i < maxlen; i++)
   {
      char c;

      if (!readBuffer(&c, 1))
         return false;

      string[i] = c;

      if (c == 0) // found null terminating byte
         break;
   }

   return true;
}

//------------------------------------------------------------------------------------
// BinReader::readUint8() reads a 1-byte unsigned integer; true is returned if
// successful; false is returned if end-of-file was reached

bool BinReader::readUint8(uint8_t& value)
{
   return readBuffer(&value, 1);
}

//------------------------------------------------------------------------------------
// BinReader::readUint16() reads a 2-byte unsigned integer; true is returned if
// successful; false is returned if end-of-file was reached

bool BinReader::readUint16(uint16_t& value)
{
   if (!readBuffer(&value, 2))
      return false;

   if (swap)
      swapBytes(&value, 2);

   return true;
}

//------------------------------------------------------------------------------------
// BinReader::readUint32() reads a 4-byte unsigned integer; true is returned if
// successful; false is returned if end-of-file was reached

bool BinReader::readUint32(uint32_t& value)
{
   if (!readBuffer(&value, 4))
      return false;

   if (swap)
      swapBytes(&value, 4);

   return true;
}

//------------------------------------------------------------------------------------
// BinReader::readUint64() reads an 8-byte unsigned integer; true is returned if
// successful; false is returned if end-of-file was reached

bool BinReader::readUint64(uint64_t& value)
{
   if (!readBuffer(&value, 8))
      return false;

   if (swap)
      swapBytes(&value, 8);

   return true;
}

//------------------------------------------------------------------------------------
// BinReader::close() closes the binary file

void BinReader::close()
{
   if (!isOpen())
      return;

   if (::close(fd) == -1)
      throw std::runtime_error("error closing " + filename);

   filename = "";
   fd       = BINARY_FILE_NOT_OPEN;
   len      = 0;
   index    = 0;
}

//------------------------------------------------------------------------------------
// BinWriter::BinWriter() allocates the internal write buffer

BinWriter::BinWriter(int bufferSize)
   : filename(""), fd(BINARY_FILE_NOT_OPEN), size(bufferSize), index(0), flushed(0)
{
   buf = new uint8_t[size];
}

//------------------------------------------------------------------------------------
// BinWriter::open() opens the named binary file for writing; if newFile is true,
// a new file is created; otherwise, an existing file is opened for writing

void BinWriter::open(const std::string& binFilename, bool newFile)
{
   if (isOpen())
      close();

   filename = binFilename;

   if (newFile)
      fd = ::open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC,
                  S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
   else
      fd = ::open(filename.c_str(), O_WRONLY);

   if (!isOpen())
      throw std::runtime_error("unable to open " + filename);

   index   = 0;
   flushed = 0;
}

//------------------------------------------------------------------------------------
// BinReader::flush() performs a disk write, flushing the internal write buffer

void BinWriter::flush()
{
   if (!isOpen())
      throw std::runtime_error("attempt to write to unopened binary file");

   if (index == 0) // nothing to flush
      return;

   ssize_t bytes = write(fd, buf, index);

   if (bytes != index)
      throw std::runtime_error("error writing to " + filename);

   flushed += index;
   index    = 0;
}

//------------------------------------------------------------------------------------
// BinWriter::writeBuffer() copies the specified number of bytes from the given buffer
// to the internal write buffer, flushing the internal buffer if it fills up

void BinWriter::writeBuffer(const void *buffer, int numBytes)
{
   const uint8_t *source = static_cast<const uint8_t *>(buffer);

   while (numBytes > 0)
   {
      if (index >= size)
         flush();

      int bytesToCopy = std::min(size - index, numBytes);

      std::memcpy(&buf[index], source, bytesToCopy);

      source   += bytesToCopy;
      index    += bytesToCopy;
      numBytes -= bytesToCopy;
   }
}

//------------------------------------------------------------------------------------
// BinWriter::writeString() copies a string, up to and including its null terminating
// byte, to the internal writer buffer, flushing the internal buffer if it fills up

void BinWriter::writeString(const char *string)
{
   writeBuffer(string, std::strlen(string) + 1);
}

//------------------------------------------------------------------------------------
// BinWriter::writeUint8() writes a 1-byte unsigned integer

void BinWriter::writeUint8(uint8_t value)
{
   writeBuffer(&value, 1);
}

//------------------------------------------------------------------------------------
// BinWriter::writeUint16() writes a 2-byte unsigned integer

void BinWriter::writeUint16(uint16_t value)
{
   writeBuffer(&value, 2);
}

//------------------------------------------------------------------------------------
// BinWriter::writeUint32() writes a 4-byte unsigned integer

void BinWriter::writeUint32(uint32_t value)
{
   writeBuffer(&value, 4);
}

//------------------------------------------------------------------------------------
// BinWriter::writeUint64() writes an 8-byte unsigned integer

void BinWriter::writeUint64(uint64_t value)
{
   writeBuffer(&value, 8);
}

//------------------------------------------------------------------------------------
// BinWriter::close() closes the binary file

void BinWriter::close()
{
   if (!isOpen())
      return;

   if (index > 0)
      flush();

   if (::close(fd) == -1)
      throw std::runtime_error("error closing " + filename);

   filename = "";
   fd       = BINARY_FILE_NOT_OPEN;
   index    = 0;
   flushed  = 0;
}

//------------------------------------------------------------------------------------
// swapBytes() swaps the byte ordering of an unsigned integer of the specified length

void swapBytes(void *buffer, int numBytes)
{
   uint8_t *byteArray = static_cast<uint8_t *>(buffer);

   int i = 0;
   int j = numBytes - 1;

   while (i < j)
   {
      uint8_t byte = byteArray[i];
      byteArray[i] = byteArray[j];
      byteArray[j] = byte;

      i++;
      j--;
   }
}
