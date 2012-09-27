/**************************************************************************
 * Output, for printing messages/errors etc.
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

class Output;

#pragma once
#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <stdio.h>
#include "multiostream.hxx"
#include <iostream>
#include <fstream>

using std::endl;

/// Class for text output to stdout and/or log file
/*!
  This class can be used to output either in C printf format:

    output.write("A string %s and number %d\n", str, i);

  or as a C++ stream buffer:

    output << "A string " << str << " and number " << i << endl;

  If a file has been opened (i.e. the processor's log file) then the string
  will be written to the file. In addition, output to stdout can be enabled
  and disabled.
*/
class Output : private multioutbuf_init<char, std::char_traits<char> >, 
               public std::basic_ostream<char, std::char_traits<char> > {

  typedef std::char_traits<char> _Tr;
  typedef ::multioutbuf_init<char, _Tr> multioutbuf_init;
  
 public:
  Output() : multioutbuf_init(), 
    std::basic_ostream<char, _Tr>(multioutbuf_init::buf()) {
    
    enable();
  }
    
  /// Specify a log file to open
  Output(const char *fname) : multioutbuf_init(), 
    std::basic_ostream<char, _Tr>(multioutbuf_init::buf()) {
    
    enable();
    open(fname);
  } 
  ~Output() {close();}
  
  void enable();  ///< Enables writing to stdout (default)
  void disable(); ///< Disables stdout

  int open(const char *fname, ...); ///< Open an output log file
  void close();                    ///< Close the log file

  void write(const char*string, ...); ///< Write a string using C printf format

  void print(const char*string, ...); ///< Same as write, but only to screen
  
  void add(std::basic_ostream<char, _Tr>& str) {
    multioutbuf_init::buf()->add(str);
  }

  void remove(std::basic_ostream<char, _Tr>& str) {
    multioutbuf_init::buf()->remove(str);
  }

  static Output *getInstance(); ///< Return pointer to instance
  static void cleanup();
 private:
  static Output *instance; ///< Default instance of this class
  
  std::ofstream file; ///< Log file stream
  char buffer[1024]; ///< Buffer used for C style output
  bool enabled;      ///< Whether output to stdout is enabled
};

/// To allow statements like "output.write(...)" or "output << ..."
#define output (*Output::getInstance())

#endif // __OUTPUT_H__
