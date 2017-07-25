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

#include "multiostream.hxx"
#include <iostream>
#include <fstream>
#include <boutexception.hxx>
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
    buffer_len=BUFFER_LEN;
    buffer=new char[buffer_len];
    enable();
  }
    
  /// Specify a log file to open
  Output(const char *fname) : multioutbuf_init(), 
    std::basic_ostream<char, _Tr>(multioutbuf_init::buf()) {
    buffer_len=BUFFER_LEN;
    buffer=new char[buffer_len];
    enable();
    open(fname);
  } 
  virtual ~Output() {close();
    delete[] buffer;}
  
  virtual void enable();  ///< Enables writing to stdout (default)
  virtual void disable(); ///< Disables stdout

  int open(const char *fname, ...); ///< Open an output log file
  void close();                    ///< Close the log file

  virtual void write(const char*string, ...); ///< Write a string using C printf format

  virtual void print(const char*string, ...); ///< Same as write, but only to screen

  virtual void vwrite(const char*string, va_list args); ///< Write a string using C vprintf format

  virtual void vprint(const char*string, va_list args); ///< Same as vwrite, but only to screen

  /// Add an output stream. All output will be sent to all streams
  void add(std::basic_ostream<char, _Tr>& str) {
    multioutbuf_init::buf()->add(str);
  }

  /// Remove an output stream
  void remove(std::basic_ostream<char, _Tr>& str) {
    multioutbuf_init::buf()->remove(str);
  }

  static Output *getInstance(); ///< Return pointer to instance
  static void cleanup();   ///< Delete the instance

 private:
  static Output *instance; ///< Default instance of this class
  
  std::ofstream file; ///< Log file stream
  static const int BUFFER_LEN=1024; ///< default length
  int buffer_len; ///< the current length
  char * buffer; ///< Buffer used for C style output
  bool enabled;      ///< Whether output to stdout is enabled
};


class DummyOutput: public Output{
public:
  void write(const char * str,...)override{};
  void print(const char * str,...)override{};
  void enable()override{throw BoutException("DummyOutput cannot be enabled.\nTry compiling with --enable-debug or be less verbose?");};
  void disable()override{};
  void enable(bool en){if (en) this->enable();};
};

class ConditionalOutput: public Output{
public:
  ConditionalOutput(Output * base_):
    base(base_), enabled(true), base_is_cond(false){};
  ConditionalOutput(ConditionalOutput * base_):
    base(base_), enabled(base_->enabled), base_is_cond(true){};
  void write(const char * str,...)override;
  void vwrite(const char * str, va_list va)override{
    if (enabled){
      base->vwrite(str,va);
    }
  }
  void print(const char * str,...)override;
  void vprint(const char * str, va_list va)override{
    if (enabled){
      base->vprint(str,va);
    }
  }
  Output * getBase(){
    if (base_is_cond){
      return dynamic_cast<ConditionalOutput*>(base)->getBase();
    } else {
      return base;
    }
  };
  void enable(bool enable_){enabled=enable_;};
  void enable()override{enabled=true;};
  void disable()override{enabled=false;};
  bool isEnabled(){return enabled && (!base_is_cond || (dynamic_cast<ConditionalOutput*>(base))->isEnabled());};

  Output * base;
  bool enabled;
private:
  bool base_is_cond;
private:
};


template<typename T>
DummyOutput & operator <<(DummyOutput& out, T const & t)
{
  return out;
}

template<typename T>
DummyOutput & operator <<(DummyOutput& out, const T * t)
{
  return out;
}

inline DummyOutput &
operator <<(DummyOutput& out, std::ostream & (*pf)(std::ostream &))
{
  return out;
}

inline ConditionalOutput &
operator <<( ConditionalOutput& out, std::ostream & (*pf)(std::ostream &))
{
  if (out.isEnabled()) {
    *out.getBase() << pf;
  }
  return out;
};

template<typename T>
ConditionalOutput & operator <<(ConditionalOutput& out, T const & t) {
  if (out.isEnabled()) {
    *out.getBase() << t;
  }
  return out;
};

template<typename T>
ConditionalOutput & operator <<(ConditionalOutput& out, const T * t)
{
  if (out.isEnabled()) {
    *out.getBase() << t;
  }
  return out;
};



/// To allow statements like "output.write(...)" or "output << ..."
/// Output for debugging
#ifdef DEBUG_ENABLED
extern Output output_debug;
#else
extern DummyOutput output_debug;
#endif
extern ConditionalOutput output_warn;
extern ConditionalOutput output_prog;
extern ConditionalOutput output_info;
extern ConditionalOutput output_error;


#endif // __OUTPUT_H__

// Allow to reinclude this file, to set output macro again ...
#ifndef output
/// the old output should not be used anymore
#define output (_Pragma("GCC warning \"DEPRECATED: use debug, info, warn or error instead of output.\"") *Output::getInstance())
/// disable all old output
//#define output (*dynamic_cast<DummyOutput*>(Output::getInstance()))
/// use the old output without warning
//#define output (*Output::getInstance())
#endif
