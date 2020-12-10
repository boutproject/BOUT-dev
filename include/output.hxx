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
#include <functional>
#include <string>

#include "bout/assert.hxx"
#include "boutexception.hxx"
#include "unused.hxx"
#include "bout/format.hxx"
#include "bout/sys/gettext.hxx"  // for gettext _() macro

#include "fmt/core.h"

using std::endl;

/// Class for text output to stdout and/or log file
/*!
  This class can be used to output either in fmt format:

    output.write("A string {:s} and number {:d}\n", str, i);

  or as a C++ stream buffer:

    output << "A string " << str << " and number " << i << endl;

  If a file has been opened (i.e. the processor's log file) then the string
  will be written to the file. In addition, output to stdout can be enabled
  and disabled.
*/
class Output : private multioutbuf_init<char, std::char_traits<char>>,
               public std::basic_ostream<char, std::char_traits<char>> {

  using _Tr = std::char_traits<char>;
  using multioutbuf_init = ::multioutbuf_init<char, _Tr>;

public:
  Output() : multioutbuf_init(), std::basic_ostream<char, _Tr>(multioutbuf_init::buf()) {
    Output::enable();
  }

  /// Specify a log file to open
  Output(const std::string& filename) {
    Output::enable();
    open(filename);
  }

  template <class S, class... Args>
  Output(const S& format, const Args&... args) : Output(fmt::format(format, args...)) {}

  ~Output() override {
    close();
  }

  virtual void enable();  ///< Enables writing to stdout (default)
  virtual void disable(); ///< Disables stdout

  /// Open an output log file
  int open(const std::string& filename);

  template <class S, class... Args>
  int open(const S& format, const Args&... args) {
    return open(fmt::format(format, args...));
  }

  /// Close the log file
  void close();

  /// Write a string using fmt format
  virtual void write(const std::string& message);

  template <class S, class... Args>
  void write(const S& format, const Args&... args) {
    write(fmt::format(format, args...));
  }
  /// Same as write, but only to screen
  virtual void print(const std::string& message);

  template <class S, class... Args>
  void print(const S& format, const Args&... args) {
    print(fmt::format(format, args...));
  }

  /// Add an output stream. All output will be sent to all streams
  void add(std::basic_ostream<char, _Tr> &str) { multioutbuf_init::buf()->add(str); }

  /// Remove an output stream
  void remove(std::basic_ostream<char, _Tr> &str) {
    multioutbuf_init::buf()->remove(str);
  }

  static Output *getInstance(); ///< Return pointer to instance

protected:
  friend class ConditionalOutput;
  virtual Output *getBase() { return this; }
  virtual bool isEnabled() { return true; }

private:
  std::ofstream file;                 ///< Log file stream
  bool enabled;                       ///< Whether output to stdout is enabled
};

/// Class which behaves like Output, but has no effect.
/// This is to allow debug outputs to be disabled at compile time
/// 
/// 
class DummyOutput : public Output {
public:
  template <class S, class... Args>
  void write(const S&, const Args&...) {};
  template <class S, class... Args>
  void print(const S&, const Args&...) {};
  void enable() override{};
  void disable() override{};
  void enable(MAYBE_UNUSED(bool enable)){};
  bool isEnabled() override { return false; }
};

/// Layer on top of Output which passes through calls to write, print etc
/// if it is enabled, but discards messages otherwise.
/// This is used to provide different levels of output
/// (info, prog, warn, error) which can be enabled and disabled at run time.
///
class ConditionalOutput : public Output {
public:
  /// @param[in] base    The Output object which will be written to if enabled
  /// @param[in] enabled Should this be enabled by default?
  ConditionalOutput(Output *base, bool enabled = true) : base(base), enabled(enabled) {};

  /// Constuctor taking ConditionalOutput. This allows several layers of conditions
  /// 
  /// @param[in] base    A ConditionalOutput which will be written to if enabled
  /// 
  ConditionalOutput(ConditionalOutput *base)
      : base(base), enabled(base->enabled) {};

  /// If enabled, writes a string using fmt formatting
  /// by calling base->write
  /// This string is then sent to log file and stdout (on processor 0)
  void write(const std::string& message) override;

  template <class S, class... Args>
  void write(const S& format, const Args&... args) {
    if (enabled) {
      ASSERT1(base != nullptr);
      base->write(fmt::format(format, args...));
    }
  }

  /// If enabled, print a string to stdout using fmt formatting
  /// note: unlike write, this is not also sent to log files
  void print(const std::string& message) override;

  template <class S, class... Args>
  void print(const S& format, const Args&... args) {
    if (enabled) {
      ASSERT1(base != nullptr);
      base->print(fmt::format(format, args...));
    }
  }

  /// Get the lowest-level Output object which is the base of this ConditionalOutput
  Output *getBase() override {
    ASSERT1(base != nullptr);
    return base->getBase();
  };

  /// Set whether this ConditionalOutput is enabled
  /// If set to false (disabled), then all print and write calls do nothing
  void enable(bool enable_) { enabled = enable_; };

  /// Turn on outputs through calls to print and write
  void enable() override { enabled = true; };

  /// Turn off outputs through calls to print and write
  /// This includes log files and stdout
  void disable() override { enabled = false; };

  /// Check if output is enabled
  bool isEnabled() override {
    ASSERT1(base != nullptr);
    return enabled && base->isEnabled();
  };

private:
  /// The lower-level Output to send output to
  Output *base;
  /// Does this instance output anything?
  bool enabled;
};

/// Catch stream outputs to DummyOutput objects. This is so that
/// statements like
///    output_debug << "debug message";
/// compile but have no effect if BOUT_USE_OUTPUT_DEBUG is false
template <typename T> DummyOutput &operator<<(DummyOutput &out, T const &UNUSED(t)) {
  return out;
}

template <typename T> DummyOutput &operator<<(DummyOutput &out, const T *UNUSED(t)) {
  return out;
}

// Function pointer so we can apply unused macro to pf in function below
using stream_manipulator = std::ostream &(*)(std::ostream &);

inline DummyOutput &operator<<(DummyOutput &out, stream_manipulator UNUSED(pf)) {
  return out;
}

inline ConditionalOutput &operator<<(ConditionalOutput &out, stream_manipulator pf) {
  if (out.isEnabled()) {
    *out.getBase() << pf;
  }
  return out;
}

template <typename T> ConditionalOutput &operator<<(ConditionalOutput &out, T const &t) {
  if (out.isEnabled()) {
    *out.getBase() << t;
  }
  return out;
}

template <typename T> ConditionalOutput &operator<<(ConditionalOutput &out, const T *t) {
  if (out.isEnabled()) {
    *out.getBase() << t;
  }
  return out;
}

/// To allow statements like "output.write(...)" or "output << ..."
/// Output for debugging
#ifdef BOUT_USE_OUTPUT_DEBUG
extern ConditionalOutput output_debug;
#else
extern DummyOutput output_debug;
#endif
extern ConditionalOutput output_warn;  ///< warnings
extern ConditionalOutput output_progress;  ///< progress
extern ConditionalOutput output_info;  ///< information 
extern ConditionalOutput output_error; ///< errors
extern ConditionalOutput output_verbose; ///< less interesting messages

/// Generic output, given the same level as output_progress
extern ConditionalOutput output;

#endif // __OUTPUT_H__
