/*!************************************************************************
* Singleton class for reading options files
*
* Uses a bridge pattern to access OptionParser classes to parse
* different file formats
*
* Handles the command-line parsing
* 
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

class OptionsReader;

#ifndef __OPTIONSREADER_H__
#define __OPTIONSREADER_H__

#include "options.hxx"
#include "bout/format.hxx"

#include "fmt/format.h"

#include <string>

/// Class to handle reading options from file
///
/// Example
/// -------
///
/// Options opt;
/// OptionsReader::getInstance()->read(&opt, "somefile.inp");
///
/// opt now contains a tree of sections and options from the input file "somefile.inp"
///
class OptionsReader {
 public:
  /// Return a pointer to the instance singleton
  static OptionsReader *getInstance();

  /// Delete the instance
  static void cleanup() {delete instance;
    instance = nullptr;
  }

  /// Read the given file, parse options into
  /// the options tree.
  ///
  /// @param[inout] options  The options section to insert values and subsections into
  /// @param[in] file  The name of the file. printf style arguments can be used to create the file name.
  void read(Options *options, const std::string& filename);

  template <class S, class... Args>
  void read(Options* options, const S& format, const Args&... args) {
    return read(options, fmt::format(format, args...));
  }

  /// Write options to file
  ///
  /// @param[in] options  The options tree to be written
  /// @param[in] file   The name of the file to (over)write
  void write(Options *options, const std::string& filename);

  template <class S, class... Args>
  void write(Options* options, const S& format, const Args&... args) {
    return write(options, fmt::format(format, args...));
  }

  /// Parse options from the command line
  ///
  /// @param[inout] options The options section to insert values and subsections into
  /// @param[in] argc   The number of command-line arguments
  /// @param[in] argv   The command line arguments
  ///
  /// Example
  /// -------
  ///
  /// int main(int argc, char** argv) {
  ///   Options opt;
  ///   OptionsReader::getInstance()->read(&opt, argc, argv);
  ///   ...
  ///   return 0;
  /// }
  void parseCommandLine(Options *options, int argc, char **argv);
  
 private:
  /// The instance of this singleton
  static OptionsReader *instance;
  
};

#endif // __OPTIONSREADER_H__

