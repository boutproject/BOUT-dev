/*!************************************************************************
* Reads in the configuration file, supplying
* an interface to get options
* 
* File is an ini file with sections
* [section]
* and variables as
* name = string ; comment
* 
* Ben Dudson, September 2007
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

class OptionFile;

#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include "bout_types.h"

#include <stdarg.h>
#include <stdio.h>
#include <map>

/// Class for reading INI style configuration files
/*!
* 
*/

using namespace std;

class OptionFile {
public:
  OptionFile();
  OptionFile(const string &filename);
  OptionFile(int &argc, char **argv, const string &filename);
  ~OptionFile();

  /// Read options from grid file
  void read(const char *filename, ...);

  void commandLineRead(int argc, char** argv);

  // New interface
  // Prints out what values are being assigned

  void setSection(const string &name); // Set the default section
  string getSection(); // Set the default section

  template <class type> void get(const map<string,string>::iterator &, type &);
  template <class type> void get(const string &, type &, const type &);
  template <class type> void get(const string &, const string &, type &, const type &);
  template <class type> void get(const string &, const string &, const string &, type &, const type &);

  void get(const string &, int &, const int &);
  void get(const string &, BoutReal &, const BoutReal &);
  void get(const string &, bool &, const bool &);
  void get(const string &, string &, const string &);

  void get(const string &, const string &, int &, const int &);
  void get(const string &, const string &, BoutReal &, const BoutReal &);
  void get(const string &, const string &, bool &, const bool &);
  void get(const string &, const string &, string &, const string &);

  void get(const string &, const string &, const string &, int &, const int &);
  void get(const string &, const string &, const string &, BoutReal &, const BoutReal &);
  void get(const string &, const string &, const string &, bool &, const bool &);
  void get(const string &, const string &, const string &, string &, const string &);
  
  void get(const string &, const string &, const string &, const string &, int &, const int &);
  void get(const string &, const string &, const string &, const string &, BoutReal &, const BoutReal &);
  void get(const string &, const string &, const string &, const string &, bool &, const bool &);
  void get(const string &, const string &, const string &, const string &, string &, const string &);

  // Set methods to pass in options manually
  template <class type> void set(const string &, const type &);
  
  void set(const string &, const int &);
  void set(const string &, const BoutReal &);
  void set(const string &, const bool &);
  void set(const string &, const string &);
  
  const string& getSectionSep();
  void setSectionSep(const string &);
  
  map<string,string>::iterator find(const string &);
  map<string, string>::iterator find(const string &, const string &);
  
  map<string,string>::iterator end();

protected:

  void add(const string &, const string &, const string &);
  void trim(string &, const string &c=" \t");
  void trimLeft(string &, const string &c=" \t");
  void trimRight(string &, const string &c=" \t");
  void trimComments(string &);
  void parse(const string &, string &, string &);
  string getNextLine(ifstream &);
  string prependSection(const string &section, const string& key);

  string def_section; // Default section heading

  string sep;

  map<string, string> options;
};

#endif // __OPTIONS_H__
