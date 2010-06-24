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

/// Class for reading INI style configuration files
/*!
 * 
 */
class OptionFile {
 public:
  OptionFile();
  OptionFile(const char *filename);
  ~OptionFile();
  
  /// Read options from grid file
  int read(const char *filename, ...);
  
  int command_line(int argc, char** argv);

  int getInt(const char *name, int &val);
  int getInt(const char *section, const char *name, int &val);

  int getReal(const char *name, BoutReal &val);
  int getReal(const char *section, const char *name, BoutReal &val);

  char* getString(const char *name);
  char* getString(const char *section, const char *name);

  int getBool(const char *name, bool &val);
  int getBool(const char *section, const char *name, bool &val);

  // New interface
  // Prints out what values are being assigned

  void setSection(const char *name); // Set the default section
 
  int get(const char *name, int &val, const int def);
  int get(const char *name, BoutReal &val, const BoutReal def);
  int get(const char *name, bool &val, const bool def);

  int get(const char *section, const char *name, int &val, const int def);
  int get(const char *section, const char *name, BoutReal &val, const BoutReal def);
  int get(const char *section, const char *name, bool &val, const bool def);

  // Set methods to pass in options manually
  int set(const char *name, int val);
  int set(const char *name, BoutReal val);
  int set(const char *name, bool val);
  int set(const char *name, const char *string);

 private:
  
  static const char COMMENT_CHAR = ';';
  static const int MAX_LINELEN = 512;

  typedef struct {
    char *name;
    int hash;
    char *string;
  }t_option;

  int noptions;
  t_option *option;
  
  void add(const char *section, const char *name, const char *string, int linenr);
  int find(const char *name);
  unsigned int hash_string(const char *string);
  int strip_space(char *string);
  int get_nextline(FILE *fp, char *buffer, int maxbuffer, int first);
  
  char *def_section; // Default section heading
  
};

#endif // __OPTIONS_H__
