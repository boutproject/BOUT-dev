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

#include <stdarg.h>
#include <stdio.h>

class OptionsReader {
 public:
  static OptionsReader *getInstance();
  static void cleanup() {delete instance; instance = NULL;}
  
  void read(Options *options, const char *file, ...);
  
  void parseCommandLine(Options *options, int argc, char **argv);
  
 private:
  static OptionsReader *instance;
  
};

#endif // __OPTIONSREADER_H__

