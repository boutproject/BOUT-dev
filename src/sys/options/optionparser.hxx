/*!************************************************************************
 * Base class for options file parsers
 *
 * This code handles most details of looking up variables and just relies on
 * the underlying library to supply a simplified interface
 *
 * Original BOUT++ inp file has the form of a Windows INI file
 * with a format
 *  [section]
 *  name = value  # comments
 *
 * This is compatible with some readers (like Python's ConfigParser), but
 * is quite limited.
 *
 * To handle more complex data types and make interchange with other
 * codes easier JSON formatted files are planned to be supported
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

class OptionParser;

#ifndef __OPTIONPARSER_H__
#define __OPTIONPARSER_H__

#include "bout_types.hxx"
#include "options.hxx"

using namespace std;

class OptionParser {
public:
  OptionParser() {}
  virtual ~OptionParser() {}

  virtual void read(Options *options, const string &filename) = 0;
private:
};

#endif // __OPTIONPARSER_H__
