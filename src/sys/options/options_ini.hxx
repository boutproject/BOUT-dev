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

class OptionINI;

#ifndef __OPTIONS_INI_H__
#define __OPTIONS_INI_H__

#include "optionparser.hxx"

#include <string>
#include <fstream>
using namespace std;

/// Class for reading INI style configuration files
/*!
 *
 */
class OptionINI : public OptionParser {
public:
  OptionINI();
  ~OptionINI();

  /// Read options from grid file
  void read(Options *options, const string &filename);

private:

  void parse(const string &, string &, string &);
  string getNextLine(ifstream &);
};

#endif // __OPTIONS_INI_H__
