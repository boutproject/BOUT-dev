/**************************************************************************
 * Reads in the configuration file, supplying
 * an interface to get options
 *
 * File is an ini file with sections
 * [section]
 * and variables as
 * name = string # comment
 *
 * [section:subsection]  # Sub-sections separated by colons. Arbitrary depth.
 *
 * ChangeLog
 * =========
 *
 * 2011-02-12 Ben Dudson <bd512@york.ac.uk>
 *    * Rearranged to implement OptionParser interface
 *
 * 2010-06 Sean Farley
 *    * Re-written to use STL maps
 *
 * 2010-02-10 Ben Dudson <bd512@york.ac.uk>
 *
 *    * Adding set methods to allow other means to control code
 *      Intended to help with integration into FACETS
 *
 * 2007-09-01 Ben Dudson <bd512@york.ac.uk>
 *
 *    * Initial version
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

#include <utils.hxx>
#include <boutexception.hxx>
#include "options_ini.hxx"

OptionINI::OptionINI() {
}

OptionINI::~OptionINI() {
}

/**************************************************************************
 * Read input file
 **************************************************************************/

void OptionINI::read(Options *options, const string &filename) {
  ifstream fin;
  fin.open(filename.c_str());

  if (!fin.good()) {
    throw BoutException("\tOptions file '%s' not found\n", filename.c_str());
  }

  Options *section = options; // Current section
  do {
    string buffer = getNextLine(fin);

    if(!buffer.empty()) {

      // Check for section
      size_t startpos, endpos;
      startpos = buffer.find_first_of("[");
      endpos   = buffer.find_last_of("]");

      if (startpos != string::npos) {
        // A section header
        if (endpos == string::npos) throw BoutException("\t'%s': Missing ']'\n\tLine: %s", filename.c_str(), buffer.c_str());

        buffer = trim(buffer, "[]");

        if(buffer.empty()) throw BoutException("\t'%s': Missing section name\n\tLine: %s", filename.c_str(), buffer.c_str());

        section = options;
        size_t scorepos;
        while((scorepos = buffer.find_first_of(":")) != string::npos) {
          // sub-section
          string sectionname = trim(buffer.substr(0,scorepos));
          buffer = trim(buffer.substr(scorepos+1));

          section = section->getSection(sectionname);
        }
        section = section->getSection(buffer);
      } else {
        // A key=value pair

        string key, value;
        // Get a key = value pair
        parse(buffer, key, value);
        // Add this to the current section
        section->set(key, value, filename);
      } // section test
    } // buffer.empty
  } while(!fin.eof());

  fin.close();
}

/**************************************************************************
 * Private functions
 **************************************************************************/

// Returns the next useful line, stripped of comments and whitespace
string OptionINI::getNextLine(ifstream &fin) {
  string line;

  getline(fin, line);
  line = lowercasequote(trim(trimComments(line))); // lowercase except for inside quotes

  return line;
}

void OptionINI::parse(const string &buffer, string &key, string &value)
{
  // A key/value pair, separated by a '='

  size_t startpos = buffer.find_first_of("=");
  size_t endpos   = buffer.find_last_of("=");

  if (startpos == string::npos) {
    // Just set a flag to true
    // e.g. "restart" or "append" on command line
    key = buffer;
    value = string("TRUE");
    return;
  }

  key = trim(buffer.substr(0, startpos), " \t\r\n\"");
  value = trim(buffer.substr(startpos+1), " \t\r\n\"");

  if(key.empty() || value.empty()) throw BoutException("\tEmpty key or value\n\tLine: %s", buffer.c_str());
}
