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
#include <msg_stack.hxx>
#include "options_ini.hxx"

#include <algorithm>

using namespace std;

/**************************************************************************
 * Read input file
 **************************************************************************/

void OptionINI::read(Options *options, const string &filename) {
  ifstream fin;
  fin.open(filename.c_str());

  if (!fin.good()) {
    throw BoutException(_("\tOptions file '{:s}' not found\n"), filename);
  }
  
  Options *section = options; // Current section
  do {
    string buffer = getNextLine(fin);
    
    if(!buffer.empty()) {

      if (buffer[0] == '[') {
        // A section header
        
        auto endpos   = buffer.find_last_of(']');
        if (endpos == string::npos) {
          throw BoutException("\t'{:s}': Missing ']'\n\tLine: {:s}", filename, buffer);
        }

        buffer = trim(buffer, "[]");

        if(buffer.empty()) {
          throw BoutException("\t'{:s}': Missing section name\n\tLine: {:s}", filename,
                              buffer);
        }
        
        section = options;
        size_t scorepos;
        while((scorepos = buffer.find_first_of(':')) != string::npos) {
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

        // Ensure that brackets '()' and '[]'are balanced

        // Count net number of opening and closing brackets
        auto count_brackets = [](const string& input) {
                                int nbracket = 0;
                                for( auto ch : input) {
                                  if ((ch == '(') || (ch == '[')) {
                                    ++nbracket;
                                  }
                                  if ((ch == ')') || (ch == ']')) {
                                    --nbracket;
                                  }
                                }
                                return nbracket;
                              };

        // Starting count
        int count = count_brackets(value);

        string firstline = value; // Store the first line for error message
        
        while (count % 2 == 1) {
          // An odd number, so read another line

          if (fin.eof()) {
            throw BoutException("\t'%s': Unbalanced brackets\n\tStarting line: %s", filename.c_str(), firstline.c_str());
          }
          
          string newline = getNextLine(fin);
          count += count_brackets(newline);
          value += newline;
        }
        // Add this to the current section
        section->set(key, value, filename);
      } // section test
    } // buffer.empty
  } while(!fin.eof());

  fin.close();
}

void OptionINI::write(Options *options, const std::string &filename) {
  TRACE("OptionsINI::write");
  
  std::ofstream fout;
  fout.open(filename, ios::out | ios::trunc);

  if (!fout.good()) {
    throw BoutException(_("Could not open output file '{:s}'\n"), filename);
  }
  
  // Call recursive function to write to file
  writeSection(options, fout);
  
  fout.close();
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

void OptionINI::parse(const string &buffer, string &key, string &value) {
   // A key/value pair, separated by a '='

  size_t startpos = buffer.find_first_of('=');

  if (startpos == string::npos) {
    // Just set a flag to true
    // e.g. "restart" or "append" on command line
    key = buffer;
    value = string("TRUE");
    return;
  }

  key = trim(buffer.substr(0, startpos), " \t\r\n\"");
  value = trim(buffer.substr(startpos+1), " \t\r\n\"");
  
  if (key.empty()) {
    throw BoutException(_("\tEmpty key\n\tLine: {:s}"), buffer);
  }

  if (key.find(':') != std::string::npos) {
    throw BoutException(_("\tKey must not contain ':' character\n\tLine: {:s}"), buffer);
  }
}

void OptionINI::writeSection(const Options *options, std::ofstream &fout) {
  string section_name = options->str();

  if (section_name.length() > 0) {
    // Print the section name at the start
    fout << "[" << section_name << "]" << endl;
  }
  // Iterate over all values
  for(const auto& it : options->getChildren()) {
    if (it.second.isValue()) {
      auto value = bout::utils::variantToString(it.second.value);
      fout << it.first << " = " << value;

      if (value.empty()) {
        // Print an empty string as ""
        fout << "\"\""; 
      }
      bool in_comment = false; // Has a '#' been printed yet?
      
      if (! it.second.valueUsed() ) {
        fout << "\t\t# not used ";
        in_comment = true;
          
        if (it.second.attributes.count("source")) {
          fout << ", from: "
               << it.second.attributes.at("source").as<std::string>();
        }
      }

      if (it.second.attributes.count("type")) {
        if (!in_comment) {
          fout << "\t\t# type: ";
          in_comment = true;
        } else {
          fout << ", type: ";
        }
        fout << it.second.attributes.at("type").as<std::string>();
      }
      
      if (it.second.attributes.count("doc")) {
        if (!in_comment) {
          fout << "\t\t# ";
          in_comment = true;
        } else {
          fout << ", doc: ";
        }
        fout << it.second.attributes.at("doc").as<std::string>();
      }
      fout << endl;
    }
  }
  
  // Iterate over sub-sections
  for(const auto& it : options->subsections()) {
    fout << endl;
    writeSection(it.second, fout);
  }
}
