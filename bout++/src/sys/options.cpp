/**************************************************************************
 * Reads in the configuration file, supplying
 * an interface to get options
 * 
 * File is an ini file with sections
 * [section]
 * and variables as
 * name = string ; comment
 * 
 * To Do / Known issues
 * ====================
 *
 *  * Currently uses its own hash function. Probably better to just
 *    use an STL map for this class
 *
 * ChangeLog
 * =========
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

#include "globals.h"
#include "options.h"
#include "utils.h"
#include "boutexception.h"

#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <stdlib.h>
#include <typeinfo>
#include <sstream>

OptionFile::OptionFile() : sep("_")
{
}

OptionFile::OptionFile(const string &filename) : sep("_")
{
  read(filename.c_str());
}

OptionFile::OptionFile(int &argc, char **argv, const string &filename) : sep("_")
{
  read(filename.c_str());
  commandLineRead(argc, argv);
}

OptionFile::~OptionFile()
{
}


/**************************************************************************
 * Read input file
 **************************************************************************/

void OptionFile::read(const char *format, ...)
{
  ifstream fin;

  string buffer;
  string section; // Current section
  
  size_t startpos, endpos;

  va_list ap;  // List of arguments
  char filename[512];
  
  if(format == (const char*) NULL) {
    throw BoutException("ERROR: OptionFile::read passed NULL filename\n");
  }

  va_start(ap, format);
  vsprintf(filename, format, ap);
  va_end(ap);

  output.write("Reading options file %s\n", filename);
  
  fin.open(filename);

  if(!fin.good()) {
    throw BoutException("\tOptions file '%s' not found\n", filename);
  }

  do {
    buffer = getNextLine(fin);
    
    if(!buffer.empty()) {

      // Check for section
      startpos = buffer.find_first_of("[");
      endpos   = buffer.find_last_of("]");

      if( startpos != string::npos ) {
        if( endpos == string::npos ) {
          throw BoutException("\t'%s': Missing ']'\n\tLine: %s", filename, buffer.c_str());
        }

        trim(buffer, "[]");

        if(!buffer.empty()) {
          section = buffer;
        }
      } else {
        
        string key, value;
        
        parse(buffer, key, value);

        add(section, key, value);
      } // section test
    } // buffer.empty
  } while(!fin.eof());

  fin.close();
  
  if(options.empty()) {
    throw BoutException("\tEmpty option file '%s'\n", filename);
  }
  
/*  for (map<string,string>::iterator it=options.begin() ; it != options.end(); it++ )
      cout << (*it).first << " => " << (*it).second << endl;*/
}

void OptionFile::commandLineRead(int argc, char** argv)
{
  output << "Checking command-line options\n";

  string buffer;
  string key, value;

  // Go through command-line arguments
  for(size_t i=1;i<argc;i++) {
    // Should contain a "key=value" string
    buffer = argv[i];

    if(buffer != "restart") {
      parse(buffer, key, value);
      add("", key, value);
    }

    /// No equality sign. Assume it's a setting which is set to true
    /// This is so that adding "restart" to the command-line works
    //  add(NULL, argv[i], "true", -1);
  }

}

/**************************************************************************
 * Functions to request options
 **************************************************************************/

void OptionFile::setSection(const string &name) // Set the default section
{
  def_section = name;
}

string OptionFile::getSection() // Set the default section
{
  return def_section;
}

void OptionFile::setSectionSep(const string &s)
{
  sep = s;
}

inline const string& OptionFile::getSectionSep()
{
  return sep; ///< Used to separate sections and keys
}

inline string OptionFile::prependSection(const string &section, const string& key)
{
  if(!section.empty())
    return lowercase(section + getSectionSep() + key);
    
  return lowercase(key);
}

template <class type>
void OptionFile::get(const map<string,string>::iterator &it, type &val)
{
  if(it != end()) {
    
    stringstream ss;

    if(typeid(type) == typeid(bool)) { // Best way (that I can find) to convert strings to bool
      char c = toupper((it->second)[0]);
      if((c == 'Y') || (c == 'T') || (c == '1')) {
        ss << "1";
        ss >> val;
      } else if((c == 'N') || (c == 'F') || (c == '0')) {
        ss << "0";
      } else {  
          output << "\tOption '" << it->first << "': Boolean expected\n";
      }

    } else {
      ss << it->second;
      output << "\tOption " << it->first << " = " << val << endl;
    }
    
    ss >> val;

  }
}

template<class type>
void OptionFile::get(const string &key, type &val, const type &def)
{
  map<string, string>::iterator it(find(key));

  if(it != end()) {
    get<type>(it, val);
    return;
  }
  
  it = find(prependSection(def_section, key));
  if(it != end()) {
    get<type>(it, val);
    return;
  }

  val = def;
  output << "\tOption " << key << " = " << def << " (default)" << endl;
}

template<class type>
void OptionFile::get(const string &section, const string &key, type &val, const type &def)
{
  if(key.empty()) {
    output.write("WARNING: NULL option requested\n");
    return;
  }
  
  get<type>(key, val, def);
  
  if(val != def) return;

  get<type>(prependSection(section, key), val, def);
}

template<class type>
void OptionFile::get(const string &section1, const string &section2, const string &key, type &val, const type &def)
{
  
  get<type>(key, val, def);
  
  if(val != def) return;
  
  get<type>(prependSection(section1, key), val, def);

  if(val != def) return;
  
  get<type>(prependSection(section2, key), val, def);
}

void OptionFile::get(const string &key, int &val, const int &def)
{
  get<int>(key, val, def);
}

void OptionFile::get(const string &key, BoutReal &val, const BoutReal &def)
{
  get<BoutReal>(key, val, def);
}

void OptionFile::get(const string &key, bool &val, const bool &def)
{
  get<bool>(key, val, def);
}

void OptionFile::get(const string &section, const string &key, int &val, const int &def)
{
  get<int>(section, key, val, def);
}

void OptionFile::get(const string &section1, const string &section2, const string &key, int &val, const int &def)
{
  get<int>(section1, section2, key, val, def);
}

void OptionFile::get(const string &section, const string &key, BoutReal &val, const BoutReal &def)
{
  get<BoutReal>(section, key, val, def);
}

void OptionFile::get(const string &section1, const string &section2, const string &key, BoutReal &val, const BoutReal &def)
{
  get<BoutReal>(section1, section2, key, val, def);
}

void OptionFile::get(const string &section, const string &key, bool &val, const bool &def)
{
  get<bool>(section, key, val, def);
}

void OptionFile::get(const string &key, string &val, const string &def)
{
  get<string>(key, val, def);
}

void OptionFile::get(const string &section, const string &key, string &val, const string &def)
{
  get<string>(section, key, val, def);
}

/**************************************************************************
 * Set methods
 *
 * Convert given options to strings
 **************************************************************************/

template<class type>
void OptionFile::set(const string &key, const type &val)
{
  stringstream ss;
  
  ss << val;
  
  options[key] = ss.str();
}

void OptionFile::set(const string &key, const int &val)
{
  set<int>(key, val);
}

void OptionFile::set(const string &key, const BoutReal &val)
{
  set<BoutReal>(key, val);
}

void OptionFile::set(const string &key, const bool &val)
{
  set<bool>(key, val);
}

void OptionFile::set(const string &key, const string &val)
{
  set<string>(key, val);
}

/**************************************************************************
 * Private functions
 **************************************************************************/

void OptionFile::add(const string &section, const string &key, const string &value)
{
  if(key.empty()) {
    throw BoutException("\tEmpty key passed to 'OptionsFile::add!'");
  }

  string sectionkey(prependSection(section, key));

  options[sectionkey] = value;
}

map<string,string>::iterator OptionFile::find(const string &key)
{
  return options.find(key);
}

map<string,string>::iterator OptionFile::find(const string &section, const string &key)
{
  return options.find(prependSection(section, key));
}

map<string,string>::iterator OptionFile::end()
{
  return options.end();
}

// Strips leading and trailing spaces from a string
void OptionFile::trim(string &s, const string &c)
{
  
  // Find the first character position after excluding leading blank spaces
  size_t startpos = s.find_first_not_of(c);
  // Find the first character position from reverse af
  size_t endpos = s.find_last_not_of(c);

  // if all spaces or empty, then return an empty string
  if(( startpos == string::npos ) || ( endpos == string::npos ))
  {
    s = "";
  }
  else
  {
    s = s.substr(startpos, endpos-startpos+1);
  }
}

void OptionFile::trimRight(string &s, const string &c)
{  
  // Find the first character position after excluding leading blank spaces
  size_t startpos = s.find_first_of(c);

  // if all spaces or empty, then return an empty string
  if(( startpos != string::npos ) )
  {
    s = s.substr(0,startpos);
  }
}

void OptionFile::trimLeft(string &s, const string &c)
{
  
  // Find the first character position from reverse af
  size_t endpos = s.find_last_of(c);

  // if all spaces or empty, then return an empty string
  if(( endpos != string::npos ) )
  {  
    s = s.substr(endpos);
  }
}

// Strips the comments from a string
void OptionFile::trimComments(string &s)
{
  trimRight(s, "#;");
}

// Returns the next useful line, stripped of comments and whitespace
string OptionFile::getNextLine(ifstream &fin)
{
  string line("");
  
  getline(fin, line);
  trimComments(line);
  trim(line);
  line = lowercasequote(line); // lowercase except for inside quotes

/*  output << "DEBUG: " << line << endl;*/
  return line;
}

void OptionFile::parse(const string &buffer, string &key, string &value)
{
   // A key/value pair, separated by a '='

  size_t startpos = buffer.find_first_of("=");
  size_t endpos   = buffer.find_last_of("=");

  if( startpos == string::npos && startpos != endpos ) {
    throw BoutException("\tMissing (or multiple) '=' sign(s)\n\tLine: %s", buffer.c_str());
  }

  key = buffer.substr(0, startpos);
  value = buffer.substr(startpos+1);

  trim(key, " \t\"");
  trim(value, " \t\"");
  
/*  output << "DEBUG2: (" << key << "," << value << ")" << endl;*/

  if(key.empty() || value.empty()) {
    throw BoutException("\tEmpty key or value\n\tLine: %s", buffer.c_str());
  }
}
