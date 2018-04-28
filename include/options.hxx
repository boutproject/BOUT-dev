/*!************************************************************************
* Option hierarchy representation
*
* The Options class represents a tree structure of key-value settings.
* Provides get and set methods on these options.
*
* Internally, all quantities are stored as strings for simplicity
* and so that option file parsers don't have to do type conversion
* without knowing the type a priori.
*
* There is a singleton object "root" which contains the top-level
* options and allows access to all sub-sections
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

class Options;

#pragma once
#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include "bout_types.hxx"
#include "unused.hxx"

#include <map>
#include <string>
using std::string;

/// Class to represent hierarchy of options
/*!
 *
 *
 * Getting and setting values
 * --------------------------
 *
 * Each Options object represents a collection of key-value pairs
 * which can be used as a map.
 *
 *     Options options;
 *     options.set("key", 1.0, "code"); // Sets a key
 *
 *     int val;
 *     options.get("key", val, 0.0); // Sets val to 1.0
 *
 * If a variable has not been set then the default value is used
 *
 *     int other;
 *     options.get("otherkey", other, 2.0); // Sets other to 2.0 because "otherkey" not found
 *
 * Internally, all values are stored as strings, so conversion is performed silently:
 *
 *     options.set("value", "2.34", "here"); // Set a string
 *
 *     BoutReal value;
 *     options.get("value", value, 0.0); // Sets value to 2.34
 *
 * If a conversion cannot be done, then an exception is thrown
 *
 * Sections
 * --------
 *
 * Each Options object can also contain any number of sections, which are
 * themselves Options objects.
 *
 *     Options *section = options.getSection("section");
 *
 * This always succeeds; if the section does not exist then it is created.
 * Options also know about their parents:
 *
 *     section->getParent() == &options // Pointer to options object
 *
 * Root options object
 * -------------------
 *
 * For convenience, to avoid having to pass Options objects around everywhere,
 * there is a global singleton Options object which can be accessed with a static function
 *
 * Options *root = Options::getRoot();
 *
 * This is used to represent all the options passed to BOUT++ either in a file or on the
 * command line.
 *
 */
class Options {
public:
  /// Constructor. This is called to create the root object
  Options();
  
  /// Constructor used to create non-root objects
  ///
  /// @param[in] parent        Parent object
  /// @param[in] sectionName   Name of the section, including path from the root
  Options(Options *parent, string sectionName) : parent(parent), sectionName(sectionName) {};
  
  /// Destructor
  ~Options();

  /// Get a pointer to the only root instance (singleton)
  static Options* getRoot();

  /*!
   * Free all memory
   */
  static void cleanup();

  // Setting options
  template<typename T> void forceSet(const string &key, T t, const string &source=""){
    set(key,t,source,true);
  }
  void set(const string &key, const int &val, const string &source="", bool force=false);
  void set(const string &key, BoutReal val, const string &source="", bool force=false);
  void set(const string &key, const bool &val, const string &source="", bool force=false);
  void set(const string &key, const char * &val, const string &source="", bool force=false) {
    _set(key,val,source,force);
  }
  void set(const string &key, const string &val, const string &source = "", bool force=false) {
    _set(key,val,source,force);
  };

  /// Set a string with a char* array. This converts to std::string
  /// rather than allow an implicit conversion to bool
  void set(const string &key, const char *val, const string &source = "") {
    set(key, string(val), source);
  }

  /*!
   * Test if a key is set to a value
   *
   */
  bool isSet(const string &key);

  // Getting options
  void get(const string &key, int &val, int def);
  void get(const string &key, BoutReal &val, BoutReal def);
  void get(const string &key, bool &val, bool def);
  void get(const string &key, string &val, const string &def);

  // This is a temporary replacement for 4-argument get
  // and will be removed in the next major release
  template<typename T, typename U>
  void get(const string &key, T &val, U def, bool UNUSED(log)) {
    get(key, val, def);
  }

  /// Creates new section if doesn't exist
  Options* getSection(const string &name);
  Options* getParent() {return parent;}

  /*!
   * Print string representation of this object and all sections in a tree structure
   */
  string str();

  /// Print the options which haven't been used
  void printUnused();

  /// clean the cache of parsed options
  static void cleanCache();

  /*!
   * Class used to store values, together with
   * information about their origin and usage
   */
  struct OptionValue {
    string value;
    string source;     // Source of the setting
    bool used;         // Set to true when used
  };

  /// Read-only access to internal options and sections
  /// to allow iteration over the tree
  const std::map<string, OptionValue>& values() const {return options;}
  const std::map<string, Options*>& subsections() const {return sections;}

 private:
  static Options *root; ///< Only instance of the root section

  Options *parent;
  string sectionName; // section name (if any), for logging only

  std::map<string, OptionValue> options;
  std::map<string, Options*> sections;

  void _set(const string &key, const string &val, const string &source, bool force);
};

/// Define for reading options which passes the variable name
#define OPTION(options, var, def)  \
  options->get(#var, var, def)

#define OPTION2(options, var1, var2, def){ \
    options->get(#var1, var1, def);  \
    options->get(#var2, var2, def);}

#define OPTION3(options, var1, var2, var3, def){  \
    options->get(#var1, var1, def);               \
    options->get(#var2, var2, def);               \
    options->get(#var3, var3, def);}

#define OPTION4(options, var1, var2, var3, var4, def){ \
    options->get(#var1, var1, def);               \
    options->get(#var2, var2, def);               \
    options->get(#var3, var3, def);               \
    options->get(#var4, var4, def);}

#define OPTION5(options, var1, var2, var3, var4, var5, def){ \
    options->get(#var1, var1, def);                      \
    options->get(#var2, var2, def);                      \
    options->get(#var3, var3, def);                      \
    options->get(#var4, var4, def);                      \
    options->get(#var5, var5, def);}

#define OPTION6(options, var1, var2, var3, var4, var5, var6, def){ \
    options->get(#var1, var1, def);                               \
    options->get(#var2, var2, def);                               \
    options->get(#var3, var3, def);                               \
    options->get(#var4, var4, def);                               \
    options->get(#var5, var5, def);                               \
    options->get(#var6, var6, def);}

#define VAROPTION(options, var, def) {					\
    if (options->isSet(#var)){						\
      options->get(#var, var, def);					\
    } else {								\
      Options::getRoot()->getSection("all")->get(#var, var, def);	\
    }}									\

#endif // __OPTIONS_H__
