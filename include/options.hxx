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
#include "output.hxx"

#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <utility>
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
 *     options["key"] = 1.0;
 *
 *     int val = options["key"]; // Sets val to 1 
 *
 * It can be useful to record where values came from by adding a source label:
 * 
 *     options["key"].setTo(1.0, "code"); // Sets a value with source "code"
 *
 * Equivalent interfaces are:
 * 
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
 *     Options &section = options["section"];
 * 
 * or
 *     Options *section = options.getSection("section");
 *
 * This always succeeds; if the section does not exist then it is created.
 * Options also know about their parents:
 *
 *     Options &parent = section.parent();
 *     
 * or
 * 
 *     section->getParent() == &options // Pointer to options object
 *
 * Root options object
 * -------------------
 *
 * For convenience, to avoid having to pass Options objects around everywhere,
 * there is a global singleton Options object which can be accessed with a static function
 *
 *    Options &root = Options::root();
 * 
 * or 
 *
 *    Options *root = Options::getRoot();
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
  Options(Options *parent_instance, string full_name)
      : parent_instance(parent_instance), full_name(std::move(full_name)){};

  /// Destructor
  ~Options();

  /// Get a reference to the only root instance
  static Options &root();
  
  /// Free all memory
  static void cleanup();

  /// Get a sub-section or value
  /// Note: Using this makes this object a section
  Options& operator[](const string &name);

  /// Assignment from any type T
  /// Note: Using this makes this object a value.
  ///
  /// Tries to stream the input to a std::stringstream
  /// and then saves the resulting string.
  template <typename T>
  const T & operator=(const T &inputvalue) {
    // Convert to string
    std::stringstream ss;
    ss << inputvalue;

    // Set the internal value. This will perform
    // checks to ensure that this Options object
    // is a value and not a section.
    _set(ss.str(), "", false);
    return inputvalue;
  }

  // Setting options
  template<typename T> void forceSet(const string &key, T t, const string &source=""){
    set(key,t,source,true);
  }
  
  void setTo(const int &val, const string &source="", bool force=false);
  void setTo(BoutReal val, const string &source="", bool force=false);
  void setTo(bool val, const string &source="", bool force=false);


  void setTo(const char *val, const string &source="", bool force=false) {
    _set(val,source,force);
  }
  void setTo(const string &val, const string &source = "", bool force=false) {
    _set(val,source,force);
  };
  
  /// Test if a key is set by the user.
  /// Values set via default values are ignored.
  bool isSet();

  // Getting options

  /// Cast operator, which allows this class to be
  /// assigned to type T
  ///
  template <typename T> operator T() { return get<T>(); }

  /// Get the value as a specified type
  /// If there is no value then an exception is thrown
  template <typename T> T get() {
    if (!isSet()) {
      throw BoutException("Option %s has no value", full_name.c_str());
    }
    
    T val;
    std::stringstream ss(value.value);
    ss >> val;
    
    // Check if the parse failed
    if (ss.fail()) {
      throw BoutException("Option %s could not be parsed", full_name.c_str());
    }

    // Check if there are characters remaining
    std::string remainder;
    std::getline(ss, remainder);
    for (const char &ch : remainder) {
      if (!std::isspace(static_cast<unsigned char>(ch))) {
        // Meaningful character not parsed
        throw BoutException("Option %s could not be parsed", full_name.c_str());
      }
    }
    
    // Mark this option as used
    value.used = true;

    output_info << "\tOption " << full_name  << " = " << val;
    if (!value.source.empty()) {
      // Specify the source of the setting
      output_info << " (" << value.source << ")";
    }
    output_info << endl;

    return val;
  }

  int get(int def);
  BoutReal get(BoutReal def);
  bool get(bool def);
  std::string get(const std::string &def);

  Options &parent() { return *parent_instance; }

  //////////////////////////////////////
  // Backward-compatible interface
  
  /// Get a pointer to the only root instance (singleton)
  static Options* getRoot() { return &root(); }

  template<typename T>
  void set(const string &key, T val, const string &source = "", bool force = false) {
    (*this)[key].setTo(val, source, force);
  }
  
  /*!
   * Test if a key is set by the user.
   * Values set via default values are ignored.
   */
  bool isSet(const string &key) { return (*this)[key].isSet(); }

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
  Options* getSection(const string &name) { return &(*this)[name]; }
  Options* getParent() {return parent_instance;}

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
  std::map<string, OptionValue> values() const;
  std::map<string, Options*> subsections() const;

 private:
  
  static Options *root_instance; ///< Only instance of the root section

  Options *parent_instance;
  string full_name; // full path name for logging only

  /// An Option object can be a section (node), a value (leaf), or empty (undecided)
  enum class OptionType { empty, section, value };

  OptionType type = OptionType::empty; ///< The type of this Option object
  
  std::string name; ///< If a section then this is the section name
  std::map<string, Options*> children; ///< If a section then has children
  OptionValue value; ///< If a value
  
  void _set(const string &val, const string &source, bool force);
};

/// Specialised assignment operator
template <>
inline const BoutReal & Options::operator=<BoutReal>(const BoutReal &inputvalue) {
  std::stringstream ss;
  // Make sure the precision is large enough to hold a BoutReal
  ss << std::scientific << std::setprecision(17) << inputvalue;
  _set(ss.str(), "", false);
  return inputvalue;
}

/// Specialised get routines
template <>
inline std::string Options::get<std::string>() {
  if (!isSet()) {
    throw BoutException("Option %s has no value", full_name.c_str());
  }
    
  // Mark this option as used
  value.used = true;
  
  output_info << "\tOption " << full_name  << " = " << value.value;
  if (!value.source.empty()) {
    // Specify the source of the setting
    output_info << " (" << value.source << ")";
  }
  output_info << endl;
  
  return value.value;
}

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
