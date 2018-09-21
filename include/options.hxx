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
 *     
 *     // Set values
 *     options["key"] = 1.0;
 *
 *     // Get values. Throws BoutException if not found
 *     int val = options["key"]; // Sets val to 1 
 *
 *     // Return as specified type. Throws BoutException if not found
 *     BoutReal var = options["key"].as<BoutReal>();
 *
 *     // A default value can be used if key is not found
 *     BoutReal value = options["pi"].withDefault(3.14);
 *    
 *     // Assign value with source label. Throws if already has a value from same source
 *     options["newkey"].assign(1.0, "some source");
 *
 *     // Force assign a new value
 *     options["newkey"].force(2.0, "some source");
 *
 * A legacy interface is also supported:
 * 
 *     options.set("key", 1.0, "code"); // Sets a key from source "code"
 *
 *     int val;
 *     options.get("key", val, 0); // Sets val to 1, default to 0 if not found
 *
 * If a variable has not been set then the default value is used
 *
 *     int other;
 *     options.get("otherkey", other, 2.0); // Sets other to 2 because "otherkey" not found
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
 * which can be nested:
 *
 *     options["section"]["subsection"]["value"] = 3;
 *
 * This always succeeds; if the section does not exist then it is created.
 *
 * The legacy interface uses pointers:
 *
 *     Options *section = options.getSection("section");
 *
 * e.g.
 *     options->getSection("section")->getSection("subsection")->set("value", 3);
 * 
 * Options also know about their parents:
 *
 *     Options &parent = section.parent();
 *     
 * or
 * 
 *     Options *parent = section->getParent();
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
  
  /// Get a reference to the only root instance
  static Options &root();
  
  /// Free all memory
  static void cleanup();

  /// Get a sub-section or value
  ///
  /// Example:
  ///
  /// Options parent;
  /// auto child  = parent["child"];
  ///
  /// parent is now a section.
  Options& operator[](const string &name);
  // Prevent ambiguous overload resolution due to operator T() below
  Options& operator[](const char *name) { return (*this)[std::string(name)]; }

  /// Get a sub-section or value
  /// If this object is not a section, or if name
  /// is not a child, then a BoutException will be thrown
  const Options& operator[](const string &name) const;
  const Options& operator[](const char *name) const { return (*this)[std::string(name)]; }

  /// Assignment from any type T
  /// Note: Using this makes this object a value.
  ///
  /// Tries to stream the input to a std::stringstream
  /// and then saves the resulting string.
  template <typename T>
  T operator=(T inputvalue) {
    // Convert to string
    std::stringstream ss;
    ss << inputvalue;

    // Set the internal value.
    _set(ss.str(), "", false);
    return inputvalue;
  }

  /// Assign a value to the option.
  /// This will throw an exception if already has a value
  ///
  /// Example:
  ///
  /// Options option;
  /// option["test"].assign(42, "some source");
  ///
  template<typename T>
  void assign(T val, const string &source="") {
    std::stringstream ss;
    ss << val;
    _set(ss.str(), source, false);
  }
  
  /// Force to a value
  /// Overwrites any existing setting
  template<typename T>
  void force(T val, const string &source = "") {
    is_value = false; // Invalidates any existing setting
    assign(val, source);
  }
  
  /// Test if a key is set by the user.
  /// Values set via default values are ignored.
  bool isSet() const;

  // Getting options

  /// Cast operator, which allows this class to be
  /// assigned to type T
  ///
  /// Example:
  ///
  /// Options option;
  /// option["test"] = 2.0;
  /// int value = option["test"];
  ///
  template <typename T> operator T() const { return as<T>(); }

  /// Get the value as a specified type
  /// If there is no value then an exception is thrown
  /// Note there are specialised versions of this template
  /// for some types.
  ///
  /// Example:
  ///
  /// Options option;
  /// option["test"] = 2.0;
  /// int value = option["test"].as<int>();
  ///
  template <typename T> T as() const {
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
    value.used = true; // Note this is mutable

    output_info << "\tOption " << full_name  << " = " << val;
    if (!value.source.empty()) {
      // Specify the source of the setting
      output_info << " (" << value.source << ")";
    }
    output_info << endl;

    return val;
  }

  /// Get the value of this option. If not found,
  /// set to the default value
  template <typename T> T withDefault(T def) {
    if (!is_value) {
      // Option not found
      assign(def, "default");
      value.used = true; // Mark the option as used

      output_info << "\tOption " << full_name << " = " << def << " (default)"
                  << std::endl;
      return def;
    }
    return as<T>();
  }

  /// Get the value of this option. If not found,
  /// return the default value but do not set
  template <typename T> T withDefault(T def) const {
    if (!is_value) {
      // Option not found
      output_info << "\tOption " << full_name << " = " << def << " (default)"
                  << std::endl;
      return def;
    }
    return as<T>();
  }

  /// Get the parent Options object
  Options &parent() {
    if (parent_instance == nullptr) {
      throw BoutException("Option %s has no parent", full_name.c_str());
    }
    return *parent_instance;
  }
  
  //////////////////////////////////////
  // Backward-compatible interface
  
  /// Get a pointer to the only root instance (singleton)
  static Options* getRoot() { return &root(); }

  template<typename T>
  void set(const string &key, T val, const string &source = "", bool force = false) {
    if (force) {
      (*this)[key].force(val, source);
    } else {
      (*this)[key].assign(val, source);
    }
  }
  
  // Setting options
  template<typename T> void forceSet(const string &key, T t, const string &source=""){
    (*this)[key].force(t,source);
  }
  
  /*!
   * Test if a key is set by the user.
   * Values set via default values are ignored.
   */
  bool isSet(const string &key) { return (*this)[key].isSet(); }

  /// Get options, passing in a reference to a variable
  /// which will be set to the requested value.
  template<typename T, typename U>
  void get(const string &key, T &val, U def, bool UNUSED(log)=false) {
    val = (*this)[key].withDefault<T>(def);
  }

  /// Creates new section if doesn't exist
  Options* getSection(const string &name) { return &(*this)[name]; }
  Options* getParent() {return parent_instance;}
  
  //////////////////////////////////////
  
  /*!
   * Print string representation of this object and all sections in a tree structure
   */
  string str() const {return full_name;}

  /// Print the options which haven't been used
  void printUnused() const;

  /// clean the cache of parsed options
  static void cleanCache();

  /*!
   * Class used to store values, together with
   * information about their origin and usage
   */
  struct OptionValue {
    string value;
    string source;     // Source of the setting
    mutable bool used = false;  // Set to true when used
  };

  /// Read-only access to internal options and sections
  /// to allow iteration over the tree
  std::map<string, OptionValue> values() const;
  std::map<string, const Options*> subsections() const;

 private:
  
  static Options *root_instance; ///< Only instance of the root section

  Options *parent_instance;
  string full_name; // full path name for logging only

  /// An Option object can be a section and/or a value, or neither (empty)

  bool is_section = false; ///< Is this Options object a section?
  std::map<string, Options> children; ///< If a section then has children

  bool is_value = false; ///< Is this Options object a value?
  OptionValue value{}; ///< If a value
  
  void _set(string val, string source, bool force);
};

/// Specialised assignment operator
template <>
inline BoutReal Options::operator=<BoutReal>(BoutReal inputvalue) {
  std::stringstream ss;
  // Make sure the precision is large enough to hold a BoutReal
  ss << std::scientific << std::setprecision(17) << inputvalue;
  _set(ss.str(), "", false);
  return inputvalue;
}

/// Specialised as routines
template <>
inline std::string Options::as<std::string>() const {
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

// Specialised assign methods
template<> void Options::assign<bool>(bool val, const string &source);
template<> void Options::assign<BoutReal>(BoutReal val, const string &source);

// Note: const char* version needed to avoid conversion to bool
template<>
inline void Options::assign<const char *>(const char *val, const string &source) {
  _set(val,source,false);
}
template<>
inline void Options::assign<std::string>(std::string val, const string &source) {
  _set(val,source,false);
};

// Specialised withDefault methods
template<> int Options::withDefault<int>(int def);
template<> BoutReal Options::withDefault<BoutReal>(BoutReal def);
template<> bool Options::withDefault<bool>(bool def);
template<> std::string Options::withDefault<std::string>(std::string def);

template<> int Options::withDefault<int>(int def) const;
template<> BoutReal Options::withDefault<BoutReal>(BoutReal def) const;
template<> bool Options::withDefault<bool>(bool def) const;
template<> std::string Options::withDefault<std::string>(std::string def) const;

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
