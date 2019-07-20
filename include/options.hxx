/*!************************************************************************
* Option hierarchy representation
*
* The Options class represents a tree structure of key-value settings.
* Provides get and set methods on these options.
*
* Internally, values are stored in a variant. Conversion of types
* is handled internally, so this is transparent to the user. If desired,
* the user can directly access the "value" member which is the variant.
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
#include "utils.hxx"
#include "bout/sys/variant.hxx"
#include "bout/sys/type_name.hxx"
#include "bout/deprecated.hxx"
#include "field2d.hxx"
#include "field3d.hxx"

#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <utility>
#include <cmath>

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
 * Conversion is performed silently:
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
  Options() = default;
  
  /// Constructor used to create non-root objects
  ///
  /// @param[in] parent        Parent object
  /// @param[in] sectionName   Name of the section, including path from the root
  Options(Options *parent_instance, std::string full_name)
      : parent_instance(parent_instance), full_name(std::move(full_name)){};

  /// Copy constructor
  Options(const Options& other);

  ~Options() = default;

  /// Get a reference to the only root instance
  static Options &root();
  
  /// Free all memory
  static void cleanup();

  /// The type used to store values
  using ValueType =
      bout::utils::variant<bool, int, BoutReal, std::string, Field2D, Field3D,
                           Array<BoutReal>, Matrix<BoutReal>, Tensor<BoutReal>>;

  /// The type used to store attributes
  /// Extends the variant class so that cast operator can be implemented
  /// and assignment operator overloaded
  ///
  /// Note: Due to default initialisation rules, if an attribute
  /// is used without being set, it will be false, 0, 0.0 and
  /// throw std::bad_cast if cast to std::string
  /// 
  class AttributeType : public bout::utils::variant<bool, int, BoutReal, std::string> {
  public:
    using Base = bout::utils::variant<bool, int, BoutReal, std::string>;

    /// Constructor
    AttributeType() = default;
    /// Copy constructor
    AttributeType(const AttributeType& other) = default;
    /// Move constructor
    AttributeType(AttributeType&& other) : Base(std::move(other)) {}

    /// Destructor
    ~AttributeType() = default;

    /// Assignment operator, including move assignment
    using Base::operator=;

    /// Assignment from const char*
    AttributeType& operator=(const char* str) {
      operator=(std::string(str));
      return *this;
    }

    /// Cast operator, which allows this class to be
    /// assigned to type T
    /// This will throw std::bad_cast if it can't be done
    template <typename T> operator T() const { return as<T>(); }

    /// Get the value as a specified type
    /// This will throw std::bad_cast if it can't be done
    template <typename T>
    T as() const {
      return bout::utils::variantStaticCastOrThrow<Base, T>(*this);
    }
  };

  /// The value stored
  ValueType value;
  
  /// A collection of attributes belonging to the value
  /// Special attributes:
  ///  - time_dimension   [string] If this is set then changes to the value
  ///                     do not need to be forced. The string will be used
  ///                     when writing the output as the name of the time
  ///                     dimension (unlimited first dimension in NetCDF files).
  ///
  ///  - source           [string] Describes where the value came from
  ///                     e.g. a file name, or "default".
  /// 
  ///  - type             [string] The type the Option is converted to
  ///                     when used.
  /// 
  ///  - doc              [string] Documentation, describing what the variable does
  ///
  std::map<std::string, AttributeType> attributes;
  
  /// Get a sub-section or value
  ///
  /// Example:
  ///
  /// Options parent;
  /// auto child  = parent["child"];
  ///
  /// parent is now a section.
  Options& operator[](const std::string &name);
  // Prevent ambiguous overload resolution due to operator T() below
  Options& operator[](const char *name) { return (*this)[std::string(name)]; }

  /// Get a sub-section or value
  /// If this object is not a section, or if name
  /// is not a child, then a BoutException will be thrown
  const Options& operator[](const std::string &name) const;
  const Options& operator[](const char *name) const { return (*this)[std::string(name)]; }

  /// Assignment from any type T
  /// Note: Using this makes this object a value.
  ///
  /// Tries to stream the input to a std::stringstream
  /// and then saves the resulting string.
  template <typename T>
  T operator=(T inputvalue) {
    assign<T>(inputvalue);
    return inputvalue;
  }

  /// Copy assignment
  ///
  /// This replaces the value, attributes and all children
  ///
  /// Note that if only the value is desired, then that can be copied using
  /// the value member directly e.g. option2.value = option1.value;
  ///
  Options& operator=(const Options& other);
  
  /// Assign a value to the option.
  /// This will throw an exception if already has a value
  ///
  /// Example:
  ///
  /// Options option;
  /// option["test"].assign(42, "some source");
  ///
  /// Note: Specialised versions for types stored in ValueType
  template<typename T>
  void assign(T val, const std::string source="") {
    std::stringstream ss;
    ss << val;
    _set(ss.str(), source, false);
  }
  
  /// Force to a value
  /// Overwrites any existing setting
  template<typename T>
  void force(T val, const std::string source = "") {
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
  /// An optional argument is an object which the result should be similar to.
  /// The main use for this is in Field2D and Field3D specialisations,
  /// where the Mesh and cell location are taken from this input.
  /// 
  template <typename T>
  T as(const T& UNUSED(similar_to) = {}) const {
    if (!is_value) {
      throw BoutException("Option %s has no value", full_name.c_str());
    }

    T val;
    
    // Try casting. This will throw std::bad_cast if it can't be done
    try {
      val = bout::utils::variantStaticCastOrThrow<ValueType, T>(value);
    } catch (const std::bad_cast &e) {
      // If the variant is a string then we may be able to parse it
      
      if (bout::utils::holds_alternative<std::string>(value)) {
        std::stringstream ss(bout::utils::get<std::string>(value));
        ss >> val;
        
        // Check if the parse failed
        if (ss.fail()) {
          throw BoutException("Option %s could not be parsed ('%s')", full_name.c_str(),
                              bout::utils::variantToString(value).c_str());
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
      } else {
        // Another type which can't be casted
        throw BoutException("Option %s could not be converted to type %s",
                            full_name.c_str(), typeid(T).name());
      }
    }
    
    // Mark this option as used
    value_used = true; // Note this is mutable

    output_info << "\tOption " << full_name  << " = " << val;
    if (attributes.count("source")) {
      // Specify the source of the setting
      output_info << " (" << bout::utils::variantToString(attributes.at("source")) << ")";
    }
    output_info << endl;

    return val;
  }

  /// Get the value of this option. If not found,
  /// set to the default value
  template <typename T> T withDefault(T def) {

    // Set the type
    attributes["type"] = bout::utils::typeName<T>();
    
    if (!is_value) {
      // Option not found
      assign(def, DEFAULT_SOURCE);
      value_used = true; // Mark the option as used

      output_info << _("\tOption ") << full_name << " = " << def << " (" << DEFAULT_SOURCE
                  << ")" << std::endl;
      return def;
    }
    T val = as<T>(def);
    // Check if this was previously set as a default option
    if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
      // Check that the default values are the same
      if (!similar(val, def)) {
        throw BoutException("Inconsistent default values for '%s': '%s' then '%s'",
                            full_name.c_str(), bout::utils::variantToString(value).c_str(), toString(def).c_str());
      }
    }
    return val;
  }

  /// Overloaded version for const char*
  /// Note: Different from template since return type is different to input
  std::string withDefault(const char* def) {
    return withDefault<std::string>(std::string(def));
  }
  
  /// Overloaded version to copy from another option
  Options& withDefault(const Options& def) {
    // if def is a section, then it does not make sense to try to use it as a default for
    // a value
    ASSERT0(def.is_value);

    if (!is_value) {
      // Option not found
      *this = def;

      output_info << _("\tOption ") << full_name << " = " << def.full_name << " ("
                  << DEFAULT_SOURCE << ")" << std::endl;
    } else {
      // Check if this was previously set as a default option
      if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
        // Check that the default values are the same
        if (!similar(bout::utils::variantToString(value),
                     bout::utils::variantToString(def.value))) {
          throw BoutException("Inconsistent default values for '%s': '%s' then '%s'",
              full_name.c_str(), bout::utils::variantToString(value).c_str(),
              bout::utils::variantToString(def.value).c_str());
        }
      }
    }
    return *this;
  }

  /// Get the value of this option. If not found,
  /// return the default value but do not set
  template <typename T> T withDefault(T def) const {
    if (!is_value) {
      // Option not found
      output_info << _("\tOption ") << full_name << " = " << def << " (" << DEFAULT_SOURCE
                  << ")" << std::endl;
      return def;
    }
    T val = as<T>(def);
    // Check if this was previously set as a default option
    if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
      // Check that the default values are the same
      if (!similar(val, def)) {
        throw BoutException("Inconsistent default values for '%s': '%s' then '%s'",
                            full_name.c_str(), bout::utils::variantToString(value).c_str(), toString(def).c_str());
      }
    }
    return val;
  }

  /// Allow the user to override defaults set later, also used by the
  /// BOUT_OVERRIDE_DEFAULT_OPTION.
  template <typename T> T overrideDefault(T def) {

    // Set the type
    attributes["type"] = bout::utils::typeName<T>();

    if (!is_value) {
      // Option not found
      assign(def, "user_default");
      value_used = true; // Mark the option as used
      is_value = true; // Prevent this default being replaced by setDefault()

      output_info << _("\tOption ") << full_name << " = " << def << " (" << "user_default"
                  << ")" << std::endl;
      return def;
    }

    // Return value of this option as type 'T'
    return as<T>(def);
  }

  /// Get the parent Options object
  Options &parent() {
    if (parent_instance == nullptr) {
      throw BoutException("Option %s has no parent", full_name.c_str());
    }
    return *parent_instance;
  }

  /// Equality operator
  /// Converts to the same type and compares
  /// This conversion may fail, throwing std::bad_cast
  template<typename T>
  bool operator==(const T& other) const {
    return as<T>() == other;
  }

  /// Overloaded equality operator for literal strings
  bool operator==(const char* other) const;
  
  /// Comparison operator
  template<typename T>
  bool operator<(const T& other) const {
    return as<T>() < other;
  }

  /// Overloaded comparison operator for literal strings
  bool operator<(const char* other) const;
  
  //////////////////////////////////////
  // Backward-compatible interface
  
  /// Get a pointer to the only root instance (singleton)
  static Options* getRoot() { return &root(); }

  template<typename T>
  void set(const std::string &key, T val, const std::string &source = "", bool force = false) {
    if (force) {
      (*this)[key].force(val, source);
    } else {
      (*this)[key].assign(val, source);
    }
  }
  
  // Setting options
  template<typename T> void forceSet(const std::string &key, T t, const std::string &source=""){
    (*this)[key].force(t,source);
  }
  
  /*!
   * Test if a key is set by the user.
   * Values set via default values are ignored.
   */
  bool isSet(const std::string &key) const {
    // Note using operator[] here would result in exception if key does not exist
    if (!is_section)
      return false;
    auto it = children.find(lowercase(key));
    if (it == children.end())
      return false;
    return it->second.isSet();
  }

  /// Get options, passing in a reference to a variable
  /// which will be set to the requested value.
  template<typename T, typename U>
  void get(const std::string &key, T &val, U def, bool UNUSED(log)=false) {
    val = (*this)[key].withDefault<T>(def);
  }
  template<typename T, typename U>
  void get(const std::string &key, T &val, U def, bool UNUSED(log)=false) const {
    val = (*this)[key].withDefault<T>(def);
  }

  /// Creates new section if doesn't exist
  Options* getSection(const std::string &name) { return &(*this)[name]; }
  const Options* getSection(const std::string &name) const { return &(*this)[name]; }
  Options* getParent() const {return parent_instance;}
  
  //////////////////////////////////////
  
  /*!
   * Print string representation of this object and all sections in a tree structure
   */
  std::string str() const {return full_name;}

  /// Print the options which haven't been used
  void printUnused() const;

  /// clean the cache of parsed options
  static void cleanCache();
  
  /*!
   * Class used to store values, together with
   * information about their origin and usage
   */
  struct OptionValue {
    std::string value;
    std::string source;     // Source of the setting
    mutable bool used = false;  // Set to true when used

    /// This constructor needed for map::emplace
    /// Can be removed in C++17 with map::insert and brace initialisation
    OptionValue(std::string value, std::string source, bool used)
        : value(std::move(value)), source(std::move(source)), used(used) {}
  };

  /// Read-only access to internal options and sections
  /// to allow iteration over the tree
  using ValuesMap = std::map<std::string, OptionValue>;
  DEPRECATED(ValuesMap values() const);
  std::map<std::string, const Options*> subsections() const;

  const std::map<std::string, Options>& getChildren() const {
    return children;
  }

  bool isValue() const {
    return is_value;
  }
  bool isSection(const std::string& name = "") const;
  
  /// If the option value has been used anywhere
  bool valueUsed() const { return value_used; }


  /// Set a documentation string as an attribute "doc"
  /// Returns a reference to this, to allow chaining
  Options& doc(const std::string& docstring) {
    attributes["doc"] = docstring;
    return *this;
  }
  
 private:
  
  /// The source label given to default values
  static const std::string DEFAULT_SOURCE;
  
  static Options *root_instance; ///< Only instance of the root section

  Options *parent_instance {nullptr};
  std::string full_name; // full path name for logging only

  /// An Option object can be a section and/or a value, or neither (empty)

  bool is_section = false; ///< Is this Options object a section?
  std::map<std::string, Options> children; ///< If a section then has children

  bool is_value = false; ///< Is this Options object a value?
  mutable bool value_used = false; ///< Record whether this value is used
  
  template <typename T>
  void _set(T val, std::string source, bool force) {
    // If already set, and not time evolving then check for changing values
    // If a variable has a "time_dimension" attribute then it is assumed
    // that updates to the value is ok and don't need to be forced.
    if (isSet() && (attributes.find("time_dimension") == attributes.end())) {
      // Check if current value the same as new value
      if (!bout::utils::variantEqualTo(value, val)) {
        if (force or !bout::utils::variantEqualTo(attributes["source"], source)) {
          output_warn << _("\tOption ") << full_name << " = "
                      << bout::utils::variantToString(value) << " ("
                      << bout::utils::variantToString(attributes["source"])
                      << _(") overwritten with:") << "\n"
                      << "\t\t" << full_name << " = " << toString(val) << " (" << source
                      << ")\n";
        } else {
          throw BoutException(
              _("Options: Setting a value from same source (%s) to new value "
                "'%s' - old value was '%s'."),
              source.c_str(), toString(val).c_str(),
              bout::utils::variantToString(value).c_str());
        }
      }
    }

    value = std::move(val);
    attributes["source"] = std::move(source);
    value_used = false;
    is_value = true;
  }
  
  /// Tests if two values are similar. 
  template <typename T> bool similar(T a, T b) const { return a == b; }
};

// Specialised assign methods for types stored in ValueType
template<> inline void Options::assign<>(bool val, const std::string source) { _set(val, source, false); }
template<> inline void Options::assign<>(int val, const std::string source) { _set(val, source, false); }
template<> inline void Options::assign<>(BoutReal val, const std::string source) { _set(val, source, false); }
template<> inline void Options::assign<>(std::string val, const std::string source) { _set(val, source, false); }
// Note: const char* version needed to avoid conversion to bool
template<> inline void Options::assign<>(const char *val, const std::string source) { _set(std::string(val), source, false);}
// Note: Field assignments don't check for previous assignment (always force)
template<> void Options::assign<>(Field2D val, const std::string source);
template<> void Options::assign<>(Field3D val, const std::string source);
template<> void Options::assign<>(Array<BoutReal> val, const std::string source);
template<> void Options::assign<>(Matrix<BoutReal> val, const std::string source);
template<> void Options::assign<>(Tensor<BoutReal> val, const std::string source);

/// Specialised similar comparison methods
template <> inline bool Options::similar<BoutReal>(BoutReal a, BoutReal b) const { return fabs(a - b) < 1e-10; }

/// Specialised as routines
template <> std::string Options::as<std::string>(const std::string& similar_to) const;
template <> int Options::as<int>(const int& similar_to) const;
template <> BoutReal Options::as<BoutReal>(const BoutReal& similar_to) const;
template <> bool Options::as<bool>(const bool& similar_to) const;
template <> Field2D Options::as<Field2D>(const Field2D& similar_to) const;
template <> Field3D Options::as<Field3D>(const Field3D& similar_to) const;

/// Define for reading options which passes the variable name
#define OPTION(options, var, def)  \
  pointer(options)->get(#var, var, def)

#define OPTION2(options, var1, var2, def){ \
    pointer(options)->get(#var1, var1, def);  \
    pointer(options)->get(#var2, var2, def);}

#define OPTION3(options, var1, var2, var3, def){  \
    pointer(options)->get(#var1, var1, def);               \
    pointer(options)->get(#var2, var2, def);               \
    pointer(options)->get(#var3, var3, def);}

#define OPTION4(options, var1, var2, var3, var4, def){ \
    pointer(options)->get(#var1, var1, def);               \
    pointer(options)->get(#var2, var2, def);               \
    pointer(options)->get(#var3, var3, def);               \
    pointer(options)->get(#var4, var4, def);}

#define OPTION5(options, var1, var2, var3, var4, var5, def){ \
    pointer(options)->get(#var1, var1, def);                      \
    pointer(options)->get(#var2, var2, def);                      \
    pointer(options)->get(#var3, var3, def);                      \
    pointer(options)->get(#var4, var4, def);                      \
    pointer(options)->get(#var5, var5, def);}

#define OPTION6(options, var1, var2, var3, var4, var5, var6, def){ \
    pointer(options)->get(#var1, var1, def);                               \
    pointer(options)->get(#var2, var2, def);                               \
    pointer(options)->get(#var3, var3, def);                               \
    pointer(options)->get(#var4, var4, def);                               \
    pointer(options)->get(#var5, var5, def);                               \
    pointer(options)->get(#var6, var6, def);}

#define VAROPTION(options, var, def) {					\
    if (pointer(options)->isSet(#var)){						\
      pointer(options)->get(#var, var, def);					\
    } else {								\
      Options::getRoot()->getSection("all")->get(#var, var, def);	\
    }}									\

#endif // __OPTIONS_H__
