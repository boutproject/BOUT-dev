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
* Copyright 2010-2024 BOUT++ contributors
*
* Contact: Ben Dudson, dudson2@llnl.gov
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
#ifndef OPTIONS_H
#define OPTIONS_H

#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldperp.hxx"
#include "bout/output.hxx"
#include "bout/sys/type_name.hxx"
#include "bout/sys/variant.hxx"
#include "bout/traits.hxx"
#include "bout/unused.hxx"
#include "bout/utils.hxx"

#include <fmt/core.h>

#include <cmath>
#include <functional>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>

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
 * Copying options
 * ---------------
 *
 * The copy constructor and copy assignment operator are deleted, so
 * this is a compile-time error:
 *
 *     Options options2 = options1["value"];
 *
 * This is because it's ambiguous what is meant, and because accidental copies
 * were a frequent source of hard-to-understand bugs. Usually a reference is
 * intended, rather than a copy:
 *
 *     Options& ref = options1["value"];
 *
 * so that changes to `ref` or its children are reflected in `options1`.
 * If the intent is to copy the value of the option, then just copy that:
 *
 *     option2.value = options1["value"].value;
 *
 * If a full deep copy of the option, its attributes and children
 * (recursively) is really desired, then use the `copy()` method:
 *
 *     Options options2 = options1["value"].copy();
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
  Options(Options* parent_instance, std::string full_name)
      : parent_instance(parent_instance), full_name(std::move(full_name)){};

  /// Initialise with a value
  /// These enable Options to be constructed using initializer lists
  template <typename T>
  Options(T value) {
    assign<T>(value);
  }

  /// The type used to store values
  using ValueType =
      bout::utils::variant<bool, int, BoutReal, std::string, Field2D, Field3D, FieldPerp,
                           Array<BoutReal>, Matrix<BoutReal>, Tensor<BoutReal>>;

  /// A tree representation with leaves containing ValueType.
  /// Used to construct Options from initializer lists.
  ///
  /// Note: Either there are children OR value is assigned
  struct CopyableOptions {
    template <typename T>
    CopyableOptions(T value) : value(std::move(value)) {}

    /// Special case for char*, which can otherwise become cast to bool
    CopyableOptions(const char* value) : value(std::string(value)) {}

    CopyableOptions(
        std::initializer_list<std::pair<std::string, CopyableOptions>> children)
        : children(children) {}
    ValueType value;
    std::initializer_list<std::pair<std::string, CopyableOptions>> children;
  };

  /// Type of initializer_list that can be used to create Options
  /// This is a workaround for initializer_lists not being movable.
  using InitializerList = std::initializer_list<std::pair<std::string, CopyableOptions>>;

  /// Construct with a nested initializer list
  /// This allows Options trees to be constructed using a mix of types.
  ///
  /// Example:  { {"key1", 42}, {"key2", field} }
  ///
  /// Note: Options doesn't have a copy constructor, and initializer lists
  ///       don't play nicely with uncopyable types. Instead, we create
  ///       a tree of CopyableOptions and then move.
  Options(InitializerList values, Options* parent_instance = nullptr,
          std::string section_name = "");

  /// Options must be explicitly copied
  ///
  ///     Option option2 = option1.copy();
  ///
  [[deprecated("Please use a reference or .copy() instead")]] Options(
      const Options& other);

  /// Copy assignment must be explicit
  ///
  ///     Option option2 = option1.copy();
  ///
  /// Note that the value can be copied using:
  ///
  ///     option2.value = option1.value;
  ///
  [[deprecated("Please use a reference or .copy() instead")]] Options&
  operator=(const Options& other); // Use a reference or .copy() method

  /// Make a deep copy of this Options,
  /// recursively copying children.
  Options copy() const;

  Options(Options&& other) noexcept;
  Options& operator=(Options&& other) noexcept;

  ~Options() = default;

  /// Get a reference to the only root instance
  static Options& root();

  /// Free all memory
  static void cleanup();

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

    AttributeType() = default;
    AttributeType(const AttributeType& other) = default;
    AttributeType(AttributeType&& other) = default;
    AttributeType& operator=(const AttributeType& other) = default;
    AttributeType& operator=(AttributeType&& other) = default;
    ~AttributeType() = default;

    /// Assignment operator, including move assignment
    using Base::operator=;

    /// Assignment from const char*
    AttributeType& operator=(const char* str) {
      operator=(std::string(str));
      return *this;
    }

    /// Initialise with a value
    /// This enables AttributeTypes to be constructed using initializer lists
    template <typename T>
    AttributeType(T value) {
      operator=(value);
    }

    /// Cast operator, which allows this class to be
    /// assigned to type T
    /// This will throw std::bad_cast if it can't be done
    template <typename T>
    operator T() const {
      return as<T>();
    }

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

  /// Return true if this value has attribute \p key
  bool hasAttribute(const std::string& key) const {
    return attributes.find(key) != attributes.end();
  }

  /// Set attributes, overwriting any already set
  ///
  /// Parameters
  /// ----------
  /// An initializer_list so that multiple attributes can be set at the same time
  ///
  /// Returns
  /// -------
  /// A reference to `this`, for method chaining
  ///
  /// Example
  /// -------
  ///
  ///     Options options;
  ///     options["value"].setAttributes({
  ///         {"units", "m/s"},
  ///         {"conversion", 10.2},
  ///         {"long_name", "some velocity"}
  ///       });
  Options& setAttributes(
      const std::initializer_list<std::pair<std::string, Options::AttributeType>>&
          attrs) {
    for (const auto& attr : attrs) {
      attributes[attr.first] = attr.second;
    }
    return *this;
  }

  /// Get a sub-section or value
  ///
  /// Example:
  ///
  /// Options parent;
  /// auto child  = parent["child"];
  ///
  /// parent is now a section.
  Options& operator[](const std::string& name);
  // Prevent ambiguous overload resolution due to operator T() below
  Options& operator[](const char* name) { return (*this)[std::string(name)]; }

  /// Get a sub-section or value
  /// If this object is not a section, or if name
  /// is not a child, then a BoutException will be thrown
  const Options& operator[](const std::string& name) const;
  const Options& operator[](const char* name) const { return (*this)[std::string(name)]; }

  /// Return type of `Options::fuzzyFind`
  struct FuzzyMatch {
    /// The matching option
    const Options& match;
    /// Edit distance from original search term
    std::string::size_type distance;
    /// Comparison operator so this works in a std::multiset
    friend bool operator<(const FuzzyMatch& lhs, const FuzzyMatch& rhs) {
      return lhs.distance < rhs.distance;
    }
  };

  /// Find approximate matches for \p name throughout the whole
  /// tree. \p distance controls the similarity of results
  ///
  /// Returns a set of possible matches ordered by similarity to \p
  /// name. A \p distance of 1 means: a single insertion, deletion,
  /// substitution, or transposition; that the case differs, for
  /// example, "key" and "KEY" match with distance 1; or that an
  /// unqualified name matches a fully-qualified name, for example
  /// "key" matches "section:key" with distance 1. Note that
  /// "first:second:key" will not (closely) match "third:fourth:key",
  /// but "key" will match both.
  std::multiset<FuzzyMatch> fuzzyFind(const std::string& name,
                                      std::string::size_type distance = 4) const;

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

  /// Assign a value to the option.
  /// This will throw an exception if already has a value
  ///
  /// Returns
  /// -------
  /// A reference to `this`, for method chaining
  ///
  /// Example
  /// -------
  ///
  /// Options option;
  /// option["test"].assign(42, "some source");
  ///
  /// Note: Specialised versions for types stored in ValueType
  template <typename T>
  Options& assign(T val, std::string source = "") {
    std::stringstream as_str;
    as_str << val;
    _set(as_str.str(), std::move(source), false);
    return *this;
  }

  /// Force to a value
  /// Overwrites any existing setting
  ///
  /// Returns
  /// -------
  /// A reference to `this`, for method chaining
  template <typename T>
  Options& force(T val, const std::string source = "") {
    is_section = true; // Invalidates any existing setting
    return assign(val, source);
  }

  /// Assign a value that is expected to vary in time.
  ///
  /// Overwrites any existing setting, and ensures the "time_dimension"
  /// attribute is set. If \p save_repeat is false, doesn't set
  /// "time_dimension". This can be useful in some generic functions
  template <typename T>
  Options& assignRepeat(T val, std::string time_dimension = "t", bool save_repeat = true,
                        std::string source = "") {
    force(val, std::move(source));
    if (save_repeat) {
      attributes["time_dimension"] = std::move(time_dimension);
    }
    return *this;
  }

  /// Test if a key is set by the user.
  /// Values set via default values are ignored.
  bool isSet() const;

  // Getting options

  /// Cast operator, which allows this class to be assigned to type
  /// T. This is only allowed for types that are members of the
  /// `ValueType` variant. For other types, please use
  /// `Options::as<T>()`
  ///
  /// Example:
  ///
  ///     Options option;
  ///     option["test"] = 2.0;
  ///     int value = option["test"];
  ///
  template <typename T, typename = typename std::enable_if_t<
                            bout::utils::isVariantMember<T, ValueType>::value>>
  operator T() const {
    return as<T>();
  }

  /// Get the value as a specified type. If there is no value then an
  /// exception is thrown. Note there are specialised versions of
  /// this template for some types.
  ///
  /// Example:
  ///
  ///     Options option;
  ///     option["test"] = 2.0;
  ///     int value = option["test"].as<int>();
  ///
  /// An optional argument is an object which the result should be similar to.
  /// The main use for this is in Field2D and Field3D specialisations,
  /// where the Mesh and cell location are taken from this input.
  ///
  /// Attributes override properties of the \p similar_to argument.
  template <typename T>
  T as(const T& UNUSED(similar_to) = {}) const {
    if (is_section) {
      throw BoutException("Option {:s} has no value", full_name);
    }

    T val;

    // Try casting. This will throw std::bad_cast if it can't be done
    try {
      val = bout::utils::variantStaticCastOrThrow<ValueType, T>(value);
    } catch (const std::bad_cast& e) {
      // If the variant is a string then we may be able to parse it

      if (bout::utils::holds_alternative<std::string>(value)) {
        std::stringstream as_str(bout::utils::get<std::string>(value));
        as_str >> val;

        // Check if the parse failed
        if (as_str.fail()) {
          throw BoutException("Option {:s} could not be parsed ('{:s}')", full_name,
                              bout::utils::variantToString(value));
        }

        // Check if there are characters remaining
        std::string remainder;
        std::getline(as_str, remainder);
        for (const unsigned char chr : remainder) {
          if (!std::isspace(chr)) {
            // Meaningful character not parsed
            throw BoutException("Option {:s} could not be parsed", full_name);
          }
        }
      } else {
        // Another type which can't be casted
        throw BoutException("Option {:s} could not be converted to type {:s}", full_name,
                            typeid(T).name());
      }
    }

    // Mark this option as used
    value_used = true; // Note this is mutable

    output_info << "\tOption " << full_name << " = " << val;
    if (attributes.count("source")) {
      // Specify the source of the setting
      output_info << " (" << bout::utils::variantToString(attributes.at("source")) << ")";
    }
    output_info << '\n';

    return val;
  }

  /// Get the value of this option. If not found,
  /// set to the default value
  template <typename T>
  T withDefault(T def) {

    // Set the type
    attributes["type"] = bout::utils::typeName<T>();

    if (is_section) {
      // Option not found
      assign(def, DEFAULT_SOURCE);
      value_used = true; // Mark the option as used

      output_info << _("\tOption ") << full_name << " = " << def << " (" << DEFAULT_SOURCE
                  << ")\n";
      return def;
    }
    T val = as<T>(def);
    // Check if this was previously set as a default option
    if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
      // Check that the default values are the same
      if (!similar(val, def)) {
        throw BoutException("Inconsistent default values for '{:s}': '{:s}' then '{:s}'",
                            full_name, bout::utils::variantToString(value),
                            toString(def));
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
    ASSERT0(def.isValue());

    if (is_section) {
      // Option not found. Copy the value from the default.
      this->_set_no_check(def.value, DEFAULT_SOURCE);

      output_info << _("\tOption ") << full_name << " = " << def.full_name << " ("
                  << DEFAULT_SOURCE << ")\n";
    } else {
      // Check if this was previously set as a default option
      if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
        // Check that the default values are the same
        if (!similar(bout::utils::variantToString(value),
                     bout::utils::variantToString(def.value))) {
          throw BoutException(
              "Inconsistent default values for '{:s}': '{:s}' then '{:s}'", full_name,
              bout::utils::variantToString(value),
              bout::utils::variantToString(def.value));
        }
      }
    }
    return *this;
  }

  /// Get the value of this option. If not found,
  /// return the default value but do not set
  template <typename T>
  T withDefault(T def) const {
    if (is_section) {
      // Option not found
      output_info << _("\tOption ") << full_name << " = " << def << " (" << DEFAULT_SOURCE
                  << ")\n";
      return def;
    }
    T val = as<T>(def);
    // Check if this was previously set as a default option
    if (bout::utils::variantEqualTo(attributes.at("source"), DEFAULT_SOURCE)) {
      // Check that the default values are the same
      if (!similar(val, def)) {
        throw BoutException("Inconsistent default values for '{:s}': '{:s}' then '{:s}'",
                            full_name, bout::utils::variantToString(value),
                            toString(def));
      }
    }
    return val;
  }

  /// Allow the user to override defaults set later, also used by the
  /// BOUT_OVERRIDE_DEFAULT_OPTION.
  template <typename T>
  T overrideDefault(T def) {

    // Set the type
    attributes["type"] = bout::utils::typeName<T>();

    if (is_section) {
      // Option not found
      assign(def, "user_default");
      is_section = false; // Prevent this default being replaced by setDefault()
      return def;
    }

    return as<T>();
  }

  /// Overloaded version for const char*
  /// Note: Different from template since return type is different to input
  std::string overrideDefault(const char* def) {
    return overrideDefault<std::string>(std::string(def));
  }

  /// Get the parent Options object
  Options& parent() {
    if (parent_instance == nullptr) {
      throw BoutException("Option {:s} has no parent", full_name);
    }
    return *parent_instance;
  }

  /// Equality operator
  /// Converts to the same type and compares
  /// This conversion may fail, throwing std::bad_cast
  template <typename T>
  bool operator==(const T& other) const {
    return as<T>() == other;
  }

  /// Overloaded equality operator for literal strings
  bool operator==(const char* other) const;

  /// Comparison operator
  template <typename T>
  bool operator<(const T& other) const {
    return as<T>() < other;
  }

  /// Overloaded comparison operator for literal strings
  bool operator<(const char* other) const;

  //////////////////////////////////////
  // Backward-compatible interface

  /// Get a pointer to the only root instance (singleton)
  static Options* getRoot() { return &root(); }

  template <typename T>
  void set(const std::string& key, T val, const std::string& source = "",
           bool force = false) {
    if (force) {
      (*this)[key].force(val, source);
    } else {
      (*this)[key].assign(val, source);
    }
  }

  // Setting options
  template <typename T>
  void forceSet(const std::string& key, T val, const std::string& source = "") {
    (*this)[key].force(val, source);
  }

  /*!
   * Test if a key is set by the user.
   * Values set via default values are ignored.
   */
  bool isSet(const std::string& key) const {
    // Note using operator[] here would result in exception if key does not exist
    if (!is_section) {
      return false;
    }
    auto child = children.find(key);
    if (child == children.end()) {
      return false;
    }
    return child->second.isSet();
  }

  /// Get options, passing in a reference to a variable
  /// which will be set to the requested value.
  template <typename T, typename U>
  void get(const std::string& key, T& val, U def, bool UNUSED(log) = false) {
    val = (*this)[key].withDefault<T>(def);
  }
  template <typename T, typename U>
  void get(const std::string& key, T& val, U def, bool UNUSED(log) = false) const {
    val = (*this)[key].withDefault<T>(def);
  }

  /// Creates new section if doesn't exist
  Options* getSection(const std::string& name) { return &(*this)[name]; }
  const Options* getSection(const std::string& name) const { return &(*this)[name]; }
  Options* getParent() const { return parent_instance; }

  //////////////////////////////////////

  /*!
   * Print string representation of this object and all sections in a tree structure
   */
  std::string str() const { return full_name; }

  /// Print just the name of this object without parent sections
  std::string name() const {
    auto pos = full_name.rfind(':');
    if (pos == std::string::npos) {
      // No parent section or sections
      return full_name;
    }
    return full_name.substr(pos + 1);
  }

  /// Return a new Options instance which contains all the values
  /// _not_ used from this instance. If an option has a "source"
  /// attribute in \p exclude_sources it is counted as having been
  /// used and so won't be included in the returned value. By default,
  /// this is "Output" and "user_default" (used by overrideDefault).
  Options getUnused(const std::vector<std::string>& exclude_sources = {
                        "Output", "user_default"}) const;

  /// Print the options which haven't been used
  void printUnused() const;

  /// Set the attribute "conditionally used" to be true for \p options
  /// and all its children/sections, causing `Options::getUnused` to
  /// assume those options have been used. This is useful to ignore
  /// options when checking for typos etc.
  void setConditionallyUsed();

  /// clean the cache of parsed options
  static void cleanCache();

  /// Read-only access to internal options and sections
  /// to allow iteration over the tree
  std::map<std::string, const Options*> subsections() const;

  const std::map<std::string, Options>& getChildren() const { return children; }

  /// Return a vector of all the full names of all the keys below this
  /// in the tree (not gauranteed to be sorted)
  std::vector<std::string> getFlattenedKeys() const;

  /// Return true if this is a value
  bool isValue() const { return not is_section; }
  /// Return true if this is a section
  bool isSection(const std::string& name = "") const;

  /// If the option value has been used anywhere
  bool valueUsed() const { return value_used; }

  /// Set a documentation string as an attribute "doc"
  /// Returns a reference to this, to allow chaining
  Options& doc(const std::string& docstring) {
    attributes["doc"] = docstring;
    return *this;
  }

  friend bool operator==(const Options& lhs, const Options& rhs) {
    if (lhs.isValue() and rhs.isValue()) {
      return lhs.value == rhs.value;
    }
    return lhs.children == rhs.children;
  }

  static std::string getDefaultSource();

  /// API for delayed loading of data from the grid file
  /// Currently only for 3D data
  using lazyLoadFunction = std::unique_ptr<std::function<Tensor<BoutReal>(
      int xstart, int xend, int ystart, int yend, int zstart, int zend)>>;
  void setLazyLoad(lazyLoadFunction func) { lazyLoad = std::move(func); }
  /// Load and get a chunk of the data
  Tensor<BoutReal> doLazyLoad(int xstart, int xend, int ystart, int yend, int zstart,
                              int zend) const {
    ASSERT1(lazyLoad != nullptr);
    return (*lazyLoad)(xstart, xend, ystart, yend, zstart, zend);
  }
  /// Some backends support to only read the data when needed.  This
  /// allows to check whether the data is loaded, or whether it needs
  /// to be loaded by doLazyLoad.
  bool is_loaded() const { return lazyLoad == nullptr; }
  /// Get the shape of the value
  std::vector<int> getShape() const;
  void setLazyShape(std::vector<int> shape) { lazy_shape = std::move(shape); }

private:
  /// The source label given to default values
  static const std::string DEFAULT_SOURCE;

  Options* parent_instance{nullptr};
  std::string full_name; // full path name for logging only

  // An Option object can either be a section or a value, defaulting to a section
  bool is_section = true;                  ///< Is this Options object a section?
  std::map<std::string, Options> children; ///< If a section then has children
  mutable bool value_used = false;         ///< Record whether this value is used

  // Function to load data
  lazyLoadFunction lazyLoad{nullptr};
  // Shape of underlying data
  std::vector<int> lazy_shape;

  template <typename T>
  void _set_no_check(T val, std::string source) {
    if (not children.empty()) {
      throw BoutException(
          "Trying to assign value to Option '{}', but it's a non-empty section",
          full_name);
    }

    value = std::move(val);
    attributes["source"] = std::move(source);
    value_used = false;
    is_section = false;
  }

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
              _("Options: Setting a value from same source ({:s}) to new value "
                "'{:s}' - old value was '{:s}'."),
              source, toString(val), bout::utils::variantToString(value));
        }
      }
    }

    _set_no_check(std::move(val), std::move(source));
  }

  /// Tests if two values are similar.
  template <typename T>
  bool similar(T lhs, T rhs) const {
    return lhs == rhs;
  }
};

// Specialised assign methods for types stored in ValueType
template <>
inline Options& Options::assign<>(bool val, std::string source) {
  _set(val, std::move(source), false);
  return *this;
}
template <>
inline Options& Options::assign<>(int val, std::string source) {
  _set(val, std::move(source), false);
  return *this;
}
template <>
inline Options& Options::assign<>(BoutReal val, std::string source) {
  _set(val, std::move(source), false);
  return *this;
}
template <>
inline Options& Options::assign<>(std::string val, std::string source) {
  _set(std::move(val), std::move(source), false);
  return *this;
}
// Note: const char* version needed to avoid conversion to bool
template <>
inline Options& Options::assign<>(const char* val, std::string source) {
  _set(std::string(val), std::move(source), false);
  return *this;
}
// Note: Field assignments don't check for previous assignment (always force)
template <>
Options& Options::assign<>(Field2D val, std::string source);
template <>
Options& Options::assign<>(Field3D val, std::string source);
template <>
Options& Options::assign<>(FieldPerp val, std::string source);
template <>
Options& Options::assign<>(Array<BoutReal> val, std::string source);
template <>
Options& Options::assign<>(Matrix<BoutReal> val, std::string source);
template <>
Options& Options::assign<>(Tensor<BoutReal> val, std::string source);

/// Specialised similar comparison methods
template <>
inline bool Options::similar<BoutReal>(BoutReal lhs, BoutReal rhs) const {
  return fabs(lhs - rhs) < 1e-10;
}

/// Specialised as routines
template <>
std::string Options::as<std::string>(const std::string& similar_to) const;
template <>
int Options::as<int>(const int& similar_to) const;
template <>
BoutReal Options::as<BoutReal>(const BoutReal& similar_to) const;
template <>
bool Options::as<bool>(const bool& similar_to) const;
template <>
Field2D Options::as<Field2D>(const Field2D& similar_to) const;
template <>
Field3D Options::as<Field3D>(const Field3D& similar_to) const;
template <>
FieldPerp Options::as<FieldPerp>(const FieldPerp& similar_to) const;
template <>
Array<BoutReal> Options::as<Array<BoutReal>>(const Array<BoutReal>& similar_to) const;
template <>
Matrix<BoutReal> Options::as<Matrix<BoutReal>>(const Matrix<BoutReal>& similar_to) const;
template <>
Tensor<BoutReal> Options::as<Tensor<BoutReal>>(const Tensor<BoutReal>& similar_to) const;

/// Convert \p value to string
std::string toString(const Options& value);

/// Output a stringified \p value to a stream
///
/// This is templated to avoid implict casting: anything is
/// convertible to an `Options`, and we want _exactly_ `Options`
template <class T, typename = bout::utils::EnableIfOptions<T>>
inline std::ostream& operator<<(std::ostream& out, const T& value) {
  return out << toString(value);
}

namespace bout {
/// Check if the global Options contains any unused keys and throw an
/// exception if so. This check can be skipped by setting
/// `input:error_on_unused_options=false` in the global Options.
void checkForUnusedOptions();
/// Check if the given \p options contains any unused keys and throw
/// an exception if so.
///
/// The error message contains helpful suggestions on possible
/// misspellings, and how to automatically fix most common errors with
/// library options. The \p data_dir and \p option_file arguments are
/// used to customise the error message for the actual input file used
void checkForUnusedOptions(const Options& options, const std::string& data_dir,
                           const std::string& option_file);
} // namespace bout

namespace bout {
namespace details {
/// Implementation of fmt::formatter<Options> in a non-template class
/// so that we can put the function definitions in the .cxx file,
/// avoiding lengthy recompilation if we change it
struct OptionsFormatterBase {
  auto parse(fmt::format_parse_context& ctx) -> fmt::format_parse_context::iterator;
  auto format(const Options& options, fmt::format_context& ctx) const
      -> fmt::format_context::iterator;

private:
  /// Include the 'doc' attribute, if present
  bool docstrings{false};
  /// If an option is unused add a comment and whether it is
  /// conditionally unused
  bool unused{false};
  /// If true, print variables as 'section:variable', rather than a
  /// section header '[section]' and plain 'variable'
  bool inline_section_names{false};
  /// Only include the key name, and not the value
  bool key_only{false};
  /// Include the 'source' attribute, if present
  bool source{false};
  /// Format string to passed down to subsections
  std::string format_string;
};
} // namespace details
} // namespace bout

/// Format `Options` to string. Format string specification is:
///
/// - 'd': include 'doc' attribute if present
/// - 'i': inline section names
/// - 'k': only print the key, not the value
/// - 's': include 'source' attribute if present
template <>
struct fmt::formatter<Options> : public bout::details::OptionsFormatterBase {};

// NOLINTBEGIN(cppcoreguidelines-macro-usage)

/// Define for reading options which passes the variable name
#define OPTION(options, var, def) pointer(options)->get(#var, var, def)

#define OPTION2(options, var1, var2, def)    \
  {                                          \
    pointer(options)->get(#var1, var1, def); \
    pointer(options)->get(#var2, var2, def); \
  }

#define OPTION3(options, var1, var2, var3, def) \
  {                                             \
    pointer(options)->get(#var1, var1, def);    \
    pointer(options)->get(#var2, var2, def);    \
    pointer(options)->get(#var3, var3, def);    \
  }

#define OPTION4(options, var1, var2, var3, var4, def) \
  {                                                   \
    pointer(options)->get(#var1, var1, def);          \
    pointer(options)->get(#var2, var2, def);          \
    pointer(options)->get(#var3, var3, def);          \
    pointer(options)->get(#var4, var4, def);          \
  }

#define OPTION5(options, var1, var2, var3, var4, var5, def) \
  {                                                         \
    pointer(options)->get(#var1, var1, def);                \
    pointer(options)->get(#var2, var2, def);                \
    pointer(options)->get(#var3, var3, def);                \
    pointer(options)->get(#var4, var4, def);                \
    pointer(options)->get(#var5, var5, def);                \
  }

#define OPTION6(options, var1, var2, var3, var4, var5, var6, def) \
  {                                                               \
    pointer(options)->get(#var1, var1, def);                      \
    pointer(options)->get(#var2, var2, def);                      \
    pointer(options)->get(#var3, var3, def);                      \
    pointer(options)->get(#var4, var4, def);                      \
    pointer(options)->get(#var5, var5, def);                      \
    pointer(options)->get(#var6, var6, def);                      \
  }

#define VAROPTION(options, var, def)                              \
  {                                                               \
    if (pointer(options)->isSet(#var)) {                          \
      pointer(options)->get(#var, var, def);                      \
    } else {                                                      \
      Options::getRoot()->getSection("all")->get(#var, var, def); \
    }                                                             \
  }

/// Define for over-riding library defaults for options, should be called in global
/// namespace so that the new default is set before main() is called.
#define BOUT_OVERRIDE_DEFAULT_OPTION(name, value)                                  \
  namespace {                                                                      \
  const auto BOUT_CONCAT(user_default,                                             \
                         __LINE__) = Options::root()[name].overrideDefault(value); \
  }

// NOLINTEND(cppcoreguidelines-macro-usage)

#endif // OPTIONS_H
