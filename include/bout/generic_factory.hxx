/// \file Base type for factories

#pragma once
#ifndef __BOUT_GENERIC_FACTORY_H__
#define __BOUT_GENERIC_FACTORY_H__

#include "boutexception.hxx"
#include "options.hxx"

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

/// Generic Factory, adapted from Modern C++ Design/Loki by A. Alexandrescu
///
/// Use with RegisterInFactory to provide a generic way of creating
/// new derived types at runtime. By default assumes the type can be
/// created with an Options*
///
/// Uses static polymorphism (via CRTP) to overload static data. This
/// is done by inheriting from this class and templating on the
/// _inherited_ class, and then providing four public data members:
///
/// Example:
///
///     class Base {
///     public:
///       Base(Options*) {}
///     };
///     class Derived : public Base {
///     public:
///       Derived(Options*) : Base({}) {}
///     };
///
///     class MyFactory : public Factory<Base, MyFactory> {
///       static constexpr auto type_name = "Base";
///       static constexpr auto section_name = "base";
///       static constexpr auto option_name = "type";
///       static constexpr auto default_type = "derived_type";
///     };
///
///     RegisterInFactory<Base, Derived, MyFactory> register("derived_type");
///     auto foo = MyFactory::getInstance().create("derived_type");
///
/// @tparam BaseType       The base class that this factory creates
/// @tparam DerivedFactory The derived factory inheriting from this class
/// @tparam TypeCreator    The function signature for creating a new BaseType
///
/// MIT Licence
template <class BaseType, class DerivedFactory,
          class TypeCreator = std::function<std::unique_ptr<BaseType>(Options*)>>
class Factory {
protected:
  // Type returned from the creation function
  using ReturnType = typename TypeCreator::result_type;

  Factory() = default;

  /// Disgusting hack to defeat linker throwing out the registration
  /// symbols. If necessary, override this and put the (empty)
  /// implementation in the same TU as the registration symbols
  static void ensureRegistered() {}

  /// Storage of the creation functions
  std::map<std::string, TypeCreator> type_map;

  /// Return either \p options or the section from root
  Options* optionsOrDefaultSection(Options* options) const {
    if (options == nullptr) {
      options = &Options::root()[DerivedFactory::section_name];
    }
    return options;
  }

public:
  virtual ~Factory() = default;

  /// Get the singleton instance
  static DerivedFactory& getInstance() {
    static DerivedFactory instance{};
    DerivedFactory::ensureRegistered();
    return instance;
  }

  /// Return the name of the default type to create
  static constexpr auto getDefaultType() { return DerivedFactory::default_type; }
  /// Return the name of the section to get from the root Options
  static constexpr auto getSectionName() { return DerivedFactory::section_name; }
  /// Return the name of the Option value that sets the type
  static constexpr auto getOptionName() { return DerivedFactory::option_name; }

  /// Add a new type \p name to the factory
  ///
  /// @param[in] name     An identifier for this type
  /// @param[in] creator  A function for creating this type
  /// @returns true if the type was successfully added
  virtual bool add(const std::string& name, TypeCreator creator) {
    return type_map.insert(std::make_pair(name, creator)).second;
  }

  /// Remove a type \p name from the factory
  ///
  /// @param[in] name  The identifier for the type to be removed
  /// @returns true if the type was successfully removed
  virtual bool remove(const std::string& name) { return type_map.erase(name) == 1; }

  /// Get the name of the type to create
  ///
  /// Looks in \p options first, then the root Options
  ///
  /// @param[in] options  Options section to look for type name
  /// @returns the name of the type to create
  std::string getType(Options* options = nullptr) const {
    options = optionsOrDefaultSection(options);
    return (*options)[DerivedFactory::option_name].withDefault(
        DerivedFactory::getDefaultType());
  }

  /// Create a new object using the type set in \p options
  ReturnType create(Options* options = nullptr) const {
    options = optionsOrDefaultSection(options);
    return create(getType(options), options);
  }

  /// Create a new object of type \p name using the root Options
  ReturnType create(const std::string& name) const {
    return create(name, &Options::root()[DerivedFactory::section_name]);
  }

  /// Create a new object of type \p name
  ///
  /// @param[in] name  The identifier for the type to be created
  /// @returns the new object
  template <typename... Args>
  ReturnType create(const std::string& name, Args&&... args) const {
    auto index = type_map.find(name);
    if (index != std::end(type_map)) {
      return index->second(std::forward<Args>(args)...);
    }
    // List available options in error
    std::string available;
    auto available_list = listAvailable();
    for (auto i : available_list) {
      available += i + "\n";
    }
    throw BoutException("Error when trying to create a {0:s}: Could not find "
                        "'{1:s}'\nAvailable {0:s}s are:\n{2:s}",
                        DerivedFactory::type_name, name, available);
  }

  /// Create a new object of the type given in options["type"]
  ///
  /// @param[in] options  The Options object to get the type to be created from
  /// @returns the new object
  template <typename... Args>
  ReturnType create(Options* options, Args&&... args) const {
    return create(getType(optionsOrDefaultSection(options)), args...);
  }

  /// List available types that can be created
  ///
  /// @returns a vector of std::string
  std::vector<std::string> listAvailable() const {
    std::vector<std::string> available;
    available.reserve(type_map.size());
    for (const auto& name : type_map) {
      available.push_back(name.first);
    }
    return available;
  }
};

/// Helper class for adding new types to Factory
///
/// See Factory for example
///
/// Adapted from
/// http://www.drdobbs.com/conversations-abstract-factory-template/184403786
///
/// @tparam BaseType       Which factory to add \p DerivedType to
/// @tparam DerivedType    The new type to add to Factory<BaseType>
template <class BaseType, class DerivedType, class DerivedFactory>
class RegisterInFactory {
public:
  RegisterInFactory(const std::string& name) {
    DerivedFactory::getInstance().add(name,
                                      [](Options* options) -> std::unique_ptr<BaseType> {
                                        return std::make_unique<DerivedType>(options);
                                      });
  }
};

#endif // __BOUT_GENERIC_FACTORY_H__
