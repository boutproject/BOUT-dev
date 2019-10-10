/// \file Base type for factories

#pragma once
#ifndef __BOUT_GENERIC_FACTORY_H__
#define __BOUT_GENERIC_FACTORY_H__

#include "boutexception.hxx"

#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

/// Generic Factory, adapted from Modern C++ Design/Loki by A. Alexandrescu
///
/// Use with RegisterInFactory to provide a generic way of creating
/// new derived types at runtime
///
/// Example:
///
///     class Base {};
///     class Derived : public Base {};
///     RegisterInFactory<Base, Derived> register("derived_type");
///     auto foo = Factory<Base>::getInstance().create("derived_type");
///
/// TODO: Use std::unique_ptr<BaseType> instead of BaseType*
///
/// @tparam BaseType       The base class that this factory creates
/// @tparam TypeCreator    The function signature for creating a new BaseType
///
/// MIT Licence
template <typename BaseType, typename TypeCreator = std::function<BaseType *()>>
class Factory {
public:
  /// Add a new type \p name to the factory
  ///
  /// @param[in] name     An identifier for this type
  /// @param[in] creator  A function for creating this type
  /// @returns true if the type was successfully added
  virtual bool add(const std::string &name, TypeCreator creator) {
    return type_map.insert(std::make_pair(name, creator)).second;
  }

  /// Remove a type \p name from the factory
  ///
  /// @param[in] name  The identifier for the type to be removed
  /// @returns true if the type was successfully removed
  virtual bool remove(const std::string &name) { return type_map.erase(name) == 1; }

  /// Create a new object of type \p name
  ///
  /// @param[in] name  The identifier for the type to be created
  /// @returns a pointer to the new object
  template<typename... Args>
  BaseType *create(const std::string &name, Args&&... args) {
    auto index = type_map.find(name);
    if (index != std::end(type_map)) {
      return index->second(std::forward<Args>(args) ...);
    }
    // List available options in error
    std::string available;
    auto available_list = listAvailable();
    for (auto i : available_list) {
      available += i + "\n";
    }
    throw BoutException("Available:\n%s\nCould not find '%s'", available.c_str(), name.c_str());
  }

  /// List available types that can be created
  ///
  /// @returns a vector of std::string
  virtual std::vector<std::string> listAvailable() {
    std::vector<std::string> available;
    for (const auto &name : type_map) {
      available.push_back(name.first);
    }
    return available;
  }

  /// Get the singleton instance
  static Factory &getInstance() {
    static Factory instance;
    return instance;
  }

protected:
  std::map<std::string, TypeCreator> type_map;
  Factory() = default;
};

/// Helper class for adding new types to Factory
///
/// Example:
///
///     class Base {};
///     class Derived : public Base {};
///     RegisterInFactory<Base, Derived> register("derived_type");
///     auto foo = Factory<Base>::getInstance().create("derived_type");
///
/// Adapted from
/// http://www.drdobbs.com/conversations-abstract-factory-template/184403786
///
/// @tparam BaseType       Which factory to add \p DerivedType to
/// @tparam DerivedType    The new type to add to Factory<BaseType>
template <typename BaseType, typename DerivedType>
class RegisterInFactory {
public:
  RegisterInFactory(const std::string &name) {
    Factory<BaseType>::getInstance().add(name,
                                         []() -> BaseType * { return new DerivedType; });
  }
};

#endif // __BOUT_GENERIC_FACTORY_H__
