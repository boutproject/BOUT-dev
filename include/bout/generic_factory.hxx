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
template <typename BaseType,
          typename TypeCreator = std::function<std::unique_ptr<BaseType>()>>
class Factory {
public:
  virtual ~Factory() = default;

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

  /// Create a new object of type \p name
  ///
  /// @param[in] name  The identifier for the type to be created
  /// @returns a pointer to the new object
  template <typename... Args>
  std::unique_ptr<BaseType> create(const std::string& name, Args&&... args) {
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
    throw BoutException("Available:\n%s\nCould not find '%s'", available.c_str(),
                        name.c_str());
  }

  /// List available types that can be created
  ///
  /// @returns a vector of std::string
  virtual std::vector<std::string> listAvailable() const {
    std::vector<std::string> available;
    for (const auto& name : type_map) {
      available.push_back(name.first);
    }
    return available;
  }

  /// Get the singleton instance
  static Factory& getInstance() {
    static Factory instance;
    return instance;
  }

protected:
  std::map<std::string, TypeCreator> type_map;
  Factory() = default;
};

template <class BaseType, class DerivedFactory,
          class TypeCreator = std::function<std::unique_ptr<BaseType>(Options*)>,
          class BaseFactory = Factory<BaseType, TypeCreator>>
class StandardFactory : public BaseFactory {
protected:
  using ReturnType = typename TypeCreator::result_type;
  using BaseFactoryType = BaseFactory;
  StandardFactory() = default;
  static void ensureRegistered() {}

public:
  static DerivedFactory& getInstance() {
    static DerivedFactory instance{};
    DerivedFactory::ensureRegistered();
    return instance;
  }

  static constexpr auto getDefaultType() { return DerivedFactory::default_type; }
  static constexpr auto getSectionName() { return DerivedFactory::section_name; }
  static constexpr auto getOptionName() { return DerivedFactory::option_name; }

  std::string getType(Options* options = nullptr) {
    if (options == nullptr) {
      options = &Options::root()[DerivedFactory::section_name];
    }

    return (*options)[DerivedFactory::option_name].withDefault(
        DerivedFactory::getDefaultType());
  }

  ReturnType create(Options* options = nullptr) {
    auto type = getType(options);
    return create(type, options);
  }

  ReturnType create(const std::string& name) {
    return create(name, &Options::root()[DerivedFactory::section_name]);
  }

  template <typename... Args>
  ReturnType create(const std::string& name, Args&&... args) {
    try {
      return static_cast<BaseFactory*>(this)->create(name, std::forward<Args>(args)...);
    } catch (const BoutException& e) {
      throw BoutException("Error when trying to create a %s: %s",
                          DerivedFactory::type_name, e.what());
    }
  }
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
  RegisterInFactory(const std::string& name) {
    Factory<BaseType>::getInstance().add(name, []() -> std::unique_ptr<BaseType> {
      return std::make_unique<DerivedType>();
    });
  }
};

template <class BaseType, class DerivedType, class DerivedFactory>
class RegisterInStandardFactory {
public:
  RegisterInStandardFactory(const std::string& name) {
    DerivedFactory::getInstance().add(
      name, [](Options* options) -> std::unique_ptr<BaseType> {
        return std::make_unique<DerivedType>(options);
      });
  }
};

#endif // __BOUT_GENERIC_FACTORY_H__
