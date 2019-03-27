
#pragma once

#ifndef TYPE_NAME_HXX
#define TYPE_NAME_HXX

#include <typeinfo>
#include <string>
#include "bout_types.hxx"

class Field2D;
class Field3D;

namespace bout {
namespace utils {
  
  template <typename T>
  std::string typeName() {
    return typeid(T).name();
  }

  template <>
  std::string typeName<bool>();

  template <>
  std::string typeName<int>();

  template <>
  std::string typeName<std::string>();
  
  // Specialised for BOUT++ types to ensure that the result is human-readable
  template <>
  std::string typeName<BoutReal>();
  
  template <>
  std::string typeName<Field2D>();

  template <>
  std::string typeName<Field3D>();
}
}

#endif //TYPE_NAME_HXX
