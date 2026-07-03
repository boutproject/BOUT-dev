
#pragma once

#ifndef TYPE_NAME_HXX
#define TYPE_NAME_HXX

#include "bout/array.hxx"
#include "bout/bout_types.hxx"

#include <string>
#include <type_traits> // Required for std::is_same_v
#include <typeinfo>

// Forward declarations
class Field2D;
class Field3D;
class FieldPerp;
template <class T>
class Matrix;
template <class T>
class Tensor;

namespace bout {
namespace utils {

template <typename T>
std::string typeName() {
  if constexpr (std::is_same_v<T, bool>)
    return "bool";
  else if constexpr (std::is_same_v<T, int>)
    return "int";
  else if constexpr (std::is_same_v<T, std::string>)
    return "string";
  // Specialised for BOUT++ types to ensure that the result is human-readable
  else if constexpr (std::is_same_v<T, BoutReal>)
    return "BoutReal";
  else if constexpr (std::is_same_v<T, Field2D>)
    return "Field2D";
  else if constexpr (std::is_same_v<T, Field3D>)
    return "Field3D";
  else if constexpr (std::is_same_v<T, FieldPerp>)
    return "FieldPerp";
  else if constexpr (std::is_same_v<T, Array<int>>)
    return "Array<int>";
  else if constexpr (std::is_same_v<T, Array<BoutReal>>)
    return "Array<BoutReal>";
  else if constexpr (std::is_same_v<T, Matrix<int>>)
    return "Matrix<int>";
  else if constexpr (std::is_same_v<T, Matrix<BoutReal>>)
    return "Matrix<BoutReal>";
  else if constexpr (std::is_same_v<T, Tensor<int>>)
    return "Tensor<int>";
  else if constexpr (std::is_same_v<T, Tensor<BoutReal>>)
    return "Tensor<BoutReal>";
  else {
    return typeid(T).name();
  }
}

} // namespace utils
} // namespace bout

#endif // TYPE_NAME_HXX
