
#pragma once

#ifndef TYPE_NAME_HXX
#define TYPE_NAME_HXX

#include "bout/array.hxx"
#include "bout/bout_types.hxx"

#include <string>
#include <type_traits>
#include <typeinfo>

// Forward declarations
class Field2D;
class Field3D;
class FieldPerp;
template <class T>
class Matrix;
template <class T>
class Tensor;

namespace bout::utils {

template <typename T>
std::string typeName() {
  if constexpr (std::is_same_v<T, bool>) {
    return "bool";
  }
  if constexpr (std::is_same_v<T, int>) {
    return "int";
  }
  if constexpr (std::is_same_v<T, std::string>) {
    return "string";
  }
  // Specialised for BOUT++ types to ensure that the result is human-readable
  if constexpr (std::is_same_v<T, BoutReal>) {
    return "BoutReal";
  }
  if constexpr (std::is_same_v<T, Field2D>) {
    return "Field2D";
  }
  if constexpr (std::is_same_v<T, Field3D>) {
    return "Field3D";
  }
  if constexpr (std::is_same_v<T, FieldPerp>) {
    return "FieldPerp";
  }
  if constexpr (std::is_same_v<T, Array<int>>) {
    return "Array<int>";
  }
  if constexpr (std::is_same_v<T, Array<BoutReal>>) {
    return "Array<BoutReal>";
  }
  if constexpr (std::is_same_v<T, Matrix<int>>) {
    return "Matrix<int>";
  }
  if constexpr (std::is_same_v<T, Matrix<BoutReal>>) {
    return "Matrix<BoutReal>";
  }
  if constexpr (std::is_same_v<T, Tensor<int>>) {
    return "Tensor<int>";
  }
  if constexpr (std::is_same_v<T, Tensor<BoutReal>>) {
    return "Tensor<BoutReal>";
  }
  return typeid(T).name();
}

} // namespace bout::utils

#endif // TYPE_NAME_HXX
