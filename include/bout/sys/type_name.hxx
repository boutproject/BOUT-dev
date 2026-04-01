
#pragma once

#ifndef TYPE_NAME_HXX
#define TYPE_NAME_HXX

#include "bout/array.hxx"
#include "bout/bout_types.hxx"

#include <string>
#include <typeinfo>

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

template <>
std::string typeName<FieldPerp>();

template <>
std::string typeName<Array<int>>();
template <>
std::string typeName<Array<BoutReal>>();
template <>
std::string typeName<Matrix<int>>();
template <>
std::string typeName<Matrix<BoutReal>>();
template <>
std::string typeName<Tensor<int>>();
template <>
std::string typeName<Tensor<BoutReal>>();
} // namespace utils
} // namespace bout

#endif //TYPE_NAME_HXX
