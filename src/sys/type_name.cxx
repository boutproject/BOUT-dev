#include "bout/sys/type_name.hxx"

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldperp.hxx"

namespace bout {
namespace utils {

template <>
std::string typeName<bool>() {
  return "bool";
}

template <>
std::string typeName<int>() {
  return "int";
}

template <>
std::string typeName<std::string>() {
  return "string";
}

template <>
std::string typeName<BoutReal>() {
  return "BoutReal";
}

template <>
std::string typeName<Field2D>() {
  return "Field2D";
}

template <>
std::string typeName<Field3D>() {
  return "Field3D";
}

template <>
std::string typeName<FieldPerp>() {
  return "FieldPerp";
}

template <>
std::string typeName<Array<int>>() {
  return "Array<int>";
}
template <>
std::string typeName<Array<BoutReal>>() {
  return "Array<BoutReal>";
}
template <>
std::string typeName<Matrix<int>>() {
  return "Matrix<int>";
}
template <>
std::string typeName<Matrix<BoutReal>>() {
  return "Matrix<BoutReal>";
}
template <>
std::string typeName<Tensor<int>>() {
  return "Tensor<int>";
}
template <>
std::string typeName<Tensor<BoutReal>>() {
  return "Tensor<BoutReal>";
}

} // namespace utils
} // namespace bout
