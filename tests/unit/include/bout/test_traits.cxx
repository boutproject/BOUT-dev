#include "gtest/gtest.h"

#include "bout/field.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldperp.hxx"
#include "bout/traits.hxx"

#include <type_traits>

TEST(BoutTraitsTest, IsField) {
  using namespace bout::utils;
  static_assert(is_Field_v<Field>, "is_Field<Field> should be true");
  static_assert(is_Field_v<Field2D>, "is_Field<Field2D> should be true");
  static_assert(is_Field_v<Field3D>, "is_Field<Field3D> should be true");
  static_assert(is_Field_v<FieldPerp>, "is_Field<FieldPerp> should be true");
}

TEST(BoutTraitsTest, IsField2D) {
  using namespace bout::utils;
  static_assert(!is_Field2D_v<Field>, "is_Field2D<Field> should be false");
  static_assert(is_Field2D_v<Field2D>, "is_Field2D<Field2D> should be true");
  static_assert(!is_Field2D_v<Field3D>, "is_Field2D<Field3D> should be false");
  static_assert(!is_Field2D_v<FieldPerp>, "is_Field2D<FieldPerp> should be false");
}

TEST(BoutTraitsTest, IsField3D) {
  using namespace bout::utils;
  static_assert(!is_Field3D_v<Field>, "is_Field3D<Field> should be false");
  static_assert(!is_Field3D_v<Field2D>, "is_Field3D<Field2D> should be false");
  static_assert(is_Field3D_v<Field3D>, "is_Field3D<Field3D> should be true");
  static_assert(!is_Field3D_v<FieldPerp>, "is_Field3D<FieldPerp> should be false");
}

TEST(BoutTraitsTest, IsFieldPerp) {
  using namespace bout::utils;
  static_assert(!is_FieldPerp_v<Field>, "is_FieldPerp<Field> should be false");
  static_assert(!is_FieldPerp_v<Field2D>, "is_FieldPerp<Field2D> should be false");
  static_assert(!is_FieldPerp_v<Field3D>, "is_FieldPerp<Field3D> should be false");
  static_assert(is_FieldPerp_v<FieldPerp>, "is_FieldPerp<FieldPerp> should be true");
}

namespace {
struct CorrectForFields {};
struct CorrectForBoutReal {};

template <class T, class U, typename = bout::utils::EnableIfField<T, U>>
auto function_overloads(T, U) -> CorrectForFields {
  return {};
}

// Commenting out the below template should break the
// EnableIfFieldOverloads test
template <class T>
auto function_overloads(T, BoutReal) -> CorrectForBoutReal {
  return {};
}

struct OnlyForField2D {};
struct OnlyForField3D {};
struct OnlyForFieldPerp {};

template <class T, class U, typename = bout::utils::EnableIfField2D<T, U>>
auto specific_function_overloads(T, U) -> OnlyForField2D {
  return {};
}

template <class T, class U, typename = bout::utils::EnableIfField3D<T, U>>
auto specific_function_overloads(T, U) -> OnlyForField3D {
  return {};
}

template <class T, class U, typename = bout::utils::EnableIfFieldPerp<T, U>>
auto specific_function_overloads(T, U) -> OnlyForFieldPerp {
  return {};
}

template <class T, class U, class ResultType = bout::utils::EnableIfField<T, U>>
auto example_function(T, U) -> ResultType {
  return {};
}

} // namespace

TEST(BoutTraitsTest, EnableIfFieldOverloads) {
  using namespace bout::utils;
  static_assert(std::is_same_v<CorrectForFields,
                               decltype(function_overloads(std::declval<Field2D>(),
                                                           std::declval<Field3D>()))>,
                "EnableIfField should enable function_overloads for two Fields");
  static_assert(
      std::is_same_v<CorrectForBoutReal,
                     decltype(function_overloads(std::declval<Field2D>(),
                                                 std::declval<BoutReal>()))>,
      "EnableIfField should disable one overload of function_overloads.\n"
      "This static_assert should fail if the `function_overloads(T, BoutReal)` template\n"
      "is commented out");
  static_assert(
      std::is_same_v<OnlyForField2D,
                     decltype(specific_function_overloads(std::declval<Field2D>(),
                                                          std::declval<Field2D>()))>,
      "EnableIfField should pick the correct version of specific_function_overloads");
  static_assert(
      std::is_same_v<OnlyForField3D,
                     decltype(specific_function_overloads(std::declval<Field3D>(),
                                                          std::declval<Field3D>()))>,
      "EnableIfField should pick the correct version of specific_function_overloads");
  static_assert(
      std::is_same_v<OnlyForFieldPerp,
                     decltype(specific_function_overloads(std::declval<FieldPerp>(),
                                                          std::declval<FieldPerp>()))>,
      "EnableIfField should pick the correct version of specific_function_overloads");
}

TEST(BoutTraitsTest, EnableIfFieldReturnType) {
  using namespace bout::utils;
  static_assert(
      std::is_same_v<Field2D, decltype(example_function(std::declval<Field2D>(),
                                                        std::declval<Field2D>()))>,
      "EnableIfField should return Field2D for two Field2Ds");
  static_assert(
      std::is_same_v<Field3D, decltype(example_function(std::declval<Field2D>(),
                                                        std::declval<Field3D>()))>,
      "EnableIfField should return Field3D for a Field2D and a Field3D");
  static_assert(
      std::is_same_v<Field3D, decltype(example_function(std::declval<Field3D>(),
                                                        std::declval<Field2D>()))>,
      "EnableIfField should return Field3D for a Field2D and a Field3D");
  static_assert(
      std::is_same_v<Field3D, decltype(example_function(std::declval<Field3D>(),
                                                        std::declval<Field3D>()))>,
      "EnableIfField should return Field3D for a Field3D and a Field3D");
}
