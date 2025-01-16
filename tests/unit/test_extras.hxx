#ifndef TEST_EXTRAS_H__
#define TEST_EXTRAS_H__

#include "gtest/gtest.h"

#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "bout/bout_types.hxx"
#include "bout/field.hxx"
#include "bout/region.hxx"

class Field2D;
class Field3D;

static constexpr BoutReal BoutRealTolerance{1e-15};
// FFTs have a slightly looser tolerance than other functions
static constexpr BoutReal FFTTolerance{1.e-12};

/// Does \p str contain \p substring?
::testing::AssertionResult IsSubString(const std::string& str,
                                       const std::string& substring);

void fillField(Field3D& f, std::vector<std::vector<std::vector<BoutReal>>> values);
void fillField(Field2D& f, std::vector<std::vector<BoutReal>> values);

using bout::utils::EnableIfField;

/// Returns a field filled with the result of \p fill_function at each point
/// Arbitrary arguments can be passed to the field constructor
template <class T, class... Args, typename = EnableIfField<T>>
T makeField(const std::function<BoutReal(typename T::ind_type&)>& fill_function,
            Args&&... args) {
  T result{std::forward<Args>(args)...};
  result.allocate();

  for (auto i : result) {
    result[i] = fill_function(i);
  }

  return result;
}

/// Teach googletest how to print SpecificInds
template <IND_TYPE N>
inline std::ostream& operator<<(std::ostream& out, const SpecificInd<N>& index) {
  return out << index.ind;
}

/// Helpers to get the type of a Field as a string
auto inline getFieldType([[maybe_unused]] const Field2D& field) -> std::string {
  return "Field2D";
}
auto inline getFieldType([[maybe_unused]] const Field3D& field) -> std::string {
  return "Field3D";
}
auto inline getFieldType([[maybe_unused]] const FieldPerp& field) -> std::string {
  return "FieldPerp";
}

/// Helpers to get the (x, y, z) index values, along with the
/// single-index of a Field index
auto inline getIndexXYZ(const Ind2D& index) -> std::string {
  std::stringstream ss;
  ss << index.x() << ", " << index.y() << "; [" << index.ind << "]";
  return ss.str();
}
auto inline getIndexXYZ(const Ind3D& index) -> std::string {
  std::stringstream ss;
  ss << index.x() << ", " << index.y() << ", " << index.z() << "; [" << index.ind << "]";
  return ss.str();
}
auto inline getIndexXYZ(const IndPerp& index) -> std::string {
  std::stringstream ss;
  ss << index.x() << ", " << index.y() << ", " << index.z() << "; [" << index.ind << "]";
  return ss.str();
}

/// Is \p field equal to \p reference, with a tolerance of \p tolerance?
template <class T, class U, typename = EnableIfField<T, U>>
auto IsFieldEqual(const T& field, const U& reference,
                  const std::string& region = "RGN_ALL",
                  BoutReal tolerance = BoutRealTolerance) -> ::testing::AssertionResult {
  for (auto i : field.getRegion(region)) {
    if (fabs(field[i] - reference[i]) > tolerance) {
      return ::testing::AssertionFailure()
             << getFieldType(field) << "(" << getIndexXYZ(i) << ") == " << field[i]
             << "; Expected: " << reference[i];
    }
  }
  return ::testing::AssertionSuccess();
}

/// Is \p field equal to \p reference, with a tolerance of \p tolerance?
/// Overload for BoutReals
template <class T, typename = EnableIfField<T>>
auto IsFieldEqual(const T& field, BoutReal reference,
                  const std::string& region = "RGN_ALL",
                  BoutReal tolerance = BoutRealTolerance) -> ::testing::AssertionResult {
  for (auto i : field.getRegion(region)) {
    if (fabs(field[i] - reference) > tolerance) {
      return ::testing::AssertionFailure()
             << getFieldType(field) << "(" << getIndexXYZ(i) << ") == " << field[i]
             << "; Expected: " << reference;
    }
  }
  return ::testing::AssertionSuccess();
}

#endif //  TEST_EXTRAS_H__
