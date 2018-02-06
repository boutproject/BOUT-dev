#include "gtest/gtest.h"
#include "test_extras.hxx"

#include <cmath>

::testing::AssertionResult IsSubString(const std::string &str,
                                       const std::string &substring) {
  if (str.find(substring) != std::string::npos) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure() << '"' << substring << "\" not found in " << str;
  }
}

::testing::AssertionResult IsField3DEqualBoutReal(const Field3D &field, BoutReal number,
                                                  BoutReal tolerance) {
  for (const auto &i : field) {
    if (fabs(field[i] - number) > tolerance) {
      return ::testing::AssertionFailure() << "Field3D(" << i.x << ", " << i.y << ", "
                                           << i.z << ") == " << field[i]
                                           << "; Expected: " << number;
    }
  }

  return ::testing::AssertionSuccess();
}

::testing::AssertionResult IsField2DEqualBoutReal(const Field2D &field, BoutReal number,
                                                  BoutReal tolerance) {
  for (const auto &i : field) {
    if (fabs(field[i] - number) > tolerance) {
      return ::testing::AssertionFailure() << "Field2D(" << i.x << ", " << i.y
                                           << ") == " << field[i]
                                           << "; Expected: " << number;
    }
  }

  return ::testing::AssertionSuccess();
}
