#include "test_extras.hxx"
#include "gtest/gtest.h"

#include <cmath>

::testing::AssertionResult IsSubString(const std::string& str,
                                       const std::string& substring) {
  if (str.find(substring) != std::string::npos) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure() << '"' << substring << "\" not found in " << str;
  }
}

void fillField(Field3D& f, std::vector<std::vector<std::vector<BoutReal>>> values) {
  f.allocate();
  Ind3D i{0};
  for (auto& x : values) {
    for (auto& y : x) {
      for (auto& z : y) {
        f[i] = z;
        ++i;
      }
    }
  }
}

void fillField(Field2D& f, std::vector<std::vector<BoutReal>> values) {
  f.allocate();
  Ind2D i{0};
  for (auto& x : values) {
    for (auto& y : x) {
      f[i] = y;
      ++i;
    }
  }
}
