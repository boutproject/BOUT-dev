#include "gtest/gtest.h"
#include "test_extras.hxx"

::testing::AssertionResult IsSubString(const std::string &str,
                                       const std::string &substring) {
  if (str.find(substring) != std::string::npos) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure() << substring << " not found in " << str;
  }
}
