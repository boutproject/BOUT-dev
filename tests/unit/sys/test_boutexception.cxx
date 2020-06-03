#include "bout/build_config.hxx"

#include "gtest/gtest.h"
#include "boutexception.hxx"
#include "test_extras.hxx"

#include <iostream>
#include <string>

TEST(BoutExceptionTest, ThrowCorrect) {
  EXPECT_THROW(throw BoutException("test"), BoutException);
}

TEST(BoutExceptionTest, What) {
  std::string test_message{"Test message"};
  try {
    throw BoutException(test_message);
  } catch (const BoutException &e) {
    EXPECT_EQ(e.what(), test_message);
  }
  try {
    throw BoutException("this is {}", "second");
  } catch (const BoutException &e) {
    std::string message(e.what());
    EXPECT_EQ(message, "this is second");
  }
}


TEST(BoutExceptionTest, GetBacktrace) {
  std::string test_message{"Test message"};
  try {
    throw BoutException(test_message);
  } catch (const BoutException &e) {
    std::string expected_1{"[bt] #1"};
    std::string expected_2{"serial_tests"};
#if BOUT_USE_BACKTRACE
    // Should be able to find something about backtrace
    EXPECT_TRUE(IsSubString(e.getBacktrace(), expected_1));
    EXPECT_TRUE(IsSubString(e.getBacktrace(), expected_2));
#else
    // Should *not* be able to find something about backtrace
    EXPECT_FALSE(IsSubString(e.getBacktrace(), expected_1));
    EXPECT_FALSE(IsSubString(e.getBacktrace(), expected_2));
#endif
  }
}


TEST(BoutRhsFailTest, ThrowCorrect) {
  EXPECT_THROW(throw BoutRhsFail("RHS Fail test"), BoutRhsFail);
}

TEST(BoutParallelRhsFailTest, ThrowCorrect) {
  EXPECT_THROW(BoutParallelThrowRhsFail(1, "RHS Fail test"), BoutRhsFail);
}

TEST(BoutParallelRhsFailTest, ThrowNoError) {
  EXPECT_NO_THROW(BoutParallelThrowRhsFail(0, "RHS Fail test"));
}

TEST(BoutIterationFailTest, ThrowCorrect) {
  EXPECT_THROW(throw BoutIterationFail("Iteration fail test"), BoutIterationFail);
}
