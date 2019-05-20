#include "gtest/gtest.h"
#include "boutexception.hxx"

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
    throw BoutException("%s", "second");
  } catch (const BoutException &e) {
    std::string message(e.what());
    EXPECT_EQ(message, "second");
  }
}


TEST(BoutExceptionTest, GetBacktrace) {
  std::string test_message{"Test message"};
  try {
    throw BoutException(test_message);
  } catch (const BoutException &e) {
    std::string expected{"[bt] #1 ./serial_tests"};
#ifdef BACKTRACE
    // Should be able to find something about backtrace
    EXPECT_NE(e.getBacktrace().find(expected), std::string::npos);
#else
    // Should *not* be able to find something about backtrace
    EXPECT_EQ(e.getBacktrace().find(expected), std::string::npos);
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
