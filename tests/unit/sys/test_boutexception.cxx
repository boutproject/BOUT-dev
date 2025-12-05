#include "test_extras.hxx"
#include "bout/boutexception.hxx"
#include "gtest/gtest.h"

#include <fmt/ranges.h>

#include <string>
#include <vector>

// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
struct BoutExceptionTest : public ::testing::Test {
  BoutExceptionTest() { BoutException::enableBacktrace(); }
  ~BoutExceptionTest() override { BoutException::disableBacktrace(); }
};

TEST_F(BoutExceptionTest, ThrowCorrect) {
  EXPECT_THROW(throw BoutException("test"), BoutException);
}

TEST_F(BoutExceptionTest, What) {
  std::string test_message{"Test message"};
  try {
    throw BoutException(test_message);
  } catch (const BoutException& e) {
    EXPECT_TRUE(IsSubString(e.what(), test_message));
  }
  try {
    throw BoutException("this is {}", "second");
  } catch (const BoutException& e) {
    std::string message(e.what());
    EXPECT_TRUE(IsSubString(message, "this is second"));
  }
}

TEST_F(BoutExceptionTest, GetBacktrace) {
  std::string test_message{"Test message"};
  try {
    throw BoutException(test_message);
  } catch (const BoutException& e) {
    std::string expected_1{"#2"};
    std::string expected_2{"bout_test_main"};
    // Should be able to find something about backtrace
    EXPECT_TRUE(IsSubString(e.what(), expected_1));
    EXPECT_TRUE(IsSubString(e.what(), expected_2));
  }
}

TEST_F(BoutExceptionTest, FmtJoin) {
  const std::vector things = {1, 2, 3, 4};
  const std::string expected = "list: 1, 2, 3, 4";
  const BoutException exception{"list: {}", fmt::join(things, ", ")};
  EXPECT_TRUE(IsSubString(std::string{exception.what()}, expected));
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
