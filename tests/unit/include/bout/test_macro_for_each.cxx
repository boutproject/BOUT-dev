#include "gtest/gtest.h"

#include "bout/macro_for_each.hxx"

#include <sstream>

#define INC(x) {++(x);}

// Test that the macro is expanded
// the correct number of times
TEST(MacroForEachTest, ExpandCount) {
  int a{0};

  MACRO_FOR_EACH(INC, a);
  EXPECT_EQ(a, 1);

  a = 0;
  MACRO_FOR_EACH(INC, a, a);
  EXPECT_EQ(a, 2);

  // At least 10 arguments supported
  a = 0;
  MACRO_FOR_EACH(INC, a, a, a, a, a, a, a, a, a, a);
  EXPECT_EQ(a, 10);
}

#define PUSH_NAME(var) ss << #var;

// Check that the macro applies in the order given
TEST(MacroForEachTest, Order) {
  std::stringstream ss;

  MACRO_FOR_EACH(PUSH_NAME, a, b, c);
  EXPECT_EQ(ss.str(), "abc");
}

// No braces are put around the expansion
// This is needed in some cases, but can
// lead to unexpected behaviour.
TEST(MacroForEachTest, NoBraces) {
  int a = 0;

  // Only the first expansion is disabled
  if(false)
    MACRO_FOR_EACH(INC, a, a, a);

  EXPECT_EQ(a, 2);
}

void inc(int &val) {
  val++;
}

TEST(MacroForEachFnTest, IncludeSemiColon) {
  int a = 0;

  // Should expand as
  // inc(a); inc(a); inc(a);
  MACRO_FOR_EACH_FN(inc, a, a, a);

  EXPECT_EQ(a, 3);
}

TEST(MacroForEachFnTest, ExpandCount) {
  int a{0};

  MACRO_FOR_EACH_FN(inc, a);
  EXPECT_EQ(a, 1);

  a = 0;
  MACRO_FOR_EACH_FN(inc, a, a);
  EXPECT_EQ(a, 2);

  // At least 10 arguments supported
  a = 0;
  MACRO_FOR_EACH_FN(inc, a, a, a, a, a, a, a, a, a, a);
  EXPECT_EQ(a, 10);
}

TEST(MacroForEachFnTest, Order) {
  std::stringstream ss;

  MACRO_FOR_EACH_FN(PUSH_NAME, a, b, c);
  EXPECT_EQ(ss.str(), "abc");
}

