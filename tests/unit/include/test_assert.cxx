#include "gtest/gtest.h"

#include <bout/assert.hxx>
#include <boutexception.hxx>

TEST(AssertTest, Assert0) {
  EXPECT_NO_THROW(ASSERT0(true));
#if CHECK >= 0
  EXPECT_THROW(ASSERT0(false), BoutException);
#else
  EXPECT_NO_THROW(ASSERT0(false));
#endif
}

TEST(AssertTest, Assert1) {
  EXPECT_NO_THROW(ASSERT1(true));
#if CHECK >= 1
  EXPECT_THROW(ASSERT1(false), BoutException);
#else
  EXPECT_NO_THROW(ASSERT1(false));
#endif
}

TEST(AssertTest, Assert2) {
  EXPECT_NO_THROW(ASSERT2(true));
#if CHECK >= 2
  EXPECT_THROW(ASSERT2(false), BoutException);
#else
  EXPECT_NO_THROW(ASSERT2(false));
#endif
}

TEST(AssertTest, Assert3) {
  EXPECT_NO_THROW(ASSERT3(true));
#if CHECK >= 3
  EXPECT_THROW(ASSERT3(false), BoutException);
#else
  EXPECT_NO_THROW(ASSERT3(false));
#endif
}
