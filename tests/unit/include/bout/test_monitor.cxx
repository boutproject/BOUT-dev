#include "gtest/gtest.h"

#include "bout/monitor.hxx"
#include "boutexception.hxx"

TEST(MonitorTest, IsMultiple) {
  EXPECT_TRUE(isMultiple(1., 4.));
  EXPECT_TRUE(isMultiple(2., 4.));
  EXPECT_TRUE(isMultiple(4., 2.));
  // Floating point accuracy is strange - so have plenty of tests
  auto lst = {1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e+0,
              1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8};
  for (auto a : lst) {
    for (auto b : lst) {
      EXPECT_TRUE(isMultiple(a, b));
    }
  }

  EXPECT_FALSE(isMultiple(3., 4.));
  EXPECT_FALSE(isMultiple(4., 3.));
  EXPECT_FALSE(isMultiple(2., 5.));
  EXPECT_FALSE(isMultiple(5., 2.));
  EXPECT_FALSE(isMultiple(10., 1e8 + 1));
}

#if CHECK > 1
TEST(MonitorTest, IsMultipleConditions) {
  EXPECT_THROW(isMultiple(0., 4.), BoutException);
  EXPECT_THROW(isMultiple(4., 0.), BoutException);
  EXPECT_THROW(isMultiple(-1., 4.), BoutException);
  EXPECT_THROW(isMultiple(4., -1.), BoutException);
}
#endif
