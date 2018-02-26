#include "gtest/gtest.h"

#include "bout/array.hxx"
#include "boutexception.hxx"

#include <iostream>

// In order to keep these tests independent, they need to use
// different sized arrays in order to not just reuse the data from
// previous tests

class ArrayTest : public ::testing::Test {
public:
  ArrayTest() { Array<double>::useStore(true); }
  // Note: Calling cleanup() disables the store
  ~ArrayTest() { }
};

TEST_F(ArrayTest, ArraySize) {
  Array<double> a(5);

  ASSERT_FALSE(a.empty());
  EXPECT_EQ(a.size(), 5);
  EXPECT_TRUE(a.unique());
}

TEST_F(ArrayTest, ArrayValues) {
  Array<double> a(10);

  int count = 0;
  for (auto &i : a) {
    i = count++;
  }

  EXPECT_DOUBLE_EQ(a[1], 1);
  EXPECT_DOUBLE_EQ(a[9], 9);
}

TEST_F(ArrayTest, CopyArray) {
  Array<double> a(15);

  int count = 0;
  for (auto &i : a) {
    i = count++;
  }

  Array<double> b(a);

  ASSERT_FALSE(b.empty());
  EXPECT_EQ(b.size(), 15);
  EXPECT_DOUBLE_EQ(b[5], 5);
  EXPECT_FALSE(a.unique());
  EXPECT_FALSE(b.unique());
}

TEST_F(ArrayTest, MakeUnique) {
  Array<double> a(20);

  int count = 0;
  for (auto &i : a) {
    i = count++;
  }

  Array<double> b(a);
  // Make both b and a unique
  b.ensureUnique();

  ASSERT_TRUE(b.unique());
  ASSERT_TRUE(a.unique());
  EXPECT_EQ(a.size(), 20);
  EXPECT_EQ(b.size(), 20);

  // Should have the same values
  for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ++ai, ++bi) {
    EXPECT_DOUBLE_EQ(*ai, *bi);
  }
}

TEST_F(ArrayTest, ReleaseData) {
  Array<double> a(25);
  Array<double> b(a);
  // Make both b and a unique
  b.ensureUnique();

  // Release the data. Should put into store.
  a.clear();
  EXPECT_TRUE(a.empty());
  EXPECT_FALSE(b.empty());
  EXPECT_EQ(a.size(), 0);
}

TEST_F(ArrayTest, RetrieveData) {
  Array<double> a(30);

  int count = 0;
  for (auto &i : a) {
    i = count++;
  }

  Array<double> b(a);
  // Make both b and a unique
  b.ensureUnique();

  // Release the data. Should put into store.
  a.clear();

  // Construct, retrieve from store, and move assign
  a = Array<double>(30);

  EXPECT_FALSE(a.empty());
  EXPECT_EQ(a.size(), 30);
  EXPECT_EQ(a[4], 4); // Test if reused data from store
  EXPECT_TRUE(a.unique());

  a.ensureUnique(); // Should have no effect

  EXPECT_EQ(a.size(), 30);
  EXPECT_EQ(a[4], 4); // Test if reused data from store
  EXPECT_TRUE(a.unique());
}

TEST_F(ArrayTest, Assignment) {
  Array<double> a(35);
  Array<double> b(35);
  // Assign
  a = b;

  EXPECT_FALSE(a.unique());
  EXPECT_FALSE(b.unique());
}

#if CHECK > 2
TEST_F(ArrayTest, OutOfBoundsThrow) {
  Array<double> a(34);
  EXPECT_NO_THROW(a[33] = 1.0);
  EXPECT_THROW(a[34] = 1.0, BoutException);
  EXPECT_THROW(a[-1] = 1.0, BoutException);
}
#endif
