#include "gtest/gtest.h"

#include "bout/array.hxx"
#include "boutexception.hxx"

#include <iostream>
#include <numeric>

// In order to keep these tests independent, they need to use
// different sized arrays in order to not just reuse the data from
// previous tests
//


class ArrayTest : public ::testing::Test {
public:
  ArrayTest() { Array<double>::useStore(true); }
  // Note: Calling cleanup() disables the store
  virtual ~ArrayTest() = default;
};

TEST_F(ArrayTest, ArraySize) {
  Array<double> a(5);

  ASSERT_FALSE(a.empty());
  EXPECT_EQ(a.size(), 5);
  EXPECT_TRUE(a.unique());
}

TEST_F(ArrayTest, ArrayValues) {
  Array<double> a(10);

  std::iota(a.begin(), a.end(), 0);

  EXPECT_DOUBLE_EQ(a[1], 1);
  EXPECT_DOUBLE_EQ(a[9], 9);
}

TEST_F(ArrayTest, CopyArrayConstructor) {
  Array<double> a{15};

  std::iota(a.begin(), a.end(), 0);

  EXPECT_TRUE(a.unique());

  Array<double> b{a};

  ASSERT_FALSE(a.empty());
  ASSERT_FALSE(b.empty());
  EXPECT_EQ(b.size(), 15);
  EXPECT_DOUBLE_EQ(b[5], 5);
  EXPECT_FALSE(a.unique());
  EXPECT_FALSE(b.unique());
}

TEST_F(ArrayTest, CopyArrayOperator) {
  Array<double> a{15};

  std::iota(a.begin(), a.end(), 0);

  EXPECT_TRUE(a.unique());

  Array<double> b;
  b = a;

  ASSERT_FALSE(a.empty());
  ASSERT_FALSE(b.empty());
  EXPECT_EQ(b.size(), 15);
  EXPECT_DOUBLE_EQ(b[5], 5);
  EXPECT_FALSE(a.unique());
  EXPECT_FALSE(b.unique());
}

TEST_F(ArrayTest, CopyArrayNonMemberFunction) {
  Array<double> a{15};

  std::iota(a.begin(), a.end(), 0);

  Array<double> b;
  b = copy(a);

  ASSERT_FALSE(b.empty());
  EXPECT_EQ(b.size(), 15);
  EXPECT_DOUBLE_EQ(b[5], 5);
  EXPECT_TRUE(a.unique());
  EXPECT_TRUE(b.unique());
}

TEST_F(ArrayTest, SwapArray) {
  Array<double> a{15};

  std::iota(a.begin(), a.end(), 0);

  Array<double> b{10};

  std::iota(b.begin(), b.end(), 5);

  EXPECT_EQ(a.size(), 15);
  EXPECT_EQ(b.size(), 10);
  EXPECT_DOUBLE_EQ(a[5], 5);
  EXPECT_DOUBLE_EQ(b[5], 10);

  swap(a, b);

  EXPECT_EQ(a.size(), 10);
  EXPECT_EQ(b.size(), 15);
  EXPECT_DOUBLE_EQ(a[5], 10);
  EXPECT_DOUBLE_EQ(b[5], 5);
}

TEST_F(ArrayTest, MoveArrayConstructor) {
  Array<double> a{15};

  std::iota(a.begin(), a.end(), 0);

  Array<double> b{std::move(a)};

  ASSERT_TRUE(a.empty());
  ASSERT_FALSE(b.empty());
  EXPECT_EQ(b.size(), 15);
  EXPECT_DOUBLE_EQ(b[5], 5);
  EXPECT_FALSE(a.unique());
  EXPECT_TRUE(b.unique());
}

TEST_F(ArrayTest, Reallocate) {
  Array<double> a{};

  ASSERT_TRUE(a.empty());

  // Reallocate from empty
  a.reallocate(15);
  std::iota(a.begin(), a.end(), 0);

  ASSERT_FALSE(a.empty());
  EXPECT_EQ(a.size(), 15);
  EXPECT_DOUBLE_EQ(a[5], 5);
  EXPECT_TRUE(a.unique());

  // Reallocate to smaller
  a.reallocate(7);
  std::iota(a.begin(), a.end(), 10);

  ASSERT_FALSE(a.empty());
  EXPECT_EQ(a.size(), 7);
  EXPECT_DOUBLE_EQ(a[5], 15);

  // Reallocate to larger
  a.reallocate(30);
  std::iota(a.begin(), a.end(), 20);

  ASSERT_FALSE(a.empty());
  EXPECT_EQ(a.size(), 30);
  EXPECT_DOUBLE_EQ(a[5], 25);
}

TEST_F(ArrayTest, MakeUnique) {
  Array<double> a(20);

  std::iota(a.begin(), a.end(), 0);

  Array<double> b(a);

  std::iota(b.begin(), b.end(), 1);

  EXPECT_FALSE(b.unique());
  EXPECT_FALSE(a.unique());
  EXPECT_EQ(a.size(), 20);
  EXPECT_EQ(b.size(), 20);

  // Should have the same values
  for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ++ai, ++bi) {
    EXPECT_DOUBLE_EQ(*ai, *bi);
  }

  // Make both b and a unique
  b.ensureUnique();

  std::iota(b.begin(), b.end(), 2);

  EXPECT_TRUE(b.unique());
  EXPECT_TRUE(a.unique());
  EXPECT_EQ(a.size(), 20);
  EXPECT_EQ(b.size(), 20);

  // Should have the same values but offset
  for (auto ai = a.begin(), bi = b.begin(); ai != a.end(); ++ai, ++bi) {
    EXPECT_DOUBLE_EQ(*ai + 1, *bi);
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

  std::iota(a.begin(), a.end(), 0);

  Array<double> b(a);
  // Make both b and a unique
  b.ensureUnique();

  // Release the data. Should put into store.
  a.clear();

  // Construct, retrieve from store, and move assign
  //  note new array design doesn't use store
  a = Array<double>(30);

  EXPECT_FALSE(a.empty());
  EXPECT_EQ(a.size(), 30);
  if(a.useStore()) {
    EXPECT_EQ(a[4], 4); // Test if reused data from store
  }
  EXPECT_TRUE(a.unique());

  a.ensureUnique(); // Should have no effect

  EXPECT_EQ(a.size(), 30);
  if(a.useStore()) {
    EXPECT_EQ(a[4], 4); // Test if reused data from store
  }
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

#if CHECK > 2 && ! BOUT_USE_CUDA
TEST_F(ArrayTest, OutOfBoundsThrow) {
  Array<double> a(34);
  EXPECT_NO_THROW(a[33] = 1.0);
  EXPECT_THROW(a[34] = 1.0, BoutException);
  EXPECT_THROW(a[-1] = 1.0, BoutException);
}
#endif
