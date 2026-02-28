#include "gtest/gtest.h"

#include "bout/arraynd.hxx"

class ArrayNDTest : public ::testing::Test {
public:
  ArrayNDTest() {}
  ~ArrayNDTest() {}
};

TEST_F(ArrayNDTest, ArrayNDSize1D) {
  ArrayND<double, 1> a(5);

  ASSERT_FALSE(a.empty());
  EXPECT_EQ(a.size(), 5);
  EXPECT_EQ(std::tuple_size<decltype(a.shape())>::value, 1);
  EXPECT_EQ(std::get<0>(a.shape()), 5);
  EXPECT_TRUE(a.unique());
}

TEST_F(ArrayNDTest, ArrayNDSize2D) {
  ArrayND<double, 2> a(5, 10);

  ASSERT_FALSE(a.empty());
  EXPECT_EQ(a.size(), 5 * 10);
  EXPECT_EQ(std::tuple_size<decltype(a.shape())>::value, 2);
  EXPECT_EQ(std::get<0>(a.shape()), 5);
  EXPECT_EQ(std::get<1>(a.shape()), 10);
  EXPECT_TRUE(a.unique());
}

TEST_F(ArrayNDTest, ArrayNDUnique2D) {
  ArrayND<double, 2> a(5, 10);
  EXPECT_TRUE(a.unique());
  auto b = a[0];
  EXPECT_FALSE(a.unique());
  EXPECT_FALSE(b.unique());
  b.ensureUnique();
  EXPECT_TRUE(b.unique());
  EXPECT_TRUE(a.unique());
  auto c = a.data[0];
  EXPECT_FALSE(c.unique());
  EXPECT_FALSE(a.unique());
  a.ensureUnique();
  EXPECT_TRUE(c.unique());
  EXPECT_TRUE(a.unique());
}
