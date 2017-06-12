#include "gtest/gtest.h"
#include "bout/sys/range.hxx"

#include <vector>

TEST(RangeTest, SimpleRangeStartEnd) {
  RangeIterator range(0, 4);

  EXPECT_EQ(range.min(), 0);
  EXPECT_EQ(range.max(), 4);
}

TEST(RangeTest, SimpleRangeFirst) {
  RangeIterator range(0, 4);

  range.first();
  EXPECT_EQ(*range, 0);
}

TEST(RangeTest, SimpleRangeAdvance) {
  RangeIterator range(0, 4);

  range.next();
  EXPECT_EQ(*range, 1);

  ++range;
  EXPECT_EQ(*range, 2);

  range++;
  EXPECT_EQ(*range, 3);
}

TEST(RangeTest, SimpleRangeDone) {
  RangeIterator range(0, 1);

  range.first();
  EXPECT_FALSE(range.isDone());

  range++;
  range++;
  EXPECT_TRUE(range.isDone());
}

TEST(RangeTest, SimpleRangeEnd) {
  RangeIterator range(0, 1);

  range.first();
  EXPECT_FALSE(range.isDone());

  range = range.end();
  EXPECT_TRUE(range.isDone());
}

TEST(RangeTest, SimpleRangeLoop) {
  int count = 0;
  RangeIterator range(0, 4);

  for (range.first(); !range.isDone(); ++range, ++count) {
    EXPECT_EQ(*range, count);
    ASSERT_TRUE(*range >= 0);
    ASSERT_TRUE(*range <= 4);
  }
}

TEST(RangeTest, SimpleRangeOptionalFirst) {
  int count = 0;

  for (RangeIterator range(0, 4); !range.isDone(); ++range, ++count) {
    EXPECT_EQ(*range, count);
    ASSERT_TRUE(*range >= 0);
    ASSERT_TRUE(*range <= 4);
  }
}

TEST(RangeTest, SimpleRangeNonZeroStart) {
  int count = 5;
  RangeIterator range(5, 9);

  EXPECT_EQ(range.min(), 5);
  EXPECT_EQ(range.max(), 9);

  for (range.first(); !range.isDone(); ++range, ++count) {
    EXPECT_EQ(*range, count);
    ASSERT_TRUE(*range >= 5);
    ASSERT_TRUE(*range <= 9);
  }
}

TEST(RangeTest, JoinRange) {
  int count = 0;
  RangeIterator range_b(5, 10);
  RangeIterator range_a(0, 4, range_b);

  for (range_a.first(); !range_a.isDone(); ++range_a, ++count) {
    EXPECT_EQ(*range_a, count);
    ASSERT_TRUE(*range_a <= 10);
  }
}

TEST(RangeTest, JoinRangeOperator) {
  int count = 0;
  RangeIterator range_a(0, 4);
  RangeIterator range_b(5, 10);
  range_a += range_b;

  for (range_a.first(); !range_a.isDone(); ++range_a, ++count) {
    EXPECT_EQ(*range_a, count);
    ASSERT_TRUE(*range_a <= 10);
  }
}

TEST(RangeTest, RemoveRangeOperator) {
  std::vector<int> full_range = {0, 1, 2, 3, 4};
  int count = 0;
  RangeIterator range_a(0, 10);
  RangeIterator range_b(5, 10);
  range_a -= range_b;

  for (range_a.first(); !range_a.isDone(); ++range_a, ++count) {
    EXPECT_EQ(*range_a, full_range[count]);
    ASSERT_TRUE(*range_a <= 4);
  }
}

TEST(RangeTest, JoinRangeWithGap) {
  std::vector<int> full_range = {0, 1, 2, 3, 4, 10, 11, 12, 13, 14};
  int count = 0;
  RangeIterator range_b(10, 14);
  RangeIterator range_a(0, 4, range_b);

  for (range_a.first(); !range_a.isDone(); ++range_a, ++count) {
    EXPECT_EQ(*range_a, full_range[count]);
  }
}

TEST(RangeTest, RangeIntersectInt) {
  RangeIterator range(0, 10);

  EXPECT_TRUE(range.intersects(4));
  EXPECT_FALSE(range.intersects(12));
  EXPECT_FALSE(range.intersects(-4));
}

TEST(RangeTest, RangeIntersectIntAll) {
  RangeIterator range_a(10, 20);
  RangeIterator range_b(0, 10, range_a);

  EXPECT_TRUE(range_a.intersects(12, true));
  EXPECT_FALSE(range_b.intersects(22, true));
}

TEST(RangeTest, RangesIntersect) {
  RangeIterator range_a(5, 10);
  RangeIterator range_b(0, 10);
  RangeIterator range_c(15, 20);

  EXPECT_TRUE(range_a.intersects(range_b));
  EXPECT_TRUE(range_b.intersects(range_a));
  EXPECT_FALSE(range_a.intersects(range_c));
  EXPECT_FALSE(range_c.intersects(range_a));
}

TEST(RangeTest, RangesIntersectAll) {
  RangeIterator range_a(5, 7);
  RangeIterator range_b(10, 20);
  RangeIterator range_c(0, 10, range_b);

  EXPECT_TRUE(range_a.intersects(range_c, true));
  EXPECT_TRUE(range_c.intersects(range_a, true));
}

TEST(RangeTest, Equality) {
  RangeIterator range_a(0, 1);
  RangeIterator range_b(0, 1);

  EXPECT_TRUE(range_a == range_a);
  EXPECT_FALSE(range_a != range_a);
  EXPECT_TRUE(range_a != range_b);
  EXPECT_FALSE(range_a == range_b);
}
