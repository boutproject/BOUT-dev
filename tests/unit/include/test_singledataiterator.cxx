#include "gtest/gtest.h"

#include "bout/singledataiterator.hxx"
#include "unused.hxx"

TEST(SIndexRange, Begin) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();
  EXPECT_EQ(iter->i, 0);
}

// Dereferencing an end() iterator is an error, so we need to test
// end() works a little indirectly. If the addition and less-than
// tests fail, this one is suspect!
TEST(SIndexRange, End) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin() + 8;
  auto iter_end = range.end();
  EXPECT_EQ(iter->i, 8);
  // iter_end is one-past the last element of region
  EXPECT_TRUE(iter < iter_end);
}

TEST(SIndexRange, PrefixIncrement) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();
  ++iter;

  EXPECT_EQ(iter->i, 1);
}

TEST(SIndexRange, PostfixIncrement) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();
  auto iter2 = iter++;

  EXPECT_EQ(iter->i, 1);
  EXPECT_EQ(iter2->i, 0);
}

TEST(SIndexRange, PrefixDecrement) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.end();
  --iter;

  EXPECT_EQ(iter->i, 8);
}

TEST(SIndexRange, PostfixDecrement) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  // end() is one-past-the-last element, so we need to decrement it in
  // order to be able to dereference it
  auto iter = --(range.end());
  auto iter2 = iter--;

  EXPECT_EQ(iter->i, 7);
  EXPECT_EQ(iter2->i, 8);
}

TEST(SIndexRange, NotEquals) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter != iter2);
  ++iter;
  EXPECT_FALSE(iter != iter2);
}

TEST(SIndexRange, Equals) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;
  ++iter;

  EXPECT_TRUE(iter == iter2);
  ++iter;
  EXPECT_FALSE(iter == iter2);
}

TEST(SIndexRange, LessThan) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter < iter2);
  ++iter;
  EXPECT_FALSE(iter < iter2);
  ++iter;
  EXPECT_FALSE(iter < iter2);
}

TEST(SIndexRange, MoreThan) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter2 > iter);
  ++iter;
  EXPECT_FALSE(iter2 > iter);
  ++iter;
  EXPECT_FALSE(iter2 > iter);
}

TEST(SIndexRange, LessThanOrEqualTo) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter <= iter2);
  ++iter;
  EXPECT_TRUE(iter <= iter2);
  ++iter;
  EXPECT_FALSE(iter <= iter2);
}

TEST(SIndexRange, MoreThanOrEqualTo) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter2 >= iter);
  ++iter;
  EXPECT_TRUE(iter2 >= iter);
  ++iter;
  EXPECT_FALSE(iter2 >= iter);
}

TEST(SIndexRange, PlusEqualsInt) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();

  iter += 2;
  EXPECT_EQ(iter->i, 2);
}

TEST(SIndexRange, PlusInt) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();

  auto iter2 = iter + 3;
  EXPECT_EQ(iter2->i, 3);
}

TEST(SIndexRange, IntPlus) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();

  auto iter2 = 5 + iter;
  EXPECT_EQ(iter2->i, 5);
}

TEST(SIndexRange, MinusEqualsInt) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.end();

  iter -= 2;
  EXPECT_EQ(iter->i, 7);
}

TEST(SIndexRange, MinusInt) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.end();

  auto iter2 = iter - 3;
  EXPECT_EQ(iter2->i, 6);
}

TEST(SIndexRange, MinusIterator) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto start = range.begin();
  auto end = range.end();

  EXPECT_EQ(end - start, region.size());
}

TEST(SIndexRange, IndexInt) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};

  auto iter = range.begin();

  EXPECT_EQ(iter[4].i, 4);
}

TEST(SIndexRange, Iteration) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};
  RegionIndices region2;

  int count = 0;
  for (auto iter2 = range.begin(); iter2 != range.end(); ++iter2) {
    ++count;
    region2.push_back(region[iter2->i]);
  }

  EXPECT_EQ(count, region.size());
  EXPECT_EQ(region2, region);
}

TEST(SIndexRange, RangeBasedForLoop) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SIndexRange range{2, 2, 2, region};
  RegionIndices region2;

  int count = 0;
  for (const auto &iter : range) {
    ++count;
    region2.push_back(region[iter.i]);
  }

  EXPECT_EQ(count, region.size());
  EXPECT_EQ(region2, region);
}
