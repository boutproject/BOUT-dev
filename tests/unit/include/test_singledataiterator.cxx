#include "gtest/gtest.h"

#include "bout/singledataiterator.hxx"
#include "unused.hxx"

TEST(SingleDataIterator, Begin) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.begin();
  EXPECT_EQ(iter.region_iter, region.begin());
  EXPECT_EQ(*(iter.region_iter), 0);
}

TEST(SingleDataIterator, End) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.end();
  EXPECT_EQ(iter.region_iter, region.end());
}

TEST(SingleDataIterator, PrefixIncrement) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.begin();
  ++iter;

  EXPECT_EQ(*(iter.region_iter), 1);
}

TEST(SingleDataIterator, PostfixIncrement) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.begin();
  iter++;

  EXPECT_EQ(*(iter.region_iter), 1);
}

TEST(SingleDataIterator, PrefixDecrement) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.end();
  --iter;

  EXPECT_EQ(*(iter.region_iter), 8);
}

TEST(SingleDataIterator, PostfixDecrement) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.end();
  iter--;

  EXPECT_EQ(*(iter.region_iter), 8);
}

TEST(SingleDataIterator, NotEquals) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.begin();
  auto foo = iter;
  ++foo;

  EXPECT_NE(iter, foo);
}

TEST(SingleDataIterator, Equals) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.begin();
  auto foo = iter;
  ++foo;
  ++iter;

  EXPECT_EQ(iter, foo);
}

TEST(SingleDataIterator, LessThan) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.begin();
  auto foo = iter;
  ++foo;

  EXPECT_TRUE(iter < foo);
  ++iter;
  EXPECT_FALSE(iter < foo);
  ++iter;
  EXPECT_FALSE(iter < foo);
}

TEST(SingleDataIterator, MoreThan) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.begin();
  auto foo = iter;
  ++foo;

  EXPECT_TRUE(foo > iter);
  ++iter;
  EXPECT_FALSE(foo > iter);
  ++iter;
  EXPECT_FALSE(foo > iter);
}

TEST(SingleDataIterator, LessThanOrEqualTo) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.begin();
  auto foo = iter;
  ++foo;

  EXPECT_TRUE(iter <= foo);
  ++iter;
  EXPECT_TRUE(iter <= foo);
  ++iter;
  EXPECT_FALSE(iter <= foo);
}

TEST(SingleDataIterator, MoreThanOrEqualTo) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  iter.begin();
  auto foo = iter;
  ++foo;

  EXPECT_TRUE(foo >= iter);
  ++iter;
  EXPECT_TRUE(foo >= iter);
  ++iter;
  EXPECT_FALSE(foo >= iter);
}

TEST(SingleDataIterator, Iteration) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  int count = 0;
  for (auto foo = iter.begin(); foo != iter.end(); ++foo) {
    ++count;
  }

  EXPECT_EQ(count, region.size());
}

TEST(SingleDataIterator, RangeBasedForLoop) {
  RegionIndices region{0, 1, 2, 3, 4, 5, 6, 7, 8};
  SingleDataIterator iter{2, 2, 2, region};

  int count = 0;
  for (auto &UNUSED(foo) : iter) {
    ++count;
  }

  EXPECT_EQ(count, region.size());
}
