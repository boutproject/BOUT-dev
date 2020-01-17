#include <type_traits>
#include <algorithm>

#include "test_extras.hxx"
#include "gtest/gtest.h"

#include <bout/region.hxx>
#include <bout/operatorstencil.hxx>

template <typename T>
class IndexOffsetStructTests : public ::testing::Test {
public:
  IndexOffsetStructTests() {
    zero = T(0, std::is_same<T, IndPerp>::value ? 1 : 5,
	     std::is_same<T, Ind2D>::value ? 1 : 7);
  }
  
  IndexOffset<T> noOffset;
  T zero;
};

template<class T>
std::ostream& operator<<(std::ostream& os, const IndexOffset<T>& offset) {
  return os << "{" << offset.dx << ", " << offset.dy << ", "
	    << offset.dz << "}";
}

using IndTypes = ::testing::Types<Ind3D, Ind2D, IndPerp>;
TYPED_TEST_SUITE(IndexOffsetStructTests, IndTypes);

TYPED_TEST(IndexOffsetStructTests, Equals) {
  IndexOffset<TypeParam> rhs1 = {0, 0, 0}, rhs2 = {1, 0, 0};
  EXPECT_TRUE(this->noOffset == rhs1);
  EXPECT_FALSE(this->noOffset == rhs2);
}

TYPED_TEST(IndexOffsetStructTests, NotEquals) {
  IndexOffset<TypeParam> rhs1 = {0, 0, 0}, rhs2 = {1, 0, 0};
  EXPECT_FALSE(this->noOffset != rhs1);
  EXPECT_TRUE(this->noOffset != rhs2);
}

TYPED_TEST(IndexOffsetStructTests, Comparison) {
  IndexOffset<TypeParam> val100 = {1, 0, 0}, val010 = {0, 1, 0}, val001 = {0, 0, 1},
                         val110 = {1, 1, 0}, val101 = {1, 0, 1}, val011 = {0, 1, 1},
                         val111 = {1, 1, 1};
  ASSERT_TRUE(this->noOffset < val100);
  ASSERT_TRUE(this->noOffset < val010);
  ASSERT_TRUE(this->noOffset < val001);
  ASSERT_TRUE(val010 < val100);
  ASSERT_TRUE(val001 < val010);
  ASSERT_TRUE(val110 < val111);
  ASSERT_TRUE(val011 < val101);
  ASSERT_FALSE(val011 < val011);
  ASSERT_FALSE(val010 < val001);
  ASSERT_FALSE(val001 < this->noOffset);
  ASSERT_FALSE(val110 < val101);
}

TYPED_TEST(IndexOffsetStructTests, XPlus) {
  auto offset = this->noOffset.xp();
  EXPECT_EQ(offset.dx, 1);
  for (int i = 0; i < 5; i++) {
    auto offset = this->noOffset.xp(i);
    EXPECT_EQ(offset.dx, i);
  }
}

TYPED_TEST(IndexOffsetStructTests, XMinus) {
  auto offset = this->noOffset.xm();
  EXPECT_EQ(offset.dx, -1);
  for (int i = 0; i < 5; i++) {
    auto offset = this->noOffset.xm(i);
    EXPECT_EQ(offset.dx, -i);
  }
}

TYPED_TEST(IndexOffsetStructTests, YPlus) {
  auto offset = this->noOffset.yp();
  EXPECT_EQ(offset.dy, 1);
  for (int i = 0; i < 5; i++) {
    auto offset = this->noOffset.yp(i);
    EXPECT_EQ(offset.dy, i);
  }
}

TYPED_TEST(IndexOffsetStructTests, YMinus) {
  auto offset = this->noOffset.ym();
  EXPECT_EQ(offset.dy, -1);
  for (int i = 0; i < 5; i++) {
    auto offset = this->noOffset.ym(i);
    EXPECT_EQ(offset.dy, -i);
  }
}

TYPED_TEST(IndexOffsetStructTests, ZPlus) {
  auto offset = this->noOffset.zp();
  EXPECT_EQ(offset.dz, 1);
  for (int i = 0; i < 5; i++) {
    auto offset = this->noOffset.zp(i);
    EXPECT_EQ(offset.dz, i);
  }
}

TYPED_TEST(IndexOffsetStructTests, ZMinus) {
  auto offset = this->noOffset.zm();
  EXPECT_EQ(offset.dz, -1);
  for (int i = 0; i < 5; i++) {
    auto offset = this->noOffset.zm(i);
    EXPECT_EQ(offset.dz, -i);
  }
}

TYPED_TEST(IndexOffsetStructTests, Add) {
  IndexOffset<TypeParam> offset1 = {1, 1, 1}, offset2 = {1, 2, 1}, sum = {2, 3, 2};

  EXPECT_EQ(this->noOffset + offset1, offset1);
  EXPECT_EQ(offset1 + offset2, sum);
}

TYPED_TEST(IndexOffsetStructTests, Subtract) {
  IndexOffset<TypeParam> offset1 = {1, 1, 1}, offset2 = {1, 2, 1}, difference = {0, 1, 0};
  EXPECT_EQ(offset1 - this->noOffset, offset1);
  EXPECT_EQ(offset2 - offset1, difference);
}

TYPED_TEST(IndexOffsetStructTests, AddToIndex) {
  IndexOffset<TypeParam> offset1 = {1, 0, 0}, offset2 = {0, 2, 0}, offset3 = {0, 0, 11},
                         offset4 = {2, 3, -2};
  EXPECT_EQ(this->zero + offset1, this->zero.xp());
  EXPECT_EQ(offset1 + this->zero, this->zero.xp());
  if (!std::is_same<TypeParam, IndPerp>::value) {
    EXPECT_EQ(this->zero + offset2, this->zero.yp(2));
    EXPECT_EQ(offset2 + this->zero, this->zero.yp(2));
  }
  if (!std::is_same<TypeParam, Ind2D>::value) {
    EXPECT_EQ(this->zero + offset3, this->zero.zp(11));
    EXPECT_EQ(offset3 + this->zero, this->zero.zp(11));
  }
  if (std::is_same<TypeParam, Ind3D>::value) {
    EXPECT_EQ(this->zero + offset4, this->zero.xp(2).yp(3).zm(2));
    EXPECT_EQ(offset4 + this->zero, this->zero.xp(2).yp(3).zm(2));
  }
}

TYPED_TEST(IndexOffsetStructTests, SubtractFromIndex) {
  IndexOffset<TypeParam> offset1 = {1, 0, 0}, offset2 = {0, 2, 0}, offset3 = {0, 0, 11},
                         offset4 = {2, 3, -2};
  EXPECT_EQ(this->zero - offset1, this->zero.xm());
  if (!std::is_same<TypeParam, IndPerp>::value) {
    EXPECT_EQ(this->zero - offset2, this->zero.ym(2));
  }
  if (!std::is_same<TypeParam, Ind2D>::value) {
    EXPECT_EQ(this->zero - offset3, this->zero.zm(11));
  }
  if (std::is_same<TypeParam, Ind3D>::value) {
    EXPECT_EQ(this->zero - offset4, this->zero.zp(2).xm(2).ym(3));
  }
}


template <typename T>
class StencilUnitTests :  public ::testing::Test {
public:
  WithQuietOutput all{output};
  StencilUnitTests() {
    zero = T(0, std::is_same<T, IndPerp>::value ? 1 : 5,
	     std::is_same<T, Ind2D>::value ? 1 : 7);
    for (int i = 0; i < static_cast<int>(sizes.size()); i++) {
      std::vector<IndexOffset<T>> part;
      for (int j = 0; j < sizes[i]; j++) {
	part.push_back(noOffset.xp(j));
      }
      stencil.add([i](T ind) -> bool { return ind.x() == i; }, part);
    }
  }

  virtual ~StencilUnitTests() {
  }

  IndexOffset<T> noOffset;
  T zero;
  OperatorStencil<T> stencil;
  std::vector<int> sizes = {3, 4, 5};
};

TYPED_TEST_SUITE(StencilUnitTests, IndTypes);


// Test get by integer
TYPED_TEST(StencilUnitTests, GetByInteger) {
  for (int i = 0; i < static_cast<int>(this->sizes.size()); i++) {
    const auto& part = this->stencil.getStencilPart(i);
    int j = 0;
    for (auto& off : part) {
      EXPECT_EQ(off, this->noOffset.xp(j));
      j++;
    }
    EXPECT_EQ(this->sizes[i], j);    
  }
}

// Test get by index
TYPED_TEST(StencilUnitTests, GetByindex) {
  for (int i = 0; i < static_cast<int>(this->sizes.size()); i++) {
    const auto& part = this->stencil.getStencilPart(this->zero.xp(i));
    int j = 0;
    for (auto& off : part) {
      EXPECT_EQ(off, this->noOffset.xp(j));
      j++;
    }
    EXPECT_EQ(this->sizes[i], j);
  }
}

// Test getting for an index without passing test
TYPED_TEST(StencilUnitTests, GetFailure) {
  EXPECT_THROW((this->stencil.getStencilPart(this->zero.xp(50))), BoutException);
}

// Test getting for index which can pass multiple tests
TYPED_TEST(StencilUnitTests, GetAmbiguous) {
  this->stencil.add([](TypeParam ind) -> bool { return ind.x() == 0; }, {this->noOffset});
  ASSERT_NE(this->sizes[0], 1);
  const auto& part = this->stencil.getStencilPart(this->zero);
  EXPECT_EQ(this->sizes[0], static_cast<int>(part.size()));
}

// Test get size by integer
TYPED_TEST(StencilUnitTests, GetSizeByInteger) {
  for (int i = 0; i < static_cast<int>(this->sizes.size()); i++) {
    int size = this->stencil.getStencilSize(i);
    EXPECT_EQ(size, this->sizes[i]);
  }
}

// Test get size by index
TYPED_TEST(StencilUnitTests, GetSizeByIndex) {
  for (int i = 0; i < static_cast<int>(this->sizes.size()); i++) {
    int size = this->stencil.getStencilSize(this->zero.xp(i));
    EXPECT_EQ(size, this->sizes[i]);
  }
}

// Test getting number of this->stencil-parts
TYPED_TEST(StencilUnitTests, GetNumStencils) {
  EXPECT_EQ(this->stencil.getNumParts(), static_cast<int>(this->sizes.size()));
}

// Test getting all cells whose stencils contain a given index
TYPED_TEST(StencilUnitTests, GetIndicesWithStencilIncluding) {
  for (int i = 0; i < static_cast<int>(this->sizes.size()); i++) {
    auto indicesList = this->stencil.getIndicesWithStencilIncluding(this->zero.xp(i));
    std::sort(indicesList.begin(), indicesList.end());
    EXPECT_EQ(indicesList.size(), i + 1);
    for (int j = 0; j < static_cast<int>(indicesList.size()); j++) {
      EXPECT_EQ(indicesList[j].x(), j);
      EXPECT_EQ(indicesList[j].y(), 0);
      EXPECT_EQ(indicesList[j].z(), 0);
    }
  }  
}

// Test iterator
TYPED_TEST(StencilUnitTests, Iterator) {
  int i = 0;
  for(auto& item : this->stencil) {
    EXPECT_EQ(item.part.size(), static_cast<int>(this->sizes[i]));
    i++;
  }
  EXPECT_EQ(i, static_cast<int>(this->sizes.size()));
}
