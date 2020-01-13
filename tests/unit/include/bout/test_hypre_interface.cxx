#include "test_extras.hxx"
#include "gtest/gtest.h"

#include "field3d.hxx"
#include "bout/hypre_interface.hxx"

#ifdef BOUT_HAS_HYPRE

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"

namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

using bout::HypreVector;

template <class T>
class HypreVectorTest : public FakeMeshFixture {
public:
  WithQuietOutput all{output};
  T field;
  HypreVectorTest() : FakeMeshFixture(), field(1.5, bout::globals::mesh) {
    GlobalIndexer::recreateGlobalInstance();
  }
  virtual ~HypreVectorTest() = default;
};

using FieldTypes = ::testing::Types<Field3D, Field2D, FieldPerp>;
TYPED_TEST_SUITE(HypreVectorTest, FieldTypes);

TYPED_TEST(HypreVectorTest, FieldConstructor) {
  BOUT_FOR(i, this->field.getRegion("RGN_ALL")) {
    this->field[i] = static_cast<BoutReal>(i.ind);
  }
  HypreVector<TypeParam> vector(this->field);
  HYPRE_BigInt jlower, jupper;
  auto hypre_vector = vector.get();
  HYPRE_IJVectorGetLocalRange(hypre_vector, &jlower, &jupper);
  const auto local_size = (jupper + 1) - jlower;
  ASSERT_EQ(local_size, this->field.getNx() * this->field.getNy() * this->field.getNz());
  const TypeParam result = vector.toField();

  EXPECT_TRUE(IsFieldEqual(this->field, result, "RGN_NOY"));
}

TYPED_TEST(HypreVectorTest, FieldAssignment) {
  HypreVector<TypeParam> vector{};

  vector = this->field;

  EXPECT_TRUE(IsFieldEqual(this->field, vector.toField(), "RGN_NOY"));
}

TYPED_TEST(HypreVectorTest, MoveConstructor) {
  HypreVector<TypeParam> vector(this->field);
  HypreVector<TypeParam> moved(std::move(vector));

  EXPECT_TRUE(IsFieldEqual(this->field, moved.toField(), "RGN_NOY"));
}

TYPED_TEST(HypreVectorTest, MoveAssignment) {
  HypreVector<TypeParam> vector{this->field};
  HypreVector<TypeParam> moved{};

  moved = std::move(vector);

  EXPECT_TRUE(IsFieldEqual(this->field, moved.toField(), "RGN_NOY"));
}

TYPED_TEST(HypreVectorTest, GetElements) {
  BOUT_FOR(i, this->field.getRegion("RGN_ALL")) {
    this->field[i] = static_cast<BoutReal>(i.ind);
  }
  HypreVector<TypeParam> vector(this->field);

  BOUT_FOR(i, this->field.getRegion("RGN_NOY")) { EXPECT_EQ(vector(i), this->field[i]); }
}

TYPED_TEST(HypreVectorTest, SetElements) {
  HypreVector<TypeParam> vector{*bout::globals::mesh};

  BOUT_FOR(i, this->field.getRegion("RGN_NOY")) {
    vector(i) = static_cast<BoutReal>(i.ind);
    // Set to identical values, but only "coincidentally"
    this->field[i] = static_cast<BoutReal>(i.ind);
  }

  EXPECT_TRUE(IsFieldEqual(this->field, vector.toField(), "RGN_NOY"));
}

#if CHECKLEVEL >= 1
TYPED_TEST(HypreVectorTest, TestGetUninitialised) {
  HypreVector<TypeParam> vector;
  typename TypeParam::ind_type index(0);
  EXPECT_THROW(vector(index), BoutException);
}

TYPED_TEST(HypreVectorTest, OutOfRange) {
  HypreVector<TypeParam> vector{this->field};
  typename TypeParam::ind_type index1(this->field.getNx() * this->field.getNy()
                                      * this->field.getNz());
  EXPECT_THROW(vector(index1), BoutException);
  typename TypeParam::ind_type index2(-1);
  EXPECT_THROW(vector(index2), BoutException);
  typename TypeParam::ind_type index3(10000000);
  EXPECT_THROW(vector(index3), BoutException);
}
#endif

TYPED_TEST(HypreVectorTest, Swap) {
  HypreVector<TypeParam> vector{this->field};
  TypeParam field2(2., bout::globals::mesh);
  HypreVector<TypeParam> vector2{field2};

  swap(vector, vector2);

  EXPECT_TRUE(IsFieldEqual(vector.toField(), field2, "RGN_NOY"));
  EXPECT_TRUE(IsFieldEqual(vector2.toField(), this->field, "RGN_NOY"));
}

#endif // BOUT_HAS_HYPRE
