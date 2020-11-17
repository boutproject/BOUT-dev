#include <cmath>
#include "test_extras.hxx"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "field3d.hxx"
#include "bout/hypre_interface.hxx"

#if BOUT_HAS_HYPRE

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"

namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

using bout::HypreMatrix;
using bout::HypreVector;

template <class T>
class HypreVectorTest : public FakeMeshFixture {
public:
  WithQuietOutput all{output};
  T field;
  IndexerPtr<T> indexer;
  HypreVectorTest()
      : FakeMeshFixture(), field(1.5, bout::globals::mesh),
        indexer(std::make_shared<GlobalIndexer<T>>(bout::globals::mesh)) {
  }
  virtual ~HypreVectorTest() = default;
};

using FieldTypes = ::testing::Types<Field3D, Field2D, FieldPerp>;
TYPED_TEST_SUITE(HypreVectorTest, FieldTypes);

TYPED_TEST(HypreVectorTest, FieldConstructor) {
  BOUT_FOR(i, this->field.getRegion("RGN_ALL")) {
    this->field[i] = static_cast<BoutReal>(i.ind);
  }
  
  HypreVector<TypeParam> vector(this->field, this->indexer);
  HYPRE_BigInt jlower, jupper;
  auto hypre_vector = vector.get();
  HYPRE_IJVectorGetLocalRange(hypre_vector, &jlower, &jupper);
  const auto local_size = (jupper + 1) - jlower;
  ASSERT_EQ(local_size, this->indexer->size());
  const TypeParam result = vector.toField();

  // Note: Indexer doesn't have a stencil, so doesn't include boundaries
  EXPECT_TRUE(IsFieldEqual(result,this->field, "RGN_NOBNDRY"));
}

TYPED_TEST(HypreVectorTest, FieldAssignmentEmptyVector) {
  HypreVector<TypeParam> vector{};
  // vector doesn't have an index set

  EXPECT_THROW(vector = this->field, BoutException);
}

TYPED_TEST(HypreVectorTest, FieldAssignment) {
  HypreVector<TypeParam> vector{this->indexer};
  vector = this->field;
  EXPECT_TRUE(IsFieldEqual(this->field, vector.toField(), "RGN_NOBNDRY"));
}

TYPED_TEST(HypreVectorTest, MoveConstructor) {
  HypreVector<TypeParam> vector(this->field, this->indexer);
  HypreVector<TypeParam> moved(std::move(vector));

  EXPECT_TRUE(IsFieldEqual(this->field, moved.toField(), "RGN_NOBNDRY"));
}

TYPED_TEST(HypreVectorTest, MoveAssignment) {
  HypreVector<TypeParam> vector{this->field, this->indexer};
  HypreVector<TypeParam> moved{};

  moved = std::move(vector);

  EXPECT_TRUE(IsFieldEqual(this->field, moved.toField(), "RGN_NOBNDRY"));
}

TYPED_TEST(HypreVectorTest, Assemble) {
  const auto& region = this->field.getRegion("RGN_NOBNDRY");
  auto field_idx = *std::begin(region);
  this->field[field_idx] = 23.;
  auto i = static_cast<HYPRE_BigInt>(this->indexer->getGlobal(field_idx));
  HypreVector<TypeParam> vector(this->field, this->indexer);
  
  vector.assemble();

  auto raw_vector = vector.get();
  HYPRE_Complex actual{-1.};
  auto status = HYPRE_IJVectorGetValues(raw_vector, 1, &i, &actual);

  if (status != 0) {
    // Not clearing the (global) error will break future calls!
    HYPRE_ClearAllErrors();
  }
  EXPECT_EQ(status, 0);
  EXPECT_EQ(actual, 23.);
}

TYPED_TEST(HypreVectorTest, GetElements) {
  BOUT_FOR(i, this->field.getRegion("RGN_ALL")) {
    this->field[i] = static_cast<BoutReal>(i.ind);
  }
  HypreVector<TypeParam> vector(this->field, this->indexer);

  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) { EXPECT_EQ(vector(i), this->field[i]); }
}

TYPED_TEST(HypreVectorTest, SetElements) {
  HypreVector<TypeParam> vector{this->indexer};

  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    vector(i) = static_cast<BoutReal>(i.ind);
    // Set to identical values, but only "coincidentally"
    this->field[i] = static_cast<BoutReal>(i.ind);
  }

  EXPECT_TRUE(IsFieldEqual(this->field, vector.toField(), "RGN_NOBNDRY"));
}

TYPED_TEST(HypreVectorTest, Swap) {
  HypreVector<TypeParam> vector{this->field, this->indexer};
  TypeParam field2(2., bout::globals::mesh);
  HypreVector<TypeParam> vector2{field2, this->indexer};

  swap(vector, vector2);

  EXPECT_TRUE(IsFieldEqual(vector.toField(), field2, "RGN_NOBNDRY"));
  EXPECT_TRUE(IsFieldEqual(vector2.toField(), this->field, "RGN_NOBNDRY"));
}

//////////////////////////////////////////////////
// HypreMatrix tests

class MockTransform : public ParallelTransformIdentity {
public:
  MockTransform(Mesh& mesh_in) : ParallelTransformIdentity(mesh_in){};
  MOCK_METHOD(std::vector<PositionsAndWeights>, getWeightsForYUpApproximation,
              (int i, int j, int k), (override));
  MOCK_METHOD(std::vector<PositionsAndWeights>, getWeightsForYDownApproximation,
              (int i, int j, int k), (override));
};

template <class T>
class HypreMatrixTest : public FakeMeshFixture {
public:
  WithQuietOutput all{output};
  T field;
  IndexerPtr<T> indexer;
  MockTransform* pt{nullptr};
  std::vector<ParallelTransform::PositionsAndWeights> yUpWeights, yDownWeights;
  using ind_type = typename T::ind_type;
  ind_type indexA, indexB, iWU0, iWU1, iWU2, iWD0, iWD1, iWD2;

  HypreMatrixTest()
      : FakeMeshFixture(), field(1.5, bout::globals::mesh),
        indexer(std::make_shared<GlobalIndexer<T>>(bout::globals::mesh, squareStencil<ind_type>(bout::globals::mesh))) {
        //indexer(std::make_shared<GlobalIndexer<T>>(bout::globals::mesh)) {

    indexA = ind_type((1 + field.getNy()) * field.getNz() + 1, field.getNy(), field.getNz());
    if (std::is_same<T, FieldPerp>::value) {
      indexB = indexA.zp();

      iWD0 = indexB.zm();
      iWD1 = indexB;
      iWD2 = indexB.zp();
    } else {
      indexB = indexA.yp();

      iWD0 = indexB.ym();
      iWD1 = indexB;
      iWD2 = indexB.yp();
    }
    iWU0 = indexB.xm();
    iWU1 = indexB;
    iWU2 = indexB.xp();

    auto transform = bout::utils::make_unique<MockTransform>(*bout::globals::mesh);
    ParallelTransform::PositionsAndWeights wUp0 = {iWU0.x(), iWU0.y(), iWU0.z(), 0.5},
                                           wUp1 = {iWU1.x(), iWU1.y(), iWU1.z(), 1.0},
                                           wUp2 = {iWU2.x(), iWU2.y(), iWU2.z(), 0.5},
                                           wDown0 = {iWD0.x(), iWD0.y(), iWD0.z(), 0.5},
                                           wDown1 = {iWD1.x(), iWD1.y(), iWD1.z(), 1.0},
                                           wDown2 = {iWD2.x(), iWD2.y(), iWD2.z(), 0.5};
    yUpWeights = {wUp0, wUp1, wUp2};
    yDownWeights = {wDown0, wDown1, wDown2};
    pt = transform.get();
    field.getCoordinates()->setParallelTransform(std::move(transform));
  }
  virtual ~HypreMatrixTest() = default;
};

using FieldTypes = ::testing::Types<Field3D, Field2D, FieldPerp>;
TYPED_TEST_SUITE(HypreMatrixTest, FieldTypes);

TYPED_TEST(HypreMatrixTest, FieldConstructor) {
  using ind_type = typename TypeParam::ind_type;

  HypreMatrix<TypeParam> matrix(this->indexer);
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  auto hypre_matrix = matrix.get();
  HYPRE_IJMatrixGetLocalRange(hypre_matrix, &ilower, &iupper, &jlower, &jupper);
  ASSERT_EQ(ilower, jlower);
  ASSERT_EQ(iupper, jupper);
  const auto local_size = (iupper + 1) - ilower;
  ASSERT_GE(std::pow(local_size,2), this->field.getRegion("RGN_ALL").size());
}

TYPED_TEST(HypreMatrixTest, MoveConstructor) {
  HypreMatrix<TypeParam> moved(this->indexer);
  HypreMatrix<TypeParam> matrix{std::move(moved)};

  EXPECT_NE(matrix.get(), nullptr);
}

TYPED_TEST(HypreMatrixTest, MoveAssignment) {
  HypreMatrix<TypeParam> moved(this->indexer);
  HypreMatrix<TypeParam> matrix;

  matrix = std::move(moved);

  EXPECT_NE(matrix.get(), nullptr);
}

TYPED_TEST(HypreMatrixTest, SetGetValue) {
  HypreMatrix<TypeParam> matrix(this->indexer);
  const auto& region = this->field.getRegion("RGN_NOBNDRY");
  auto i = static_cast<HYPRE_BigInt>(this->indexer->getGlobal(*std::begin(region)));
  auto idx = *std::begin(region);
  BoutReal value = 23.;
  BoutReal actual = -1.;

  matrix.setVal(idx, idx, value);
  actual = matrix.getVal(idx, idx);
  EXPECT_EQ(actual, value);
}

TYPED_TEST(HypreMatrixTest, SetElements) {
  HypreMatrix<TypeParam> matrix(this->indexer);

  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    matrix.setVal(i, i, static_cast<BoutReal>(this->indexer->getGlobal(i)));
  }

  matrix.assemble();

  auto raw_matrix = matrix.get();

  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    BOUT_FOR_SERIAL(j, this->field.getRegion("RGN_NOBNDRY")) {
      auto i_index = static_cast<HYPRE_BigInt>(this->indexer->getGlobal(i));
      auto j_index = static_cast<HYPRE_BigInt>(this->indexer->getGlobal(j));
      HYPRE_Int ncolumns{1};
      HYPRE_Complex value;
      BOUT_OMP(critical) {
        HYPRE_IJMatrixGetValues(raw_matrix, 1, &ncolumns, &i_index, &j_index, &value);
      }
      if (i == j) {
        EXPECT_EQ(static_cast<BoutReal>(value), static_cast<BoutReal>(this->indexer->getGlobal(i)));
      } else {
        EXPECT_EQ(value, 0.0);
      }
    }
  }
}

TYPED_TEST(HypreMatrixTest, GetElements) {
  HypreMatrix<TypeParam> matrix(this->indexer);

  auto hypre_matrix = matrix.get();
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  HYPRE_IJMatrixGetLocalRange(hypre_matrix, &ilower, &iupper, &jlower, &jupper);

  HYPRE_Int ncolumns{1};
  for (HYPRE_Int i = ilower; i <= iupper; ++i) {
    for (HYPRE_Int j = jlower; j <= jupper; ++j) {
      HYPRE_Complex value = (i == j) ? static_cast<HYPRE_Complex>(i) : HYPRE_Complex{0.0};
      matrix.setVal(i,j,value);
    }
  }
  matrix.assemble();

  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    BOUT_FOR_SERIAL(j, this->field.getRegion("RGN_NOBNDRY")) {
      if (i == j) {
        EXPECT_EQ(matrix.getVal(i, j), static_cast<BoutReal>(this->indexer->getGlobal(i)));
      } else {
        EXPECT_EQ(matrix.getVal(i, j), 0.0);
      }
    }
  }
}

template <class T>
auto IsHypreMatrixEqual(const HypreMatrix<T>& matrix, const HypreMatrix<T>& reference)
    -> ::testing::AssertionResult {
  using namespace ::testing;

  auto hypre_matrix = matrix.get();
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  HYPRE_IJMatrixGetLocalRange(hypre_matrix, &ilower, &iupper, &jlower, &jupper);

  auto reference_matrix = reference.get();
  HYPRE_BigInt ilower_ref, iupper_ref, jlower_ref, jupper_ref;
  HYPRE_IJMatrixGetLocalRange(reference_matrix, &ilower_ref, &iupper_ref, &jlower_ref,
                              &jupper_ref);

  if (ilower != ilower_ref and iupper != iupper_ref and jlower != jlower_ref
      and jupper != jupper_ref) {
    return AssertionFailure() << "HypreMatrix is wrong size:\n  expected: " << ilower_ref
                              << ":" << iupper_ref << " x " << jlower_ref << ":"
                              << jupper_ref << "\n  got: " << ilower << ":" << iupper
                              << " x " << jlower << ":" << jupper;
  }

  for (HYPRE_BigInt i = ilower; i <= iupper; ++i) {
    for (HYPRE_BigInt j = jlower; j <= jupper; ++j) {
      if (matrix(i, j) != reference(i, j)) {
        return AssertionFailure()
               << "HypreMatrix not equal at (" << i << ", " << j
               << ")\n expected: " << reference(i, j) << "\n  got: " << matrix(i, j);
      }
    }
  }

  return AssertionSuccess();
}

TYPED_TEST(HypreMatrixTest, YUp) {
  using namespace ::testing;

  HypreMatrix<TypeParam> matrix(this->indexer);

  if (std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.yup(), BoutException);
    return;
  }

  HypreMatrix<TypeParam> expected(this->indexer);
  MockTransform* transform = this->pt;
  const BoutReal value = 42.0;     

  if (std::is_same<TypeParam, Field2D>::value) {
    expected(this->indexA, this->indexB) = value;
  } else {
    EXPECT_CALL(*transform, getWeightsForYUpApproximation(
                            this->indexB.x(), this->indexA.y(), this->indexB.z()))
               .WillOnce(Return(this->yUpWeights));
    expected(this->indexA, this->iWU0) = this->yUpWeights[0].weight * value;
    expected(this->indexA, this->iWU1) = this->yUpWeights[1].weight * value;
    expected(this->indexA, this->iWU2) = this->yUpWeights[2].weight * value;
  }

  matrix.yup()(this->indexA, this->indexB) = value;

  expected.assemble();
  matrix.assemble();

  EXPECT_TRUE(IsHypreMatrixEqual(matrix, expected));
}

TYPED_TEST(HypreMatrixTest, YDown) {
  using namespace ::testing;

  HypreMatrix<TypeParam> matrix(this->indexer);

  if (std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.yup(), BoutException);
    return;
  }

  HypreMatrix<TypeParam> expected(this->indexer);
  MockTransform* transform = this->pt;
  const BoutReal value = 42.0;

  if (std::is_same<TypeParam, Field2D>::value) {
    expected(this->indexB, this->indexA) = value;
  } else {
    EXPECT_CALL(*transform, getWeightsForYDownApproximation(
                                this->indexA.x(), this->indexB.y(), this->indexA.z()))
        .WillOnce(Return(this->yDownWeights));
    expected(this->indexB, this->iWD0) = this->yDownWeights[0].weight * value;
    expected(this->indexB, this->iWD1) = this->yDownWeights[1].weight * value;
    expected(this->indexB, this->iWD2) = this->yDownWeights[2].weight * value;
  }

  matrix.ydown()(this->indexB, this->indexA) = value;

  expected.assemble();
  matrix.assemble();

  EXPECT_TRUE(IsHypreMatrixEqual(matrix, expected));
}

TYPED_TEST(HypreMatrixTest, YNext0) {
  HypreMatrix<TypeParam> matrix(this->indexer);
  HypreMatrix<TypeParam> expected(this->indexer);

  const BoutReal value = 42.0;

  matrix.ynext(0)(this->indexA, this->indexB) = value;
  expected(this->indexA, this->indexB) = value;

  expected.assemble();
  matrix.assemble();

  EXPECT_TRUE(IsHypreMatrixEqual(matrix, expected));
}

#endif // BOUT_HAS_HYPRE
