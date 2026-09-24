#include "bout/build_defines.hxx"

#if BOUT_HAS_HYPRE

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"

#include <array>
#include <cmath>
#include <set>
#include <vector>

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include "gtest/gtest.h"

#include "bout/field3d.hxx"
#include "bout/hypre_interface.hxx"

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
        indexer(std::make_shared<GlobalIndexer<T>>(bout::globals::mesh)) {}
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
  EXPECT_TRUE(IsFieldEqual(result, this->field, "RGN_NOBNDRY"));
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

  // HYPRE_IJVectorGetValues when using CUDA requires indices and values use device
  // compatible memory
#if BOUT_HAS_CUDA && defined(__CUDACC__)
  HYPRE_BigInt* um_i;
  HYPRE_Complex* um_actual;
  cudaMallocManaged(&um_i, sizeof(HYPRE_BigInt));
  cudaMallocManaged(&um_actual, sizeof(HYPRE_Complex));
  *um_i = i;
  *um_actual = actual;
  auto status = HYPRE_IJVectorGetValues(raw_vector, 1, um_i, um_actual);
  actual = *um_actual;
  cudaFree(um_i);
  cudaFree(um_actual);
#else
  auto status = HYPRE_IJVectorGetValues(raw_vector, 1, &i, &actual);
#endif

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

  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    EXPECT_EQ(vector(i), this->field[i]);
  }
}

TYPED_TEST(HypreVectorTest, SetElements) {
  HypreVector<TypeParam> vector{this->field, this->indexer};

  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    vector(i) = static_cast<BoutReal>(i.ind);
    this->field[i] = static_cast<BoutReal>(i.ind);
  }

  vector.assemble();

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

class TestParallelTransform : public ParallelTransformIdentity {
public:
  explicit TestParallelTransform(Mesh& mesh_in) : ParallelTransformIdentity(mesh_in) {}

  std::vector<PositionsAndWeights> getWeightsForYUpApproximation(int i, int j,
                                                                 int k) override {
    last_up = {i, j, k};
    up_calls += 1;
    return y_up_weights;
  }

  std::vector<PositionsAndWeights> getWeightsForYDownApproximation(int i, int j,
                                                                   int k) override {
    last_down = {i, j, k};
    down_calls += 1;
    return y_down_weights;
  }

  std::vector<PositionsAndWeights> y_up_weights;
  std::vector<PositionsAndWeights> y_down_weights;
  std::array<int, 3> last_up{-1, -1, -1};
  std::array<int, 3> last_down{-1, -1, -1};
  int up_calls{0};
  int down_calls{0};
};

template <class T>
class HypreMatrixTest : public FakeMeshFixture {
public:
  WithQuietOutput all{output};
  T field;
  IndexerPtr<T> indexer;
  TestParallelTransform* pt{nullptr};
  std::vector<ParallelTransform::PositionsAndWeights> yUpWeights, yDownWeights;
  using ind_type = typename T::ind_type;
  ind_type indexA, indexB, iWU0, iWU1, iWU2, iWD0, iWD1, iWD2;

  HypreMatrixTest()
      : FakeMeshFixture(), field(1.5, bout::globals::mesh),
        indexer(std::make_shared<GlobalIndexer<T>>(
            bout::globals::mesh, squareStencil<ind_type>(bout::globals::mesh))) {
    // indexer(std::make_shared<GlobalIndexer<T>>(bout::globals::mesh)) {

    indexA =
        ind_type((1 + field.getNy()) * field.getNz() + 1, field.getNy(), field.getNz());
    if constexpr (std::is_same_v<T, FieldPerp>) {
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

    auto transform =
        bout::utils::make_unique<TestParallelTransform>(*bout::globals::mesh);
    ParallelTransform::PositionsAndWeights wUp0 = {iWU0.x(), iWU0.y(), iWU0.z(), 0.5},
                                           wUp1 = {iWU1.x(), iWU1.y(), iWU1.z(), 1.0},
                                           wUp2 = {iWU2.x(), iWU2.y(), iWU2.z(), 0.5},
                                           wDown0 = {iWD0.x(), iWD0.y(), iWD0.z(), 0.5},
                                           wDown1 = {iWD1.x(), iWD1.y(), iWD1.z(), 1.0},
                                           wDown2 = {iWD2.x(), iWD2.y(), iWD2.z(), 0.5};
    yUpWeights = {wUp0, wUp1, wUp2};
    yDownWeights = {wDown0, wDown1, wDown2};
    transform->y_up_weights = yUpWeights;
    transform->y_down_weights = yDownWeights;
    bout::globals::mesh->getCoordinates()->setParallelTransform(std::move(transform));
    pt = dynamic_cast<TestParallelTransform*>(
        &bout::globals::mesh->getCoordinates()->getParallelTransform());
    if (pt == nullptr) {
      throw BoutException("Failed to install TestParallelTransform in HypreMatrixTest");
    }
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
  ASSERT_GE(std::pow(local_size, 2), this->field.getRegion("RGN_ALL").size());
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
  auto idx = *std::begin(region);
  BoutReal value = 23.;
  BoutReal actual = -1.;

  matrix.setVal(idx, idx, value);
  actual = matrix.getVal(idx, idx);
  EXPECT_EQ(actual, value);
}

TYPED_TEST(HypreMatrixTest, AddGetValue) {
  HypreMatrix<TypeParam> matrix(this->indexer);
  const auto& region = this->field.getRegion("RGN_NOBNDRY");
  auto idx = *std::begin(region);
  BoutReal value = 23.;
  BoutReal actual = -1.;

  // addVal should work like setVal if the element does not already exist
  matrix.addVal(idx, idx, value);
  actual = matrix.getVal(idx, idx);
  EXPECT_EQ(actual, value);
}

TYPED_TEST(HypreMatrixTest, SetAddGetValue) {
  HypreMatrix<TypeParam> matrix(this->indexer);
  const auto& region = this->field.getRegion("RGN_NOBNDRY");
  auto i = static_cast<HYPRE_BigInt>(this->indexer->getGlobal(*std::begin(region)));
  auto idx = *std::begin(region);
  BoutReal value = 23.;
  BoutReal actual = -1.;

  matrix.setVal(idx, idx, value);
  actual = matrix.getVal(idx, idx);
  EXPECT_EQ(actual, value);

  BoutReal second_value = 37.;
  matrix.addVal(idx, idx, second_value);
  actual = matrix.getVal(idx, idx);
  EXPECT_EQ(actual, value + second_value);
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
      BOUT_OMP_SAFE(critical)
      {
        HYPRE_IJMatrixGetValues(raw_matrix, 1, &ncolumns, &i_index, &j_index, &value);
      }
      if (i == j) {
        EXPECT_EQ(static_cast<BoutReal>(value),
                  static_cast<BoutReal>(this->indexer->getGlobal(i)));
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
      matrix.setVal(i, j, value);
    }
  }
  matrix.assemble();

  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    BOUT_FOR_SERIAL(j, this->field.getRegion("RGN_NOBNDRY")) {
      if (i == j) {
        EXPECT_EQ(matrix.getVal(i, j),
                  static_cast<BoutReal>(this->indexer->getGlobal(i)));
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
  HypreMatrix<TypeParam> matrix(this->indexer);

  if constexpr (std::is_same_v<TypeParam, FieldPerp>) {
    EXPECT_THROW(matrix.yup(), BoutException);
    return;
  }

  HypreMatrix<TypeParam> expected(this->indexer);
  const BoutReal value = 42.0;

  if constexpr (std::is_same_v<TypeParam, Field2D>) {
    expected(this->indexA, this->indexB) = value;
  } else {
    expected(this->indexA, this->iWU0) = this->yUpWeights[0].weight * value;
    expected(this->indexA, this->iWU1) = this->yUpWeights[1].weight * value;
    expected(this->indexA, this->iWU2) = this->yUpWeights[2].weight * value;
  }

  matrix.yup()(this->indexA, this->indexB) = value;

  if constexpr (std::is_same_v<TypeParam, Field3D>) {
    EXPECT_EQ(this->pt->up_calls, 1);
    EXPECT_EQ(this->pt->last_up[0], this->indexB.x());
    EXPECT_EQ(this->pt->last_up[1], this->indexA.y());
    EXPECT_EQ(this->pt->last_up[2], this->indexB.z());
  }

  expected.assemble();
  matrix.assemble();

  EXPECT_TRUE(IsHypreMatrixEqual(matrix, expected));
}

TYPED_TEST(HypreMatrixTest, YDown) {
  HypreMatrix<TypeParam> matrix(this->indexer);

  if constexpr (std::is_same_v<TypeParam, FieldPerp>) {
    EXPECT_THROW(matrix.yup(), BoutException);
    return;
  }

  HypreMatrix<TypeParam> expected(this->indexer);
  const BoutReal value = 42.0;

  if constexpr (std::is_same_v<TypeParam, Field2D>) {
    expected(this->indexB, this->indexA) = value;
  } else {
    expected(this->indexB, this->iWD0) = this->yDownWeights[0].weight * value;
    expected(this->indexB, this->iWD1) = this->yDownWeights[1].weight * value;
    expected(this->indexB, this->iWD2) = this->yDownWeights[2].weight * value;
  }

  matrix.ydown()(this->indexB, this->indexA) = value;

  if constexpr (std::is_same_v<TypeParam, Field3D>) {
    EXPECT_EQ(this->pt->down_calls, 1);
    EXPECT_EQ(this->pt->last_down[0], this->indexA.x());
    EXPECT_EQ(this->pt->last_down[1], this->indexB.y());
    EXPECT_EQ(this->pt->last_down[2], this->indexA.z());
  }

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

namespace {

struct RawBoundaryEliminationSystem {
  std::vector<HYPRE_Int> ncols;
  std::vector<HYPRE_BigInt> rows;
  std::vector<HYPRE_BigInt> cols;
  std::vector<HYPRE_Complex> values;
  std::vector<HYPRE_Int> boundary_rows;
  HYPRE_Int* row_indexes{nullptr};
  std::unique_ptr<bout::BoundaryElimination> elimination;

  RawBoundaryEliminationSystem(std::initializer_list<HYPRE_Int> ncols_in,
                               std::initializer_list<HYPRE_BigInt> rows_in,
                               std::initializer_list<HYPRE_BigInt> cols_in,
                               std::initializer_list<HYPRE_Complex> values_in,
                               std::initializer_list<HYPRE_Int> boundary_rows_in)
      : ncols(ncols_in), rows(rows_in), cols(cols_in), values(values_in),
        boundary_rows(boundary_rows_in) {
    elimination = std::make_unique<bout::BoundaryElimination>(
        static_cast<HYPRE_Int>(rows.size()), ncols.data(), rows.data(), &row_indexes,
        cols.data(), values.data(), static_cast<HYPRE_Int>(boundary_rows.size()),
        boundary_rows.data());
  }

  ~RawBoundaryEliminationSystem() {
    if (row_indexes != nullptr) {
      HypreFree(row_indexes);
    }
  }
};

RawBoundaryEliminationSystem makeSingleBoundarySystem() {
  return {{2, 3, 2},
          {0, 1, 2},
          {0, 1, 0, 1, 2, 1, 2},
          {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0},
          {0}};
}

RawBoundaryEliminationSystem makeMultiCoupledBoundarySystem() {
  return {{2, 3, 3, 2},
          {0, 1, 2, 3},
          {0, 1, 0, 1, 3, 0, 1, 2, 2, 3},
          {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0},
          {0}};
}

RawBoundaryEliminationSystem makeBackwardCoupledRowSystem() {
  return {{2, 3, 3, 2},
          {0, 1, 2, 3},
          {0, 2, 1, 2, 3, 0, 2, 3, 1, 3},
          {2.0, 3.0, 31.0, 37.0, 29.0, 5.0, 7.0, 11.0, 23.0, 19.0},
          {0, 3}};
}

} // namespace

TEST(BoundaryEliminationTest, TransformsSingleBoundarySystem) {
  auto system = makeSingleBoundarySystem();

  EXPECT_EQ(system.ncols[0], 1);
  EXPECT_EQ(system.values[0], -1.0);
  EXPECT_EQ(system.values[1], 0.0);
  EXPECT_EQ(system.values[2], 0.0);
  EXPECT_DOUBLE_EQ(system.values[3], -0.5);
  EXPECT_EQ(system.values[4], 11.0);
  EXPECT_EQ(system.values[5], 13.0);
  EXPECT_EQ(system.values[6], 17.0);

  ASSERT_NE(system.elimination, nullptr);
  EXPECT_EQ(system.elimination->binum_array[0], 0);
  EXPECT_EQ(system.elimination->bjnum_array[0], 1);
  EXPECT_EQ(system.elimination->bdep_array[0], -1);
  EXPECT_EQ(system.elimination->bii_array[0], 2.0);
  EXPECT_EQ(system.elimination->bij_array[0], 3.0);
  EXPECT_EQ(system.elimination->aoffset_array[0], 0);
  EXPECT_EQ(system.elimination->aoffset_array[1], 1);
  EXPECT_EQ(system.elimination->aknum_array[0], 1);
  EXPECT_EQ(system.elimination->aki_array[0], 5.0);
}

TEST(BoundaryEliminationTest, ReducesRightHandSideForSingleBoundarySystem) {
  auto system = makeSingleBoundarySystem();
  std::array<HYPRE_Complex, 3> rhs{{19.0, 23.0, 29.0}};

  auto brhs = system.elimination->reduceRightHandSideInPlace(rhs.data());

  ASSERT_NE(brhs, nullptr);
  EXPECT_EQ(brhs->data[0], 19.0);
  EXPECT_EQ(rhs[0], 19.0);
  EXPECT_DOUBLE_EQ(rhs[1], -24.5);
  EXPECT_EQ(rhs[2], 29.0);
}

TEST(BoundaryEliminationTest, ExpandsSolutionForSingleBoundarySystem) {
  auto system = makeSingleBoundarySystem();
  auto brhs = std::make_shared<bout::HypreComplexArray>(1);
  brhs->data[0] = 19.0;
  std::array<HYPRE_Complex, 3> solution{{0.0, 2.0, 3.0}};

  system.elimination->expandSolutionInPlace(brhs, solution.data());

  EXPECT_DOUBLE_EQ(solution[0], 6.5);
  EXPECT_EQ(solution[1], 2.0);
  EXPECT_EQ(solution[2], 3.0);
}

TEST(BoundaryEliminationTest, ReconstructsMatvecForSingleBoundarySystem) {
  auto system = makeSingleBoundarySystem();
  std::array<HYPRE_Complex, 3> x{{1.0, 2.0, 3.0}};
  std::array<HYPRE_Complex, 3> reduced_result{{-1.0, 32.0, 77.0}};
  auto boundary_values = system.elimination->evaluateBoundaryEquations(x.data());

  system.elimination->expandMatvecResultInPlace(boundary_values, boundary_values,
                                                reduced_result.data());

  EXPECT_EQ(reduced_result[0], 8.0);
  EXPECT_EQ(reduced_result[1], 52.0);
  EXPECT_EQ(reduced_result[2], 77.0);
}

TEST(BoundaryEliminationTest, EliminatesBoundaryCouplingsFromAllAffectedRows) {
  auto system = makeMultiCoupledBoundarySystem();

  EXPECT_EQ(system.ncols[0], 1);
  EXPECT_EQ(system.values[0], -1.0);
  EXPECT_EQ(system.values[1], 0.0);
  EXPECT_EQ(system.values[2], 0.0);
  EXPECT_DOUBLE_EQ(system.values[3], -0.5);
  EXPECT_EQ(system.values[4], 11.0);
  EXPECT_EQ(system.values[5], 0.0);
  EXPECT_DOUBLE_EQ(system.values[6], -2.5);
  EXPECT_EQ(system.values[7], 19.0);
  EXPECT_EQ(system.values[8], 23.0);
  EXPECT_EQ(system.values[9], 29.0);

  EXPECT_EQ(system.elimination->na, 2);
  EXPECT_EQ(system.elimination->aoffset_array[0], 0);
  EXPECT_EQ(system.elimination->aoffset_array[1], 2);
  EXPECT_EQ(system.elimination->aknum_array[0], 1);
  EXPECT_EQ(system.elimination->aknum_array[1], 2);
  EXPECT_EQ(system.elimination->aki_array[0], 5.0);
  EXPECT_EQ(system.elimination->aki_array[1], 13.0);
}

TEST(BoundaryEliminationTest, ReducesRightHandSideForMultipleAffectedRows) {
  auto system = makeMultiCoupledBoundarySystem();
  std::array<HYPRE_Complex, 4> rhs{{19.0, 23.0, 31.0, 37.0}};

  auto brhs = system.elimination->reduceRightHandSideInPlace(rhs.data());

  ASSERT_NE(brhs, nullptr);
  EXPECT_EQ(brhs->data[0], 19.0);
  EXPECT_EQ(rhs[0], 19.0);
  EXPECT_DOUBLE_EQ(rhs[1], -24.5);
  EXPECT_DOUBLE_EQ(rhs[2], -92.5);
  EXPECT_EQ(rhs[3], 37.0);
}

TEST(BoundaryEliminationTest, ReconstructsMatvecForMultipleAffectedRows) {
  auto system = makeMultiCoupledBoundarySystem();
  std::array<HYPRE_Complex, 4> x{{1.0, 2.0, 3.0, 4.0}};
  std::array<HYPRE_Complex, 4> reduced_result{{-1.0, 32.0, 52.0, 185.0}};
  auto boundary_values = system.elimination->evaluateBoundaryEquations(x.data());

  system.elimination->expandMatvecResultInPlace(boundary_values, boundary_values,
                                                reduced_result.data());

  EXPECT_EQ(reduced_result[0], 8.0);
  EXPECT_EQ(reduced_result[1], 52.0);
  EXPECT_EQ(reduced_result[2], 104.0);
  EXPECT_EQ(reduced_result[3], 185.0);
}

TEST(BoundaryEliminationTest, HandlesBackwardCoupledRows) {
  auto system = makeBackwardCoupledRowSystem();

  EXPECT_EQ(system.ncols[0], 1);
  EXPECT_EQ(system.ncols[3], 1);
  EXPECT_EQ(system.values[0], -1.0);
  EXPECT_EQ(system.values[1], 0.0);
  EXPECT_DOUBLE_EQ(system.values[2], 31.0 - (29.0 * 23.0 / 19.0));
  EXPECT_EQ(system.values[4], 0.0);
  EXPECT_EQ(system.values[5], 0.0);
  EXPECT_DOUBLE_EQ(system.values[6], -0.5);
  EXPECT_EQ(system.values[8], 0.0);
  EXPECT_EQ(system.values[9], -1.0);

  EXPECT_EQ(system.elimination->aoffset_array[0], 0);
  EXPECT_EQ(system.elimination->aoffset_array[1], 1);
  EXPECT_EQ(system.elimination->aoffset_array[2], 2);
  EXPECT_EQ(system.elimination->aknum_array[0], 2);
  EXPECT_EQ(system.elimination->aknum_array[1], 1);
}

TEST(BoundaryEliminationTest, DifferentEliminationObjectsProduceDifferentReducedSystems) {
  auto system_a = makeSingleBoundarySystem();
  RawBoundaryEliminationSystem system_b{{2, 3, 2},
                                        {0, 1, 2},
                                        {0, 1, 0, 1, 2, 1, 2},
                                        {4.0, -1.0, 6.0, 9.0, 10.0, 13.0, 17.0},
                                        {0}};
  std::array<HYPRE_Complex, 3> rhs_a{{19.0, 23.0, 29.0}};
  std::array<HYPRE_Complex, 3> rhs_b = rhs_a;

  auto brhs_a = system_a.elimination->reduceRightHandSideInPlace(rhs_a.data());
  auto brhs_b = system_b.elimination->reduceRightHandSideInPlace(rhs_b.data());

  EXPECT_EQ(brhs_a->data[0], 19.0);
  EXPECT_EQ(brhs_b->data[0], 19.0);
  EXPECT_DOUBLE_EQ(rhs_a[1], -24.5);
  EXPECT_DOUBLE_EQ(rhs_b[1], -5.5);
  EXPECT_NE(rhs_a[1], rhs_b[1]);
}

#endif // BOUT_HAS_HYPRE
