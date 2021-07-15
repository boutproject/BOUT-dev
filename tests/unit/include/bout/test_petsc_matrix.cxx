#include "bout/build_config.hxx"

#include <memory>
#include <utility>

#include "test_extras.hxx"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "field2d.hxx"
#include "field3d.hxx"
#include "fieldperp.hxx"
#include "bout/operatorstencil.hxx"
#include "bout/petsc_interface.hxx"
#include "bout/region.hxx"

#if BOUT_HAS_PETSC

#include <petscconf.h>

namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;
using ::testing::Return;

class MockTransform : public ParallelTransformIdentity {
public:
  MockTransform(Mesh& mesh_in) : ParallelTransformIdentity(mesh_in){};
  MOCK_METHOD(std::vector<PositionsAndWeights>, getWeightsForYUpApproximation,
              (int i, int j, int k), (override));
  MOCK_METHOD(std::vector<PositionsAndWeights>, getWeightsForYDownApproximation,
              (int i, int j, int k), (override));
};

// Reuse the "standard" fixture for FakeMesh
template <typename F>
class PetscMatrixTest : public FakeMeshFixture {
public:
  WithQuietOutput all{output};
  F field;
  using ind_type = typename F::ind_type;
  OperatorStencil<ind_type> stencil;
  IndexerPtr<F> indexer;
  MockTransform* pt;
  std::vector<ParallelTransform::PositionsAndWeights> yUpWeights, yDownWeights;
  typename F::ind_type indexA, indexB, iWU0, iWU1, iWU2, iWD0, iWD1, iWD2;

  PetscMatrixTest()
      : FakeMeshFixture(), field(bout::globals::mesh),
        stencil(squareStencil<ind_type>(bout::globals::mesh)),
        indexer(std::make_shared<GlobalIndexer<F>>(bout::globals::mesh, stencil)) {
    indexA = ind_type(field.getNy() * field.getNz() + 1, field.getNy(), field.getNz());
    if (std::is_same<F, FieldPerp>::value) {
      indexB = indexA.zp();
    } else {
      indexB = indexA.yp();
    }
    iWU0 = indexB.xm();
    iWU1 = indexB;
    iWU2 = indexB.xp();
    if (std::is_same<F, FieldPerp>::value) {
      iWD0 = indexB.zm();
      iWD1 = indexB;
      iWD2 = indexB.zp();
    } else {
      iWD0 = indexB.ym();
      iWD1 = indexB;
      iWD2 = indexB.yp();
    }
    std::unique_ptr<MockTransform> transform =
        bout::utils::make_unique<MockTransform>(*bout::globals::mesh);
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
    PetscErrorPrintf = PetscErrorPrintfNone;
  }

  virtual ~PetscMatrixTest() {
    PetscErrorPrintf = PetscErrorPrintfDefault;
  }
};

using FieldTypes = ::testing::Types<Field3D, Field2D, FieldPerp>;
TYPED_TEST_SUITE(PetscMatrixTest, FieldTypes);

void testMatricesEqual(Mat* m1, Mat* m2) {
  PetscInt n11, n12, n21, n22;
  MatGetLocalSize(*m1, &n11, &n12);
  MatGetLocalSize(*m2, &n21, &n22);
  PetscScalar val1, val2;
  ASSERT_EQ(n11, n21);
  ASSERT_EQ(n12, n22);
  for (PetscInt i = 0; i < n11; i++) {
    for (PetscInt j = 0; j < n12; j++) {
      MatGetValues(*m1, 1, &i, 1, &j, &val1);
      MatGetValues(*m2, 1, &i, 1, &j, &val2);
      EXPECT_EQ(val1, val2);
    }
  }
}

// Test copy constructor
TYPED_TEST(PetscMatrixTest, CopyConstructor) {
  SCOPED_TRACE("CopyConstructor");
  PetscMatrix<TypeParam> matrix(this->indexer);
  Mat* rawmat = matrix.get();
  const PetscInt i = 4, j = 1;
  const PetscScalar r = 3.141592;
  MatSetValues(*rawmat, 1, &i, 1, &j, &r, INSERT_VALUES);
  MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
  PetscMatrix<TypeParam> copy(matrix);
  Mat *matrixPtr = matrix.get(), *copyPtr = copy.get();
  EXPECT_NE(*matrixPtr, *copyPtr);
  testMatricesEqual(matrixPtr, copyPtr);
}

// Test move constructor
TYPED_TEST(PetscMatrixTest, MoveConstructor) {
  PetscMatrix<TypeParam> matrix(this->indexer);
  Mat* matrixPtr = matrix.get();
  EXPECT_NE(matrixPtr, nullptr);
  PetscMatrix<TypeParam> moved(std::move(matrix));
  Mat* movedPtr = moved.get();
  EXPECT_EQ(*matrixPtr, *movedPtr);
}

// Test copy assignment
TYPED_TEST(PetscMatrixTest, CopyAssignment) {
  SCOPED_TRACE("CopyAssignment");
  PetscMatrix<TypeParam> matrix(this->indexer);
  Mat* rawmat = matrix.get();
  const PetscInt i = 4, j = 1;
  const PetscScalar r = 3.141592;
  MatSetValues(*rawmat, 1, &i, 1, &j, &r, INSERT_VALUES);
  MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
  PetscMatrix<TypeParam> copy = matrix;
  Mat* copyPtr = copy.get();
  EXPECT_NE(*rawmat, *copyPtr);
  testMatricesEqual(rawmat, copyPtr);
}

// Test move assignment
TYPED_TEST(PetscMatrixTest, MoveAssignment) {
  PetscMatrix<TypeParam> matrix(this->indexer);
  Mat* matrixPtr = matrix.get();
  EXPECT_NE(matrixPtr, nullptr);
  PetscMatrix<TypeParam> moved = std::move(matrix);
  Mat* movedPtr = moved.get();
  EXPECT_EQ(*matrixPtr, *movedPtr);
}

// Test getting elements
TYPED_TEST(PetscMatrixTest, TestGetElements) {
  PetscMatrix<TypeParam> matrix(this->indexer);
  BOUT_FOR_SERIAL(i, this->field.getRegion("RGN_NOY")) {
    matrix(i, i) = static_cast<BoutReal>(i.ind);
  }
  Mat* rawmat = matrix.get();
  MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
  BOUT_FOR(i, this->field.getRegion("RGN_NOY")) {
    BOUT_FOR_SERIAL(j, this->field.getRegion("RGN_NOY")) {
      int i_ind = this->indexer->getGlobal(i);
      int j_ind = this->indexer->getGlobal(j);
      PetscScalar matContents;
      BOUT_OMP(critical) MatGetValues(*rawmat, 1, &i_ind, 1, &j_ind, &matContents);
      if (i == j) {
        EXPECT_EQ(matContents, static_cast<BoutReal>(i.ind));
      } else {
        EXPECT_EQ(matContents, 0.0);
      }
    }
  }
}

// Test getting constant elements
TYPED_TEST(PetscMatrixTest, TestGetElementsConst) {
  PetscMatrix<TypeParam> matrix(this->indexer);
  Mat* rawmat = matrix.get();
  typename TypeParam::ind_type ind = *std::begin(this->field.getRegion("RGN_NOBNDRY"));
  const PetscInt i = this->indexer->getGlobal(ind),
                 j = this->indexer->getGlobal(ind.xp());
  const PetscScalar r = 3.141592;
  MatSetValues(*rawmat, 1, &i, 1, &j, &r, INSERT_VALUES);
  MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
  const PetscMatrix<TypeParam> constMatrix = matrix;
  const BoutReal matContents = constMatrix(ind, ind.xp());
  EXPECT_EQ(matContents, r);
}

// Test assemble
TYPED_TEST(PetscMatrixTest, TestAssemble) {
  PetscMatrix<TypeParam> matrix(this->indexer);
  Mat* rawmat = matrix.get();
  const PetscInt i = 4, j = 1;
  const PetscScalar r = 3.141592;
  MatSetValues(*rawmat, 1, &i, 1, &j, &r, INSERT_VALUES);
  matrix.assemble();
  PetscScalar matContents;
  MatGetValues(*rawmat, 1, &i, 1, &j, &matContents);
  ASSERT_EQ(matContents, r);
}

#ifdef PETSC_USE_DEBUG

#if CHECKLEVEL >= 3
// Test trying to get an element that is out of bounds
TYPED_TEST(PetscMatrixTest, TestGetOutOfBounds) {
  PetscMatrix<TypeParam> matrix(this->indexer);
  typename TypeParam::ind_type indexa(-1), indexb(1), indexc(100000);
  typename TypeParam::ind_type index1(this->field.getNx() * this->field.getNy()
                                      * this->field.getNz());
  EXPECT_THROW((matrix(index1, indexa)), BoutException);
  EXPECT_THROW((matrix(index1, indexb)), BoutException);
  EXPECT_THROW((matrix(index1, indexc)), BoutException);
  typename TypeParam::ind_type index2(-1);
  EXPECT_THROW((matrix(index2, indexa)), BoutException);
  EXPECT_THROW((matrix(index2, indexb)), BoutException);
  EXPECT_THROW((matrix(index2, indexc)), BoutException);
  typename TypeParam::ind_type index3(10000000);
  EXPECT_THROW((matrix(index3, indexa)), BoutException);
  EXPECT_THROW((matrix(index3, indexb)), BoutException);
  EXPECT_THROW((matrix(index3, indexc)), BoutException);
}
#endif // CHECKLEVEL >= 3

// Test trying to use both INSERT_VALUES and ADD_VALUES
TYPED_TEST(PetscMatrixTest, TestMixedSetting) {
  PetscMatrix<TypeParam> matrix(this->indexer);
  typename TypeParam::ind_type i = *(this->field.getRegion("RGN_NOBNDRY").begin());
  typename TypeParam::ind_type j(i.ind + 1);
  matrix(i, i) = 1.0;
  EXPECT_THROW(matrix(i, j) += 1.1, BoutException);
}

// Test destroy
TYPED_TEST(PetscMatrixTest, TestDestroy) {
  PetscMatrix<TypeParam> matrix(this->indexer);
  Mat oldMat = *matrix.get();
  Mat newMat;
  PetscErrorCode err;
  matrix.destroy();
  err = MatDuplicate(oldMat, MAT_COPY_VALUES, &newMat);
  ASSERT_NE(err, 0); // If original matrix was destroyed, should not
                     // be able to duplicate it.
}

#endif // PETSC_USE_DEBUG

// Test getting yup
TYPED_TEST(PetscMatrixTest, TestYUp) {
  PetscMatrix<TypeParam> matrix(this->indexer, false),
      expected(this->indexer, false);
  MockTransform* transform = this->pt;
  SCOPED_TRACE("YUp");
  if (std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.yup(), BoutException);
  } else {
    BoutReal val = 3.141592;
    if (std::is_same<TypeParam, Field2D>::value) {
      expected(this->indexA, this->indexB) = val;
    } else if (std::is_same<TypeParam, Field3D>::value) {
      EXPECT_CALL(*transform, getWeightsForYUpApproximation(
                                  this->indexB.x(), this->indexA.y(), this->indexB.z()))
          .WillOnce(Return(this->yUpWeights));
      expected(this->indexA, this->iWU0) = this->yUpWeights[0].weight * val;
      expected(this->indexA, this->iWU1) = this->yUpWeights[1].weight * val;
      expected(this->indexA, this->iWU2) = this->yUpWeights[2].weight * val;
    }
    matrix.yup()(this->indexA, this->indexB) = val;
    Mat *rawmat = matrix.get(), *rawexp = expected.get();
    MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(*rawexp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawexp, MAT_FINAL_ASSEMBLY);
    testMatricesEqual(rawmat, rawexp);
  }
}

// Test getting ydown
TYPED_TEST(PetscMatrixTest, TestYDown) {
  PetscMatrix<TypeParam> matrix(this->indexer, false),
      expected(this->indexer, false);
  BoutReal val = 3.141592;
  MockTransform* transform = this->pt;
  SCOPED_TRACE("YDown");
  if (std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.ydown(), BoutException);
  } else {
    if (std::is_same<TypeParam, Field2D>::value) {
      expected(this->indexB, this->indexA) = val;
    } else if (std::is_same<TypeParam, Field3D>::value) {
      EXPECT_CALL(*transform, getWeightsForYDownApproximation(
                                  this->indexA.x(), this->indexB.y(), this->indexA.z()))
          .WillOnce(Return(this->yDownWeights));
      expected(this->indexB, this->iWD0) = this->yDownWeights[0].weight * val;
      expected(this->indexB, this->iWD1) = this->yDownWeights[1].weight * val;
      expected(this->indexB, this->iWD2) = this->yDownWeights[2].weight * val;
    }
    matrix.ydown()(this->indexB, this->indexA) = val;
    Mat *rawmat = matrix.get(), *rawexp = expected.get();
    MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(*rawexp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawexp, MAT_FINAL_ASSEMBLY);
    testMatricesEqual(rawmat, rawexp);
  }
}

// Test getting ynext(0)
TYPED_TEST(PetscMatrixTest, TestYNext0) {
  PetscMatrix<TypeParam> matrix(this->indexer),
      expected(this->indexer);
  BoutReal val = 3.141592;
  SCOPED_TRACE("YNext0");
  matrix.ynext(0)(this->indexA, this->indexB) = val;
  expected(this->indexA, this->indexB) = val;
  Mat rawmat = *matrix.get(), rawexp = *expected.get();
  MatAssemblyBegin(rawmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(rawmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(rawexp, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(rawexp, MAT_FINAL_ASSEMBLY);
  testMatricesEqual(&rawmat, &rawexp);
}

// Test getting ynext(1)
TYPED_TEST(PetscMatrixTest, TestYNextPos) {
  PetscMatrix<TypeParam> matrix(this->indexer, false),
      expected(this->indexer, false);
  BoutReal val = 3.141592;
  MockTransform* transform = this->pt;
  SCOPED_TRACE("YNextPos");
  if (std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.ynext(1), BoutException);
  } else {
    if (std::is_same<TypeParam, Field3D>::value) {
      EXPECT_CALL(*transform, getWeightsForYUpApproximation(
                                  this->indexB.x(), this->indexA.y(), this->indexB.z()))
          .Times(2)
          .WillRepeatedly(Return(this->yDownWeights));
    }
    matrix.ynext(1)(this->indexA, this->indexB) = val;
    expected.yup()(this->indexA, this->indexB) = val;
    Mat *rawmat = matrix.get(), *rawexp = expected.get();
    MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(*rawexp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawexp, MAT_FINAL_ASSEMBLY);
    testMatricesEqual(rawmat, rawexp);
  }
}

// Test getting ynext(-1)
TYPED_TEST(PetscMatrixTest, TestYNextNeg) {
  PetscMatrix<TypeParam> matrix(this->indexer, false),
      expected(this->indexer, false);
  BoutReal val = 3.141592;
  MockTransform* transform = this->pt;
  SCOPED_TRACE("YNextNeg");
  if (std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.ynext(-1), BoutException);
  } else {
    if (std::is_same<TypeParam, Field3D>::value) {
      EXPECT_CALL(*transform, getWeightsForYDownApproximation(
                                  this->indexA.x(), this->indexB.y(), this->indexA.z()))
          .Times(2)
          .WillRepeatedly(Return(this->yDownWeights));
    }
    matrix.ynext(-1)(this->indexB, this->indexA) = val;
    expected.ydown()(this->indexB, this->indexA) = val;
    Mat *rawmat = matrix.get(), *rawexp = expected.get();
    MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(*rawexp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawexp, MAT_FINAL_ASSEMBLY);
    testMatricesEqual(rawmat, rawexp);
  }
}

// Test swap
TYPED_TEST(PetscMatrixTest, TestSwap) {
  PetscMatrix<TypeParam> lhs(this->indexer), rhs(this->indexer);
  Mat l0 = *lhs.get(), r0 = *rhs.get();
  EXPECT_NE(l0, nullptr);
  EXPECT_NE(r0, nullptr);
  swap(lhs, rhs);
  Mat l1 = *lhs.get(), r1 = *rhs.get();
  EXPECT_NE(l0, l1);
  EXPECT_NE(r0, r1);
  EXPECT_EQ(l0, r1);
  EXPECT_EQ(r0, l1);
}

// Test matrix/vector multiplication (Identity)
TYPED_TEST(PetscMatrixTest, TestMatrixVectorMultiplyIdentity) {
  PetscMatrix<TypeParam> matrix(this->indexer);
  this->field.allocate();
  BOUT_FOR(i, this->field.getRegion("RGN_NOY")) {
    this->field[i] = static_cast<BoutReal>(i.ind);
    matrix(i, i) = 1.0;
  }
  PetscVector<TypeParam> vector(this->field, this->indexer);
  vector.assemble();
  matrix.assemble();
  PetscVector<TypeParam> product = matrix * vector;
  TypeParam prodField = product.toField();
  BOUT_FOR(i, prodField.getRegion("RGN_NOY")) {
    EXPECT_NEAR(prodField[i], this->field[i], 1.e-10);
  }
}

// Test matrix/vector multiplication (Ones)
TYPED_TEST(PetscMatrixTest, TestMatrixVectorMultiplyOnes) {
  PetscMatrix<TypeParam> matrix(this->indexer, false);
  this->field.allocate();
  PetscVector<TypeParam> vector(this->field, this->indexer);
  BoutReal total = 0.0;
  BOUT_FOR_OMP(i, this->field.getRegion("RGN_NOY"),
     	       parallel for reduction(+:total) schedule(BOUT_OPENMP_SCHEDULE)) {
    vector(i) = static_cast<BoutReal>(i.ind);
    total += i.ind;
    BOUT_FOR_SERIAL(j, this->field.getRegion("RGN_NOY")) { matrix(i, j) = 1.0; }
  }
  vector.assemble();
  matrix.assemble();
  PetscVector<TypeParam> product = matrix * vector;
  TypeParam prodField = product.toField();
  BOUT_FOR(i, prodField.getRegion("RGN_NOY")) {
    EXPECT_NEAR(prodField[i], total, 1.e-10);
  }
}

#endif // BOUT_HAS_PETSC
