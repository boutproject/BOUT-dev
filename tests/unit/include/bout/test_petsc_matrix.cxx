#include <utility>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "test_extras.hxx"

#include "bout/petsc_interface.hxx"
#include "bout/region.hxx"
#include "field3d.hxx"
#include "field2d.hxx"
#include "fieldperp.hxx"

#ifdef BOUT_HAS_PETSC

namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;
using ::testing::Return;

class MockTransform : public ParallelTransformIdentity {
public:
  MockTransform(Mesh& mesh_in) : ParallelTransformIdentity(mesh_in) {};
  MOCK_METHOD(std::vector<positionsAndWeights>, getWeightsForYUpApproximation,
	      (int i, int j, int k), (override));
  MOCK_METHOD(std::vector<positionsAndWeights>, getWeightsForYDownApproximation,
	      (int i, int j, int k), (override));
};

// Reuse the "standard" fixture for FakeMesh
template <typename F>
class PetscMatrixTest : public FakeMeshFixture {
public:
  F field;
  MockTransform* pt;
  std::vector<ParallelTransform::positionsAndWeights> yUpWeights, yDownWeights;
  typename F::ind_type indexA, indexB, iWU0, iWU1, iWU2, iWD0, iWD1, iWD2;
  PetscMatrixTest() : FakeMeshFixture(), field(bout::globals::mesh) {
    indexA = (typename F::ind_type)(field.getNy() * field.getNz() + 1, field.getNy(), field.getNz());
    indexB = indexA.yp();
    iWU0 = indexB.xm();
    iWU1 = indexB;
    iWU2 = indexB.xp();
    iWD0 = indexB.ym();
    iWD0 = indexB;
    iWD0 = indexB.yp();
    unique_ptr<MockTransform> transform = bout::utils::make_unique<MockTransform>(*bout::globals::mesh);
    ParallelTransform::positionsAndWeights wUp0 = {iWU0.x(), iWU0.y(), iWU0.z(), 0.5},
      wUp1 = {iWU1.x(), iWU1.y(), iWU1.z(), 1.0},
      wUp2 = {iWU2.x(), iWU2.y(), iWU2.z(), 0.5},
      wDown0 = {iWD0.x(), iWD0.y(), iWD0.z(), 0.5},
      wDown1 = {iWD1.x(), iWD1.y(), iWD1.z(), 1.0},
      wDown2 = {iWD2.x(), iWD2.y(), iWD2.z(), 0.5};
    yUpWeights = {wUp0, wUp1, wUp2};
    yDownWeights = {wDown0, wDown1, wDown2};
    pt = transform.get();
    field.getCoordinates()->setParallelTransform(std::move(transform));
    //    field.getCoordinates()->setParallelTransform(bout::utils::make_unique<MockTransform>(*bout::globals::mesh));
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
      ASSERT_EQ(val1, val2);
    }
  }
}

// Test copy constructor
TYPED_TEST(PetscMatrixTest, CopyConstructor) {
  SCOPED_TRACE("CopyConstructor");
  PetscMatrix<TypeParam> matrix(this->field);
  PetscMatrix<TypeParam> copy(matrix);
  Mat *matrixPtr = matrix.getMatrixPointer(),
    *copyPtr = copy.getMatrixPointer();
  EXPECT_NE(matrixPtr, copyPtr);
  testMatricesEqual(matrixPtr, copyPtr);
}

// Test move constructor
TYPED_TEST(PetscMatrixTest, MoveConstructor) {
  PetscMatrix<TypeParam> matrix(this->field);
  Mat* matrixPtr = matrix.getMatrixPointer();
  EXPECT_NE(matrixPtr, nullptr);
  PetscMatrix<TypeParam> moved(std::move(matrix));
  Mat *movedPtr = moved.getMatrixPointer();
  EXPECT_EQ(matrixPtr, movedPtr);
  EXPECT_EQ(matrix.getMatrixPointer(), nullptr);
}

// Test copy assignment
TYPED_TEST(PetscMatrixTest, CopyAssignment) {
  SCOPED_TRACE("CopyAssignment");
  PetscMatrix<TypeParam> matrix(this->field);
  PetscMatrix<TypeParam> copy = matrix;
  Mat *matrixPtr = matrix.getMatrixPointer(),
    *copyPtr = copy.getMatrixPointer();
  EXPECT_NE(matrixPtr, copyPtr);
  testMatricesEqual(matrixPtr, copyPtr);
}

// Test move assignment
TYPED_TEST(PetscMatrixTest, MoveAssignment) {
  PetscMatrix<TypeParam> matrix(this->field);
  Mat* matrixPtr = matrix.getMatrixPointer();
  EXPECT_NE(matrixPtr, nullptr);
  PetscMatrix<TypeParam> moved = std::move(matrix);
  Mat *movedPtr = moved.getMatrixPointer();
  EXPECT_EQ(matrixPtr, movedPtr);
  EXPECT_EQ(matrix.getMatrixPointer(), nullptr);
}

// Test getting elements
TYPED_TEST(PetscMatrixTest, TestGetElements) {
  PetscMatrix<TypeParam> matrix(this->field);
  for (auto i: this->field) {
    matrix(i, i) = 2.5*i.ind*i.ind;
  }
  Mat *rawmat = matrix.getMatrixPointer();
  PetscInt m, n;
  PetscScalar matContents;
  MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
  MatGetLocalSize(*rawmat, &m, &n);
  for (PetscInt i = 0; i < n; i++) {
    for (PetscInt j = 0; j < m; j++) {
      MatGetValues(*rawmat, 1, &i, 1, &j, &matContents);
      if (i==j) {
	ASSERT_EQ(matContents, 2.5*i*i);
      } else {
	ASSERT_EQ(matContents, 0.0);
      }
    }
  }
}

// Test trying to get an element that is out of bounds
TYPED_TEST(PetscMatrixTest, TestGetOutOfBounds) {
  PetscMatrix<TypeParam> matrix(this->field);
  typename TypeParam::ind_type indexa(-1), indexb(1), indexc(100000);
  typename TypeParam::ind_type index1(this->field.getNx() * this->field.getNy() * this->field.getNy());
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

// Test assemble
TYPED_TEST(PetscMatrixTest, TestAssemble) {
  PetscMatrix<TypeParam> matrix(this->field);
  Mat *rawmat = matrix.getMatrixPointer();
  const PetscInt i = 4, j = 1;
  const PetscScalar r = 3.141592;
  MatSetValues(*rawmat, 1, &i, 1, &j, &r, INSERT_VALUES); 
  matrix.assemble();
  PetscScalar matContents;
  MatGetValues(*rawmat, 1, &i, 1, &j, &matContents);
  ASSERT_EQ(matContents, r);
}

// Test trying to use both INSERT_VALUES and ADD_VALUES
TYPED_TEST(PetscMatrixTest, TestMixedSetting) {
  PetscMatrix<TypeParam> matrix(this->field);
  Mat *rawmat = matrix.getMatrixPointer();
  const PetscInt i = 4, j1 = 1, j2 = 2;
  const PetscScalar r = 3.141592;
  MatSetValues(*rawmat, 1, &i, 1, &j1, &r, INSERT_VALUES); 
  MatSetValues(*rawmat, 1, &i, 1, &j2, &r, ADD_VALUES); 
  EXPECT_THROW(matrix.assemble(), BoutException);
}

// Test destroy
TYPED_TEST(PetscMatrixTest, TestDestroy) {
  PetscMatrix<TypeParam> matrix(this->field);
  Mat newMat;
  PetscErrorCode err;
  matrix.destroy();
  err = MatDuplicate(*matrix.getMatrixPointer(), MAT_COPY_VALUES, &newMat);
  ASSERT_NE(err, 0); // If original matrix was destroyed, should not
		     // be able to duplicate it.
}

// Test getting yup
TYPED_TEST(PetscMatrixTest, TestYUp) {
  PetscMatrix<TypeParam> matrix(this->field), expected(this->field);
  BoutReal val = 3.141592;
  MockTransform* transform = this->pt;
  SCOPED_TRACE("YUp");
  if constexpr(std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.yup(), BoutException);
  } else {
    EXPECT_CALL(*transform, getWeightsForYDownApproximation(this->indexA.x(),
                this->indexA.y(), this->indexA.z()))
      .WillOnce(Return(this->yDownWeights));

    matrix.yup()(this->indexA, this->indexB) = val;
    if constexpr(std::is_same<TypeParam, Field2D>::value) {
      expected(this->indexA, this->indexB) = val;
    } else if constexpr(std::is_same<TypeParam, Field2D>::value) {
      expected(this->indexA, this->iWU0) = this->yUpWeights[0].weight * val;
      expected(this->indexA, this->iWU1) = this->yUpWeights[1].weight * val;
      expected(this->indexA, this->iWU2) = this->yUpWeights[2].weight * val;
    }
    Mat *rawmat = matrix.getMatrixPointer(),
      *rawexp = expected.getMatrixPointer();
    MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(*rawexp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawexp, MAT_FINAL_ASSEMBLY);
    testMatricesEqual(rawmat, rawexp);
  }
}

// Test getting ydown
TYPED_TEST(PetscMatrixTest, TestYDown) {
  PetscMatrix<TypeParam> matrix(this->field), expected(this->field);
  BoutReal val = 3.141592;
  MockTransform* transform = this->pt;
  SCOPED_TRACE("YDown");
  if constexpr(std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.ydown(), BoutException);
  } else {
    EXPECT_CALL(*transform, getWeightsForYDownApproximation(this->indexA.x(),
                this->indexA.y(), this->indexA.z()))
      .WillOnce(Return(this->yDownWeights));

    matrix.ydown()(this->indexA, this->indexB) = val;
    if constexpr(std::is_same<TypeParam, Field2D>::value) {
      expected(this->indexA, this->indexB) = val;
    } else if constexpr(std::is_same<TypeParam, Field2D>::value) {
      expected(this->indexA, this->iWD0) = this->yDownWeights[0].weight * val;
      expected(this->indexA, this->iWD1) = this->yDownWeights[1].weight * val;
      expected(this->indexA, this->iWD2) = this->yDownWeights[2].weight * val;
    }
    Mat *rawmat = matrix.getMatrixPointer(),
      *rawexp = expected.getMatrixPointer();
    MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(*rawexp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawexp, MAT_FINAL_ASSEMBLY);
    testMatricesEqual(rawmat, rawexp);
  }
}

// Test getting ynext(0)
TYPED_TEST(PetscMatrixTest, TestYNext0) {
  PetscMatrix<TypeParam> matrix(this->field), expected(this->field);
  BoutReal val = 3.141592;
  MockTransform* transform = this->pt;
  SCOPED_TRACE("YNext0");
  EXPECT_CALL(*transform, getWeightsForYDownApproximation(this->indexA.x(),
              this->indexA.y(), this->indexA.z())).Times(0);
  EXPECT_CALL(*transform, getWeightsForYUpApproximation(this->indexA.x(),
              this->indexA.y(), this->indexA.z())).Times(0);
  matrix.ynext(0)(this->indexA, this->indexB) = val;
  expected(this->indexA, this->indexB) = val;
  Mat *rawmat = matrix.getMatrixPointer(),
      *rawexp = expected.getMatrixPointer();
  MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(*rawexp, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*rawexp, MAT_FINAL_ASSEMBLY);
  testMatricesEqual(rawmat, rawexp);
}

// Test getting ynext(1)
TYPED_TEST(PetscMatrixTest, TestYNextPos) {
  PetscMatrix<TypeParam> matrix(this->field), expected(this->field);
  BoutReal val = 3.141592;
  MockTransform* transform = this->pt;
  SCOPED_TRACE("YNextPos");
  if constexpr(std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.ynext(1), BoutException);
  } else {
    EXPECT_CALL(*transform, getWeightsForYDownApproximation(this->indexA.x(),
                this->indexA.y(), this->indexA.z()))
      .WillOnce(Return(this->yDownWeights));

    matrix.ynext(1)(this->indexA, this->indexB) = val;
    expected.yup()(this->indexA, this->indexB) = val;
    Mat *rawmat = matrix.getMatrixPointer(),
      *rawexp = expected.getMatrixPointer();
    MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(*rawexp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawexp, MAT_FINAL_ASSEMBLY);
    testMatricesEqual(rawmat, rawexp);
  }
}

// Test getting ynext(-1)
TYPED_TEST(PetscMatrixTest, TestYNextNeg) {
  PetscMatrix<TypeParam> matrix(this->field), expected(this->field);
  BoutReal val = 3.141592;
  MockTransform* transform = this->pt;
  SCOPED_TRACE("YNextNeg");
  if constexpr(std::is_same<TypeParam, FieldPerp>::value) {
    EXPECT_THROW(matrix.ynext(-1), BoutException);
  } else {
    EXPECT_CALL(*transform, getWeightsForYDownApproximation(this->indexA.x(),
                this->indexA.y(), this->indexA.z()))
      .WillOnce(Return(this->yDownWeights));
    matrix.ynext(-1)(this->indexA, this->indexB) = val;
    expected.ydown()(this->indexA, this->indexB) = val;
    Mat *rawmat = matrix.getMatrixPointer(),
      *rawexp = expected.getMatrixPointer();
    MatAssemblyBegin(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(*rawexp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*rawexp, MAT_FINAL_ASSEMBLY);
    testMatricesEqual(rawmat, rawexp);
  }
}

// Test swap
TYPED_TEST(PetscMatrixTest, TestSwap) {
  PetscMatrix<TypeParam> lhs(this->field), rhs(this->field);
  Mat *l0 = lhs.getMatrixPointer(), *r0 = rhs.getMatrixPointer();
  EXPECT_NE(l0, nullptr);
  EXPECT_NE(r0, nullptr);
  swap(lhs, rhs);
  Mat *l1 = lhs.getMatrixPointer(), *r1 = rhs.getMatrixPointer();
  EXPECT_NE(l0, l1);
  EXPECT_NE(r0, r1);
  EXPECT_EQ(l0, r1);
  EXPECT_EQ(r0, l1);
}


#endif // BOUT_HAS_PETSC
