#include "bout/build_config.hxx"

#include <vector>

#include "test_extras.hxx"
#include "gtest/gtest.h"

#include "bout/petsc_interface.hxx"
#include "bout/region.hxx"

#if BOUT_HAS_PETSC

///////////////// Test PetscVector::Element /////////////////

class PetscVectorElementTest : public ::testing::Test {
public:
  WithQuietOutput all{output};
  Vec v;
  PetscInt n = 10;
  PetscScalar defaultVal = 1.;
  PetscLib lib;

  PetscVectorElementTest() {
    VecCreateMPI(MPI_COMM_WORLD, n, PETSC_DETERMINE, &v);
    VecSet(v, defaultVal);
    VecAssemblyBegin(v);
    VecAssemblyEnd(v);
  }

  virtual ~PetscVectorElementTest() { VecDestroy(&v); }
};

TEST_F(PetscVectorElementTest, AssignInsert) {
  PetscVector<Field3D>::Element v1(&v, 1), v2(&v, 2), v3(&v, 3), v3b(&v, 3), v9(&v, 9);
  v1 = 1.5;
  v2 = 2.5;
  v3 = 3.5;
  v3b = 4.0; // Should overwrite previous call
  v9 = 9.5;
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
  PetscScalar* vecContents;
  VecGetArray(v, &vecContents);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[0]);
  EXPECT_DOUBLE_EQ(1.5, vecContents[1]);
  EXPECT_DOUBLE_EQ(2.5, vecContents[2]);
  EXPECT_DOUBLE_EQ(4.0, vecContents[3]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[4]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[5]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[6]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[7]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[8]);
  EXPECT_DOUBLE_EQ(9.5, vecContents[9]);
}

TEST_F(PetscVectorElementTest, AssignElement) {
  PetscVector<Field3D>::Element v1(&v, 1), v2(&v, 2);
  v2 = v1 = 1.5;
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
  PetscScalar* vecContents;
  VecGetArray(v, &vecContents);
  EXPECT_DOUBLE_EQ(1.5, vecContents[1]);
  EXPECT_DOUBLE_EQ(1.5, vecContents[2]);
}

TEST_F(PetscVectorElementTest, AssignAdd) {
  PetscVector<Field3D>::Element v1(&v, 1), v2(&v, 2), v3(&v, 3), v3b(&v, 3), v9(&v, 9);
  v1 += 1.5;
  v2 += 2.5;
  v3 += 3.5;
  v3b += 4.0;
  v9 += 9.5;
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
  PetscScalar* vecContents;
  VecGetArray(v, &vecContents);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[0]);
  EXPECT_DOUBLE_EQ(2.5, vecContents[1]);
  EXPECT_DOUBLE_EQ(3.5, vecContents[2]);
  EXPECT_DOUBLE_EQ(8.5, vecContents[3]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[4]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[5]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[6]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[7]);
  EXPECT_DOUBLE_EQ(defaultVal, vecContents[8]);
  EXPECT_DOUBLE_EQ(10.5, vecContents[9]);
}

TEST_F(PetscVectorElementTest, ConvertToBoutReal) {
  PetscVector<Field3D>::Element v1(&v, 1);
  BoutReal val = v1;
  EXPECT_DOUBLE_EQ(defaultVal, val);
}

///////////////// Test PetscMatrixElement /////////////////

class PetscMatrixElementTest : public ::testing::Test {
public:
  WithQuietOutput all{output};
  Mat m;
  Vec b, x;
  PetscInt n1 = 10, n2 = 15;
  PetscScalar defaultVal = 1.;
  std::vector<PetscInt> positions;
  std::vector<BoutReal> weights;
  PetscLib lib;

  PetscMatrixElementTest() {
    PetscInt low, high;
    MatCreate(MPI_COMM_WORLD, &m);
    MatSetType(m, MATSEQDENSE);
    MatSetSizes(m, n1, n2, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetUp(m);
    MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
    VecCreateMPI(MPI_COMM_WORLD, n1, PETSC_DETERMINE, &b);
    VecCreateMPI(MPI_COMM_WORLD, n2, PETSC_DETERMINE, &x);
    VecGetOwnershipRange(x, &low, &high);
    for (int i = 0; i < n2; i++) {
      PetscScalar v = i;
      PetscInt iglobal = low + i;
      VecSetValues(x, 1, &iglobal, &v, INSERT_VALUES);
    }
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    positions.push_back(12);
    weights.push_back(-1.0);
    positions.push_back(13);
    weights.push_back(2.0);
    positions.push_back(14);
    weights.push_back(-1.0);
  }

  virtual ~PetscMatrixElementTest() {
    MatDestroy(&m);
    VecDestroy(&b);
    VecDestroy(&x);
  }
};

TEST_F(PetscMatrixElementTest, AssignInsert) {
  PetscMatrix<Field3D>::Element v1_1(&m, 1, 1), v2_3(&m, 2, 3),
      v3_13(&m, 3, 13, positions, weights), v9_6(&m, 9, 6),
      v2_13(&m, 2, 13, positions, weights), v2_14(&m, 2, 14), v4_11(&m, 4, 11),
      v4_11b(&m, 4, 11);
  v1_1 = 1.5;
  v2_3 = 1.0;
  v3_13 = -3.5;
  v9_6 = -1.0;
  v2_14 = 1.0;
  v2_13 = 0.5;
  v4_11 = 7.0;
  v4_11b = 0.5; // Should overwrite previous call
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
  MatMult(m, x, b);
  PetscScalar *bconts, *xconts;
  VecGetArray(b, &bconts);
  VecGetArray(x, &xconts);
  EXPECT_DOUBLE_EQ(0.0, bconts[0]);
  EXPECT_DOUBLE_EQ(1.5 * xconts[1], bconts[1]);
  EXPECT_DOUBLE_EQ(1.0 * xconts[3] + -1.0 * 0.5 * xconts[12] + 2.0 * 0.5 * xconts[13]
                       - 1.0 * 0.5 * xconts[14],
                   bconts[2]);
  EXPECT_DOUBLE_EQ(-1.0 * -3.5 * xconts[12] + 2.0 * -3.5 * xconts[13]
                       - 1.0 * -3.5 * xconts[14],
                   bconts[3]);
  EXPECT_DOUBLE_EQ(0.5 * xconts[11], bconts[4]);
  EXPECT_DOUBLE_EQ(0.0, bconts[5]);
  EXPECT_DOUBLE_EQ(0.0, bconts[6]);
  EXPECT_DOUBLE_EQ(0.0, bconts[7]);
  EXPECT_DOUBLE_EQ(0.0, bconts[8]);
  EXPECT_DOUBLE_EQ(-1.0 * xconts[6], bconts[9]);
}

TEST_F(PetscMatrixElementTest, AssignElement) {
  PetscMatrix<Field3D>::Element v1_1(&m, 1, 1), v1_2(&m, 1, 2);
  v1_2 = v1_1 = 1.5;
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
  MatMult(m, x, b);
  PetscScalar *bconts, *xconts;
  VecGetArray(b, &bconts);
  VecGetArray(x, &xconts);
  EXPECT_DOUBLE_EQ(1.5 * (xconts[1] + xconts[2]), bconts[1]);
}

TEST_F(PetscMatrixElementTest, AssignAdd) {
  PetscMatrix<Field3D>::Element v1_1(&m, 1, 1), v2_3(&m, 2, 3),
      v3_13(&m, 3, 13, positions, weights), v9_6(&m, 9, 6),
      v2_13(&m, 2, 13, positions, weights), v2_14(&m, 2, 14), v4_11(&m, 4, 11),
      v4_11b(&m, 4, 11);
  v1_1 += 1.5;
  v2_3 += 1.0;
  v3_13 += -3.5;
  v9_6 += -1.0;
  v2_14 += 1.0;
  v2_13 += 0.5;
  v4_11 += 7.0;
  v4_11b += 0.5; // Should overwrite previous call
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
  MatMult(m, x, b);
  PetscScalar *bconts, *xconts;
  VecGetArray(b, &bconts);
  VecGetArray(x, &xconts);
  EXPECT_DOUBLE_EQ(0.0, bconts[0]);
  EXPECT_DOUBLE_EQ(1.5 * xconts[1], bconts[1]);
  EXPECT_DOUBLE_EQ(1.0 * xconts[3] + -1.0 * 0.5 * xconts[12] + 2.0 * 0.5 * xconts[13]
                       - 1.0 * 0.5 * xconts[14] + 1.0 * xconts[14],
                   bconts[2]);
  EXPECT_DOUBLE_EQ(-1.0 * -3.5 * xconts[12] + 2.0 * -3.5 * xconts[13]
                       - 1.0 * -3.5 * xconts[14],
                   bconts[3]);
  EXPECT_DOUBLE_EQ(7.5 * xconts[11], bconts[4]);
  EXPECT_DOUBLE_EQ(0.0, bconts[5]);
  EXPECT_DOUBLE_EQ(0.0, bconts[6]);
  EXPECT_DOUBLE_EQ(0.0, bconts[7]);
  EXPECT_DOUBLE_EQ(0.0, bconts[8]);
  EXPECT_DOUBLE_EQ(-1.0 * xconts[6], bconts[9]);
}

TEST_F(PetscMatrixElementTest, ConvertToBoutReal) {
  PetscMatrix<Field3D>::Element m1_1(&m, 1, 1);
  m1_1 = 5.;
  BoutReal val = m1_1;
  EXPECT_DOUBLE_EQ(5., val);
}

#endif // BOUT_HAS_PETSC
