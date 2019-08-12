#include <vector>

#include "gtest/gtest.h"
#include "test_extras.hxx"

#include "bout/petsc_interface.hxx"
#include "bout/region.hxx"

#ifdef BOUT_HAS_PETSC

///////////////// Test PetscVectorElement /////////////////

class PetscVectorElementTest : public ::testing::Test {
public:
  Vec v;
  PetscInt n = 10;
  PetscScalar defaultVal = 1.;
  PetscLogEvent USER_EVENT = 0;

  PetscVectorElementTest() {
    char help[] = "BOUT++: Uses finite difference methods to solve "
      "plasma fluid problems in curvilinear coordinates";
    int *pargc = nullptr;
    char ***pargv = nullptr;
    PetscInitialize(pargc,pargv,PETSC_NULL,help);
    PetscLogEventRegister("Total BOUT++",0,&USER_EVENT);
    PetscLogEventBegin(USER_EVENT,0,0,0,0);

    VecCreateMPI(MPI_COMM_WORLD, n, PETSC_DETERMINE, &v);
    VecSet(v, defaultVal);
  }

  ~PetscVectorElementTest() {
    VecDestroy(&v);
    PetscLogEventEnd(USER_EVENT,0,0,0,0);
    PetscFinalize();
  }
};

TEST_F(PetscVectorElementTest, AssignInsert) {
  PetscVectorElement& v1 = PetscVectorElement::newElement(&v, 1),
                    & v2 = PetscVectorElement::newElement(&v, 2),
                    & v3 = PetscVectorElement::newElement(&v, 3),
                    & v3b = PetscVectorElement::newElement(&v, 3),
                    & v9 = PetscVectorElement::newElement(&v, 9);
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
  free(vecContents);
}

TEST_F(PetscVectorElementTest, AssignAdd) {
  PetscVectorElement& v1 = PetscVectorElement::newElement(&v, 1),
                    & v2 = PetscVectorElement::newElement(&v, 2),
                    & v3 = PetscVectorElement::newElement(&v, 3),
                    & v3b = PetscVectorElement::newElement(&v, 3),
                    & v9 = PetscVectorElement::newElement(&v, 9);
  v1 += 1.5;
  v2 += 2.5;
  v3 += 3.5;
  v3b = 4.0;
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
  free(vecContents);
}

TEST_F(PetscVectorElementTest, NoMixedSetting) {
  
}

///////////////// Test PetscMatrixElement /////////////////

class PetscMatrixElementTest : public ::testing::Test {
public:
  Mat m;
  Vec b, x;
  PetscInt n1 = 10, n2 = 15;
  PetscScalar defaultVal = 1.;
  std::vector<PetscInt> positions;
  std::vector<BoutReal> weights;
  PetscLogEvent USER_EVENT = 0;
  
  PetscMatrixElementTest() {
    PetscInt low, high;
    char help[] = "BOUT++: Uses finite difference methods to solve "
      "plasma fluid problems in curvilinear coordinates";
    int *pargc = nullptr;
    char ***pargv = nullptr;
    PetscInitialize(pargc,pargv,PETSC_NULL,help);
    PetscLogEventRegister("Total BOUT++",0,&USER_EVENT);
    PetscLogEventBegin(USER_EVENT,0,0,0,0);

    MatCreate(MPI_COMM_WORLD, &m);
    MatSetSizes(m, n1, n2, PETSC_DETERMINE, PETSC_DETERMINE);
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

  ~PetscMatrixElementTest() {
    MatDestroy(&m);
    VecDestroy(&b);
    VecDestroy(&x);
    PetscLogEventEnd(USER_EVENT,0,0,0,0);
    PetscFinalize();
  }
};

TEST_F(PetscMatrixElementTest, AssignInsert) {
  PetscMatrixElement& v1_1 = PetscMatrixElement::newElement(&m, 1, 1),
                    & v2_3 = PetscMatrixElement::newElement(&m, 2, 3),
                    & v3_13 = PetscMatrixElement::newElement(&m, 3, 13,
							     positions,
							     weights),
                    & v9_6 = PetscMatrixElement::newElement(&m, 9, 6),
                    & v2_13 = PetscMatrixElement::newElement(&m, 2, 13,
							     positions,
							     weights),
                    & v2_14 = PetscMatrixElement::newElement(&m, 2, 14),
                    & v4_11 = PetscMatrixElement::newElement(&m, 4, 11),
                    & v4_11b = PetscMatrixElement::newElement(&m, 4, 11);
  v1_1 = 1.5;
  v2_3 = 1.0;
  v3_13 = -3.5;
  v9_6 = -1.0;
  v2_14 = 1.0;
  v2_13 = 0.5;
  v4_11 = 7.0;
  v4_11b = 0.5; // Should overwrite previous call
  MatMult(m, x, b);
  PetscScalar *bconts, *xconts;
  VecGetArray(b, &bconts);
  VecGetArray(x, &xconts);
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
  EXPECT_DOUBLE_EQ(0.0, bconts[0]);
  EXPECT_DOUBLE_EQ(1.5*xconts[1], bconts[1]);
  EXPECT_DOUBLE_EQ(1.0*xconts[3] + -1.0*0.5*xconts[12] + 2.0*0.5*xconts[13]
		   - 1.0*0.5*xconts[14], bconts[2]);
  EXPECT_DOUBLE_EQ(-1.0*-3.5*xconts[12] + 2.0*-3.5*xconts[13]
		   - 1.0*-3.5*xconts[14], bconts[3]);
  EXPECT_DOUBLE_EQ(0.5*xconts[11], bconts[4]);
  EXPECT_DOUBLE_EQ(0.0, bconts[5]);
  EXPECT_DOUBLE_EQ(0.0, bconts[6]);
  EXPECT_DOUBLE_EQ(0.0, bconts[7]);
  EXPECT_DOUBLE_EQ(0.0, bconts[8]);
  EXPECT_DOUBLE_EQ(-1.0*xconts[6], bconts[9]);
  free(bconts);
  free(xconts);
}

TEST_F(PetscMatrixElementTest, AssignAdd) {
  PetscMatrixElement& v1_1 = PetscMatrixElement::newElement(&m, 1, 1),
                    & v2_3 = PetscMatrixElement::newElement(&m, 2, 3),
                    & v3_13 = PetscMatrixElement::newElement(&m, 3, 13,
							     positions,
							     weights),
                    & v9_6 = PetscMatrixElement::newElement(&m, 9, 6),
                    & v2_13 = PetscMatrixElement::newElement(&m, 2, 13,
							     positions,
							     weights),
                    & v2_14 = PetscMatrixElement::newElement(&m, 2, 14),
                    & v4_11 = PetscMatrixElement::newElement(&m, 4, 11),
                    & v4_11b = PetscMatrixElement::newElement(&m, 4, 11);
  v1_1 += 1.5;
  v2_3 += 1.0;
  v3_13 += -3.5;
  v9_6 += -1.0;
  v2_14 += 1.0;
  v2_13 += 0.5;
  v4_11 += 7.0;
  v4_11b += 0.5; // Should overwrite previous call
  MatMult(m, x, b);
  PetscScalar *bconts, *xconts;
  VecGetArray(b, &bconts);
  VecGetArray(x, &xconts);
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
  EXPECT_DOUBLE_EQ(0.0, bconts[0]);
  EXPECT_DOUBLE_EQ(1.5*xconts[1], bconts[1]);
  EXPECT_DOUBLE_EQ(1.0*xconts[3] + -1.0*0.5*xconts[12] + 2.0*0.5*xconts[13]
		   - 1.0*0.5*xconts[14]  + 1.0*xconts[14], bconts[2]);
  EXPECT_DOUBLE_EQ(-1.0*-3.5*xconts[12] + 2.0*-3.5*xconts[13]
		   - 1.0*-3.5*xconts[14], bconts[3]);
  EXPECT_DOUBLE_EQ(7.5*xconts[11], bconts[4]);
  EXPECT_DOUBLE_EQ(0.0, bconts[5]);
  EXPECT_DOUBLE_EQ(0.0, bconts[6]);
  EXPECT_DOUBLE_EQ(0.0, bconts[7]);
  EXPECT_DOUBLE_EQ(0.0, bconts[8]);
  EXPECT_DOUBLE_EQ(-1.0*xconts[6], bconts[9]);
  free(bconts);
  free(xconts);
}

#endif // BOUT_HAS_PETSC
