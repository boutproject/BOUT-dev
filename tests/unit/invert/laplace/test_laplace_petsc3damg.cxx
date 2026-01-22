#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include <math.h>

#include "../../../../src/invert/laplace/impls/petsc3damg/petsc3damg.hxx"

#include "test_extras.hxx"
#include "test_laplace_common.hxx"

#include "bout/invert_laplace.hxx"
#include "gtest/gtest.h"

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/griddata.hxx"
#include "bout/mesh.hxx"
#include "bout/options.hxx"
#include "bout/petsc_interface.hxx"

// The unit tests use the global mesh
using namespace bout::globals;
using namespace bout::testing;

using Petsc3dAmgTest = LaplaceTest<LaplacePetsc3dAmg>;

INSTANTIATE_TEST_SUITE_P(LaplacePetsc3dAmgTest, Petsc3dAmgTest,
                         testing::Values(Neumann{false, false, false, false},
                                         Neumann{false, false, false, true},
                                         Neumann{false, false, true, false},
                                         Neumann{false, true, false, false},
                                         Neumann{true, false, false, false}));

TEST_P(Petsc3dAmgTest, TestMatrixConstruction3D) {
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefA_2D) {
  solver.setCoefA(coef2);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.a = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefA_3D) {
  solver.setCoefA(coef3);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.a = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefC_2D) {
  solver.setCoefC(coef2);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.c1 = coef2;
  forward.c2 = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefC_3D) {
  solver.setCoefC(coef3);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.c1 = coef3;
  forward.c2 = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefC1_2D) {
  solver.setCoefC1(coef2);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.c1 = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefC1_3D) {
  solver.setCoefC1(coef3);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.c1 = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefC2_2D) {
  solver.setCoefC2(coef2);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.c2 = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefC2_3D) {
  solver.setCoefC2(coef3);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.c2 = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefD_2D) {
  solver.setCoefD(coef2);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.d = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefD_3D) {
  solver.setCoefD(coef3);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.d = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefEx_2D) {
  solver.setCoefEx(coef2);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.ex = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefEx_3D) {
  solver.setCoefEx(coef3);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.ex = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefEz_2D) {
  solver.setCoefEz(coef2);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.ez = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSetCoefEz_3D) {
  solver.setCoefEz(coef3);
  PetscMatrix<Field3D>& matrix = solver.getMatrix3D();
  PetscVector<Field3D> vector(f3, solver.getIndexer());
  forward.ez = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSolve3D) {
  Field3D expected = f3;
  const Field3D actual = solver.solve(forward(f3));
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSolve3DGuess) {
  const Field3D expected = f3;
  const Field3D guess = f3 * 1.01;
  const Field3D actual = solver.solve(forward(f3), guess);
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(Petsc3dAmgTest, TestSolvePerp) {
  FieldPerp f(1.0);
  EXPECT_THROW(solver.solve(f), BoutException);
}

#endif // BOUT_HAS_PETSC
