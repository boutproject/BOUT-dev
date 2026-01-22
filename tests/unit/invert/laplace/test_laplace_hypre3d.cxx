#include "bout/build_defines.hxx"

#if BOUT_HAS_HYPRE

#include <math.h>

#include "../../../../src/invert/laplace/impls/hypre3d/hypre3d_laplace.hxx"
#include "test_extras.hxx"
#include "test_laplace_common.hxx"

#include "bout/invert_laplace.hxx"
#include "gtest/gtest.h"

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/griddata.hxx"
#include "bout/hypre_interface.hxx"
#include "bout/mesh.hxx"
#include "bout/options.hxx"

// The unit tests use the global mesh
using namespace bout::globals;
using namespace bout::testing;

using LaplaceHypre3dTest = LaplaceTest<LaplaceHypre3d>;

INSTANTIATE_TEST_SUITE_P(LaplaceHypre3d, LaplaceHypre3dTest,
                         testing::Values(Neumann{false, false, false, false},
                                         Neumann{false, false, false, true},
                                         Neumann{false, false, true, false},
                                         Neumann{false, true, false, false},
                                         Neumann{true, false, false, false}));

TEST_P(LaplaceHypre3dTest, TestMatrixConstruction3D) {
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefA_2D) {
  solver.setCoefA(coef2);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.a = coef2;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefA_3D) {
  solver.setCoefA(coef3);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.a = coef3;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefC_2D) {
  solver.setCoefC(coef2);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.c1 = coef2;
  forward.c2 = coef2;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefC_3D) {
  solver.setCoefC(coef3);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.c1 = coef3;
  forward.c2 = coef3;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefC1_2D) {
  solver.setCoefC1(coef2);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.c1 = coef2;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefC1_3D) {
  solver.setCoefC1(coef3);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.c1 = coef3;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefC2_2D) {
  solver.setCoefC2(coef2);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.c2 = coef2;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefC2_3D) {
  solver.setCoefC2(coef3);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.c2 = coef3;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefD_2D) {
  solver.setCoefD(coef2);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.d = coef2;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefD_3D) {
  solver.setCoefD(coef3);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.d = coef3;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefEx_2D) {
  solver.setCoefEx(coef2);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.ex = coef2;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefEx_3D) {
  solver.setCoefEx(coef3);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.ex = coef3;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefEz_2D) {
  solver.setCoefEz(coef2);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.ez = coef2;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSetCoefEz_3D) {
  solver.setCoefEz(coef3);
  bout::HypreMatrix<Field3D>& matrix = solver.getMatrix3D();
  bout::HypreVector<Field3D> vector(f3, solver.getIndexer());
  bout::HypreVector<Field3D> result(0.0, solver.getIndexer());
  forward.ez = coef3;
  Field3D expected = forward(f3);
  matrix.computeAx(vector, result);
  Field3D actual = result.toField();
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSolve3D) {
  Field3D expected = f3;
  const Field3D actual = solver.solve(forward(f3));
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSolve3DGuess) {
  const Field3D expected = f3;
  const Field3D guess = f3 * 1.01;
  const Field3D actual = solver.solve(forward(f3), guess);
  EXPECT_TRUE(IsFieldEqual(expected, actual, "RGN_NOZ", tol));
}

TEST_P(LaplaceHypre3dTest, TestSolvePerp) {
  FieldPerp f(1.0);
  EXPECT_THROW(solver.solve(f), BoutException);
}

#endif // BOUT_HAS_HYPRE
