#include <math.h>
#include <tuple>

#include "gtest/gtest.h"
#include "test_extras.hxx"
#include "invert_laplace.hxx"
#include "../../../../src/invert/laplace/impls/petsc3damg/petsc3damg.hxx"

#include "bout/mesh.hxx"
#include "bout/griddata.hxx"
#include "options.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "derivs.hxx"
#include "difops.hxx"
#include "vecops.hxx"
#include "bout/petsc_interface.hxx"


#ifdef BOUT_HAS_PETSC

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

class ForwardOperator {
public:
  ForwardOperator() {}
  ForwardOperator(bool xin_neumann, bool xout_neumann, bool ydown_neumann,
		  bool yup_neumann) : a(0.0), c1(1.0), c2(1.0), d(1.0),
		                      ex(0.0), ez(0.0) {
    coords = mesh->getCoordinates(CELL_CENTER);
    inner_x_neumann = xin_neumann;
    outer_x_neumann = xout_neumann;
    lower_y_neumann = ydown_neumann;
    upper_y_neumann = yup_neumann;
  }

  const Field3D operator()(Field3D &f) {
    auto result = d * Laplace_perp(f) + (Grad(f)*Grad(c2)-DDY(c2)*DDY(f)/coords->g_22)/c1
      + a*f + ex*DDX(f) + ez*DDZ(f);
    applyBoundaries(result, f);
    return result;
  }

  Field3D a, c1, c2, d, ex, ez;
  Coordinates* coords;

private:
  bool inner_x_neumann, outer_x_neumann,  // If false then use Dirichlet conditions
    lower_y_neumann, upper_y_neumann;

  void applyBoundaries(Field3D& newF, Field3D &f) {
    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_INNER_X")) {
      if (inner_x_neumann) {
	newF[i] = (f[i.xp()] - f[i])/ coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
	newF[i] = 0.5 * (f[i] + f[i.xp()]);
      }
    }

    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_OUTER_X")) {
      if (outer_x_neumann) {
	newF[i] = (f[i] - f[i.xm()])/ coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
	newF[i] = 0.5 * (f[i.xm()] + f[i]);
      }
    }

    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_LOWER_Y")) {
      if (lower_y_neumann) {
	newF[i] = (f[i.yp()] - f[i])/ coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
	newF[i] = 0.5 * (f[i] + f[i.yp()]);
      }
    }

    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_UPPER_Y")) {
      if (upper_y_neumann) {
	newF[i] = (f[i] - f[i.ym()])/ coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
	newF[i] = 0.5 * (f[i.ym()] + f[i]);
      }
    }
  }
};


class Petsc3dAmgTest : public FakeMeshFixture,
		       public testing::WithParamInterface<std::tuple<bool, bool,
								     bool, bool > > {
public:
  Petsc3dAmgTest() : FakeMeshFixture() {
    int nx = mesh->GlobalNx,
        ny = mesh->GlobalNy,
        nz = mesh->GlobalNz;
    BoutReal x, y, z;
    static_cast<FakeMesh*>(bout::globals::mesh)->setGridDataSource(new GridFromOptions(Options::getRoot()));
    bout::globals::mesh->getCoordinates()->geometry();
    f3.allocate();
    coef2.allocate();
    coef3.allocate();

    BOUT_FOR(i, mesh->getRegion2D("RGN_ALL")) {
      x = i.x()/(BoutReal)nx - 0.5;
      y = i.y()/(BoutReal)ny - 0.5;
      coef2[i] = x + y;
    }
    BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
      x = i.x()/(BoutReal)nx - 0.5;
      y = i.y()/(BoutReal)ny - 0.5;
      z = i.z()/(BoutReal)nz - 0.5;
      f3[i] = 1e3*exp(-0.5*sqrt(x*x + y*y + z*z)/sigmasq);
      coef3[i] = x + y + sin(2*3.14159265358979323846*z);
    }
    auto param = GetParam();
    forward = ForwardOperator(std::get<0>(param), std::get<1>(param),
			      std::get<2>(param), std::get<3>(param));
    Options *options = Options::getRoot()->getSection("laplace");
    (*options)["type"] = "petsc3damg";
    (*options)["inner_boundary_flags"] = (std::get<0>(param) ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    (*options)["outer_boundary_flags"] = (std::get<1>(param) ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    (*options)["lower_boundary_flags"] = (std::get<2>(param) ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    (*options)["upper_boundary_flags"] = (std::get<3>(param) ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    (*options)["fourth_order"] = false;
    (*options)["atol"] = tol/10; // Need to specify smaller than desired tolerance to
    (*options)["rtol"] = tol/10; // ensure it is satisfied for every element.
    //    (*options)["dtol"] = 1e5;
    //    (*options)["maxits"] = 100000;
    //    (*options)["richardson_damping_factor"] = 1.0;
    //    (*options)["chebyshev_max"] = 100.0;
    //    (*options)["chebyshev_min"] = 0.01;
    //    (*options)["gmres_max_steps"] = 30;
    solver = static_cast<LaplacePetsc3dAmg*>(Laplacian::create());
  }

  ~Petsc3dAmgTest() {
    delete solver;
    Options::cleanup();
  }

  const BoutReal sigmasq = 0.02;
  LaplacePetsc3dAmg *solver;
  Field2D coef2;
  Field3D f3, coef3;
  const BoutReal tol = 1e-8;
  ForwardOperator forward;
};

INSTANTIATE_TEST_SUITE_P(LaplacePetsc3dAmgTest, Petsc3dAmgTest,
			 testing::Values(std::make_tuple(false, false, false, false),
					 std::make_tuple(false, false, false, true),
					 std::make_tuple(false, false, true, false),
					 std::make_tuple(false, true, false, false),
					 std::make_tuple(true, false, false, false)));

TEST_P(Petsc3dAmgTest, TestMatrixConstruction3D){
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefA_2D){
  solver->setCoefA(coef2);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.a = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefA_3D){
  solver->setCoefA(coef3);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.a = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefC_2D){
  solver->setCoefC(coef2);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.c1 = coef2;
  forward.c2 = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefC_3D){
  solver->setCoefC(coef3);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.c1 = coef3;
  forward.c2 = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }

}

TEST_P(Petsc3dAmgTest, TestSetCoefC1_2D){
  solver->setCoefC1(coef2);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.c1 = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefC1_3D){
  solver->setCoefC1(coef3);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.c1 = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }

}

TEST_P(Petsc3dAmgTest, TestSetCoefC2_2D){
  solver->setCoefC2(coef2);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.c2 = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefC2_3D){
  solver->setCoefC2(coef3);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.c2 = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefD_2D){
  solver->setCoefD(coef2);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.d = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefD_3D){
  solver->setCoefD(coef3);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.d = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefEx_2D){
  solver->setCoefEx(coef2);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.ex = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefEx_3D){
  solver->setCoefEx(coef3);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.ex = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefEz_2D){
  solver->setCoefEz(coef2);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.ez = coef2;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSetCoefEz_3D){
  solver->setCoefEz(coef3);
  PetscMatrix<Field3D> &matrix = solver->getMatrix3D();
  PetscVector<Field3D> vector(f3);
  forward.ez = coef3;
  Field3D expected = forward(f3);
  Field3D actual = (matrix * vector).toField();
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSolve3D){
  Field3D expected = f3;
  const Field3D actual = solver->solve(forward(f3));
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSolve3DGuess){
  Field3D expected = f3, guess = f3;// + 0.1;
  const Field3D actual = solver->solve(forward(f3), guess);
  BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
    EXPECT_NEAR(expected[i], actual[i], tol);
  }
}

TEST_P(Petsc3dAmgTest, TestSolvePerp){
  FieldPerp f(1.0);
  EXPECT_THROW(solver->solve(f), BoutException);
}

#endif // BOUT_HAS_PETSC
