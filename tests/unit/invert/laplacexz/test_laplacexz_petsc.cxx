#include "bout/build_config.hxx"

#include <math.h>
#include <tuple>

#include "gtest/gtest.h"
#include "test_extras.hxx"
#include "invert_laplace.hxx"
#include "../../../../src/invert/laplacexz/impls/petsc/laplacexz-petsc.hxx"

#include "bout/mesh.hxx"
#include "bout/griddata.hxx"
#include "options.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "derivs.hxx"
#include "difops.hxx"
#include "vecops.hxx"
#include "bout/petsc_interface.hxx"

#if BOUT_HAS_PETSC

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

class ForwardOperatorXZ {
public:
  ForwardOperatorXZ() {}
  ForwardOperatorXZ(Mesh* mesh, bool xin_neumann, bool xout_neumann)
    : A(1.0, mesh), B(0.0, mesh), coords(mesh->getCoordinates(CELL_CENTER)),
        inner_x_neumann(xin_neumann), outer_x_neumann(xout_neumann) {}

  Field3D operator()(Field3D& f) const {
    const auto AJ = A * coords->J;
    const auto xx_coef = AJ * coords->g11;
    const auto zz_coef = AJ * coords->g33;
    const auto xz_coef = AJ * coords->g13;
    const auto ddx_f = DDX(f);
    const auto ddz_f = DDZ(f);

    const auto xx = (DDX(xx_coef) * ddx_f) + (xx_coef * D2DX2(f));
    const auto zz = (DDZ(zz_coef) * ddz_f) + (zz_coef * D2DZ2(f));
    const auto xz = (DDX(xz_coef) * ddz_f) + (xz_coef * D2DXDZ(f));
    const auto zx = (DDZ(xz_coef) * ddx_f) + (xz_coef * D2DXDZ(f));

    auto result = ((xx + zz + xz + zx) / coords->J) + (B * f);
    applyBoundaries(result, f);
    return result;
  }

  Field3D A, B;
  Coordinates* coords;

private:
  bool inner_x_neumann, outer_x_neumann;  // If false then use Dirichlet conditions
    // lower_y_neumann, upper_y_neumann;

  void applyBoundaries(Field3D& newF, Field3D& f) const {
    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_INNER_X")) {
      if (inner_x_neumann) {
        newF[i] = (f[i.xp()] - f[i]) / coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
        newF[i] = 0.5 * (f[i] + f[i.xp()]);
      }
    }

    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_OUTER_X")) {
      if (outer_x_neumann) {
        newF[i] = (f[i] - f[i.xm()]) / coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
        newF[i] = 0.5 * (f[i.xm()] + f[i]);
      }
    }
  }
};

class LaplaceXZPetscTest : public FakeMeshFixture,
                           public testing::WithParamInterface<std::tuple<bool, bool>> {
public:
  WithQuietOutput info{output_info}, warn{output_warn}, progress{output_progress}, all{output};
  LaplaceXZPetscTest()
      : FakeMeshFixture(), solver(bout::globals::mesh, getOptions(GetParam())),
        forward(bout::globals::mesh, std::get<0>(GetParam()), std::get<1>(GetParam())) {
    PetscErrorPrintf = PetscErrorPrintfNone;
    const BoutReal nx = mesh->GlobalNx;
    const BoutReal ny = mesh->GlobalNy;
    const BoutReal nz = mesh->GlobalNz;
    static_cast<FakeMesh*>(bout::globals::mesh)
        ->setGridDataSource(new GridFromOptions(Options::getRoot()));

    auto* coords = bout::globals::mesh->getCoordinates();

    coords->geometry();
    f3.allocate();
    A.allocate();
    B.allocate();

    BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
      const BoutReal x = i.x() / nx - 0.5;
      const BoutReal y = i.y() / ny - 0.5;
      const BoutReal z = i.z() / nz - 0.5;
      f3[i] = 1e3 * exp(-0.5 * sqrt(x * x + y * y + z * z) / sigmasq);
      A[i] = x + y + sin(2 * 3.14159265358979323846 * z);
      B[i] = 1.0;
    }
  }

  ~LaplaceXZPetscTest() {
    Options::cleanup();
    PetscErrorPrintf = PetscErrorPrintfDefault;
  }

  LaplaceXZpetsc solver;
  Field3D f3, A, B;
  static constexpr BoutReal sigmasq = 0.02;
  static constexpr BoutReal tol = 1e-8;
  ForwardOperatorXZ forward;

private:

  static Options* getOptions(std::tuple<bool, bool> param) {
    Options *options = Options::getRoot()->getSection("laplacexz");
    (*options)["type"] = "petsc";
    (*options)["inner_boundary_flags"] = (std::get<0>(param) ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    (*options)["outer_boundary_flags"] = (std::get<1>(param) ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    (*options)["fourth_order"] = false;
    (*options)["atol"] = tol/30; // Need to specify smaller than desired tolerance to
    (*options)["rtol"] = tol/30; // ensure it is satisfied for every element.
    return options;
  }
};

INSTANTIATE_TEST_SUITE_P(LaplaceXZTest, LaplaceXZPetscTest,
                         testing::Values(std::make_tuple(true, false)));

TEST_P(LaplaceXZPetscTest, TestSolve3D){
  Field3D expected = f3;
  solver.setCoefs(A, B);
  forward.A = A;
  forward.B = B;
  const Field3D actual = solver.solve(forward(f3), 0.0);
  EXPECT_TRUE(IsFieldEqual(actual, expected, "RGN_NOBNDRY", tol));
}

TEST_P(LaplaceXZPetscTest, TestSolve3DGuess){
  Field3D expected = f3, guess = f3*1.01;
  solver.setCoefs(A, B);
  forward.A = A;
  forward.B = B;
  const Field3D actual = solver.solve(forward(f3), guess);
  EXPECT_TRUE(IsFieldEqual(actual, expected, "RGN_NOBNDRY", tol));
}

#endif // BOUT_HAS_PETSC
