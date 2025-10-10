#include "bout/build_defines.hxx"

#if not BOUT_USE_METRIC_3D

#include <math.h>
#include <tuple>

#include "../../../../src/invert/laplace/impls/cyclic/cyclic_laplace.hxx"

#include "bout/invert_laplace.hxx"
#include "gtest/gtest.h"

#include "bout/derivs.hxx"
#include "bout/difops.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/griddata.hxx"
#include "bout/mesh.hxx"
#include "bout/options.hxx"
#include "bout/vecops.hxx"

#include "fake_mesh_fixture.hxx"

// The unit tests use the global mesh
using namespace bout::globals;

class CyclicForwardOperator {
public:
  CyclicForwardOperator(bool xin_neumann, bool xout_neumann)
      : inner_x_neumann(xin_neumann), outer_x_neumann(xout_neumann),

        a(0.0), c1(1.0), c2(1.0), d(1.0), ex(0.0), ez(0.0) {
    coords = mesh->getCoordinates(CELL_CENTER);
  }

  const Field3D operator()(Field3D& f) {
    auto result = d * Delp2(f)
                  + (coords->g11 * DDX(f) + coords->g13 * DDZ(f)) * DDX(c2) / c1 + a * f
                  + ex * DDX(f) + ez * DDZ(f);
    applyBoundaries(result, f);
    return result;
  }

private:
  CyclicForwardOperator();
  bool inner_x_neumann, outer_x_neumann; // If false then use Dirichlet conditions

  void applyBoundaries(Field3D& newF, const Field3D& f) {
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

public:
  Field2D a, c1, c2, d, ex, ez;
  Coordinates* coords;
};

class CyclicTest : public FakeMeshFixture,
                   public testing::WithParamInterface<std::tuple<bool, bool, bool>> {
public:
  WithQuietOutput info{output_info}, warn{output_warn}, progress{output_progress},
      all{output};
  CyclicTest()
      : FakeMeshFixture(), param(GetParam()), solver(getOptions(param)),
        forward(std::get<1>(param), std::get<2>(param)) {
    int nx = mesh->GlobalNx, ny = mesh->GlobalNy, nz = mesh->GlobalNz;

    static_cast<FakeMesh*>(bout::globals::mesh)
        ->setGridDataSource(new GridFromOptions(Options::getRoot()));
    bout::globals::mesh->getCoordinates()->geometry();
    f3.allocate();
    coef2.allocate();
    coef3.allocate();

    BOUT_FOR(i, mesh->getRegion2D("RGN_ALL")) {
      BoutReal x = i.x() / (BoutReal)nx - 0.5;
      BoutReal y = i.y() / (BoutReal)ny - 0.5;
      coef2[i] = x + y;
    }
    BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
      BoutReal x = i.x() / (BoutReal)nx - 0.5;
      BoutReal y = i.y() / (BoutReal)ny - 0.5;
      BoutReal z = i.z() / (BoutReal)nz - 0.5;
      f3[i] = 1e3 * exp(-0.5 * sqrt(x * x + y * y + z * z) / sigmasq);
      coef3[i] = x + y + sin(2 * 3.14159265358979323846 * z);
    }
  }

  ~CyclicTest() { Options::cleanup(); }

private:
  std::tuple<bool, bool, bool> param;

  Options* getOptions(std::tuple<bool, bool, bool> param) {
    Options::root()["mesh"]["periodicX"] = std::get<0>(param);
    Options& options = Options::root()["laplace"];
    options["type"] = "cyclic";
    options["inner_boundary_flags"] =
        (std::get<1>(param) ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    options["outer_boundary_flags"] =
        (std::get<2>(param) ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    options["fourth_order"] = false;
    options["atol"] = tol / 30; // Need to specify smaller than desired tolerance to
    options["rtol"] = tol / 30; // ensure it is satisfied for every element.
    return &options;
  }

public:
  const BoutReal sigmasq = 0.02;
  LaplaceCyclic solver;
  Field2D coef2;
  Field3D f3, coef3;
  static constexpr BoutReal tol = 1e-8;
  CyclicForwardOperator forward;
};

INSTANTIATE_TEST_SUITE_P(LaplaceCyclicTest, CyclicTest,
                         testing::Values(std::make_tuple(false, false, false),
                                         std::make_tuple(false, false, true),
                                         std::make_tuple(false, true, false),
                                         std::make_tuple(true, false, false)));

TEST_P(CyclicTest, DummyTest) {
  // No test yet
}

#endif // BOUT_USE_METRIC_3D
