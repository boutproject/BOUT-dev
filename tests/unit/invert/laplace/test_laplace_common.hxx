#pragma once

#include "fake_mesh.hxx"
#include "fake_mesh_fixture.hxx"

#include "bout/bout_types.hxx"
#include "bout/build_defines.hxx"
#include "bout/derivs.hxx"
#include "bout/difops.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/invert_laplace.hxx"

#include <ostream>

#include "fmt/format.h"

#if BOUT_HAS_PETSC
#include "bout/petsc_interface.hxx"
#endif

namespace bout::testing {
// If false then use Dirichlet conditions on the corresponding boundary
struct Neumann {
  bool inner_x;
  bool outer_x;
  bool lower_y;
  bool upper_y;
};

// Teach googletest how to print a `Neumann`
inline std::ostream& operator<<(std::ostream& out, const Neumann& neumann) {
  return out << fmt::format(
             "Neumann{{inner_x: {}, outer_x: {}, lower_y: {}, upper_y: {}}}",
             neumann.inner_x, neumann.outer_x, neumann.lower_y, neumann.upper_y);
}

class ForwardOperator {
public:
  ForwardOperator(Mesh& mesh, Neumann boundaries)
      : a(0.0), c1(1.0), c2(1.0), d(1.0), ex(0.0), ez(0.0),
        coords(mesh.getCoordinates(CELL_CENTER)), neumann(boundaries) {}

  Field3D operator()(Field3D& f) {
    auto result = d * Laplace_perp(f, CELL_DEFAULT, "free", "RGN_NOY")
                  + (Grad(f) * Grad(c2) - DDY(c2) * DDY(f) / coords->g_22) / c1 + a * f
                  + ex * DDX(f) + ez * DDZ(f);
    applyBoundaries(result, f);
    return result;
  }

  Field3D a, c1, c2, d, ex, ez;

private:
  Coordinates* coords;

  Neumann neumann;

  void applyBoundaries(Field3D& newF, Field3D& f) {
    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_INNER_X")) {
      if (neumann.inner_x) {
        newF[i] = (f[i.xp()] - f[i]) / coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
        newF[i] = 0.5 * (f[i] + f[i.xp()]);
      }
    }

    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_OUTER_X")) {
      if (neumann.outer_x) {
        newF[i] = (f[i] - f[i.xm()]) / coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
        newF[i] = 0.5 * (f[i.xm()] + f[i]);
      }
    }

    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_LOWER_Y")) {
      if (neumann.lower_y) {
        newF[i] = (f[i.yp()] - f[i]) / coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
        newF[i] = 0.5 * (f[i] + f[i.yp()]);
      }
    }

    BOUT_FOR(i, f.getMesh()->getRegion3D("RGN_UPPER_Y")) {
      if (neumann.upper_y) {
        newF[i] = (f[i] - f[i.ym()]) / coords->dx[i] / sqrt(coords->g_11[i]);
      } else {
        newF[i] = 0.5 * (f[i.ym()] + f[i]);
      }
    }
  }
};

template <class LaplaceSolver>
class LaplaceTest : public FakeMeshFixture,
                    public ::testing::WithParamInterface<Neumann> {
public:
  WithQuietOutput info{output_info}, warn{output_warn}, progress{output_progress},
      all{output};
  LaplaceTest()
      : solver(&getOptions(GetParam())), forward(*bout::globals::mesh, GetParam()) {
#if BOUT_HAS_PETSC
    PetscErrorPrintf = PetscErrorPrintfNone;
#endif
    using bout::globals::mesh;
    const auto nx = static_cast<BoutReal>(mesh->GlobalNx);
    const auto ny = static_cast<BoutReal>(mesh->GlobalNy);
    const auto nz = static_cast<BoutReal>(mesh->GlobalNz);
    dynamic_cast<FakeMesh*>(bout::globals::mesh)
        ->setGridDataSource(new GridFromOptions(Options::getRoot()));
    bout::globals::mesh->getCoordinates()->geometry();
    f3.allocate();
    coef2.allocate();
    coef3.allocate();

    BOUT_FOR(i, mesh->getRegion2D("RGN_ALL")) {
      const BoutReal x = (i.x() / nx) - 0.5;
      const BoutReal y = (i.y() / ny) - 0.5;
      coef2[i] = x + y;
    }
    BOUT_FOR(i, mesh->getRegion3D("RGN_ALL")) {
      const BoutReal x = (i.x() / nx) - 0.5;
      const BoutReal y = (i.y() / ny) - 0.5;
      const BoutReal z = (i.z() / nz) - 0.5;
      f3[i] = 1e3 * exp(-0.5 * sqrt((x * x) + (y * y) + (z * z)) / sigmasq);
      coef3[i] = x + y + sin(2 * 3.14159265358979323846 * z);
    }
    // Communicate to fill in Z guards (other sides are boundaries)
    mesh->communicate(f3, coef2, coef3);
  }

  LaplaceTest(const LaplaceTest&) = delete;
  LaplaceTest(LaplaceTest&&) = delete;
  LaplaceTest& operator=(const LaplaceTest&) = delete;
  LaplaceTest& operator=(LaplaceTest&&) = delete;
  ~LaplaceTest() override {
    Options::cleanup();
#if BOUT_HAS_PETSC
    PetscErrorPrintf = PetscErrorPrintfDefault;
#endif
  }

  const BoutReal sigmasq = 0.02;
  LaplaceSolver solver;
  Field2D coef2;
  Field3D f3, coef3;
  static constexpr BoutReal tol = 1e-8;
  ForwardOperator forward;

private:
  static Options& getOptions(Neumann param) {
    Options& options = Options::root()["laplace"];
    options["type"] = "hypre3d";
    options["inner_boundary_flags"] = (param.inner_x ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    options["outer_boundary_flags"] = (param.outer_x ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    options["lower_boundary_flags"] = (param.lower_y ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    options["upper_boundary_flags"] = (param.upper_y ? INVERT_AC_GRAD : 0) + INVERT_RHS;
    options["fourth_order"] = false;
    options["atol"] = tol / 30; // Need to specify smaller than desired tolerance to
    options["rtol"] = tol / 30; // ensure it is satisfied for every element.

#if BOUT_HAS_PETSC
    auto& petsc_options = options["petsc"];
    petsc_options["mg_levels_ksp_max_it"] = 4;
    petsc_options["mg_levels_pc_type"] = "sor";
#endif

    return options;
  }
};

} // namespace bout::testing
