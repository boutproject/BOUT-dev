/**************************************************************************
 * Testing Perpendicular Laplacian inversion using PETSc solvers
 *
 **************************************************************************
 * Copyright 2013 J.T. Omotani
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/constants.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/invert_laplace.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"
#include "bout/traits.hxx"

#include "fmt/core.h"

#include <cmath>
#include <string_view>

BoutReal max_error_at_ystart(const Field3D& error);
void apply_flat_boundary(Field3D& bcoef);

template <class T, class U>
void check_laplace(int test_num, std::string_view test_name, Laplacian& invert,
                   int inner_flags, int outer_flags, const T& acoef, const T& ccoef,
                   const T& dcoef, const U& bcoef, const Field3D& field, int ystart,
                   Options& dump) {
  static_assert(bout::utils::is_Field_v<T>, "check_laplace requires Field2D or Field3D");
  static_assert(bout::utils::is_Field_v<U>, "check_laplace requires Field2D or Field3D");

  invert.setInnerBoundaryFlags(inner_flags);
  invert.setOuterBoundaryFlags(outer_flags);
  invert.setCoefA(acoef);
  invert.setCoefC(ccoef);
  invert.setCoefD(dcoef);

  checkData(bcoef);

  Field3D sol;
  Field3D error;
  Field3D abs_error;
  BoutReal max_error = -1;

  try {
    sol = invert.solve(sliceXZ(bcoef, ystart));
    error = (field - sol) / field;
    abs_error = field - sol;
    max_error = max_error_at_ystart(abs(abs_error));
  } catch (BoutException& err) {
    output.write("BoutException occured in invert->solve(b1): {}\n", err.what());
  }

  output.write("\nTest {}: {}\n", test_num, test_name);
  output.write("Magnitude of maximum absolute error is {}\n", max_error);

  dump[fmt::format("a{}", test_num)] = acoef;
  dump[fmt::format("b{}", test_num)] = bcoef;
  dump[fmt::format("c{}", test_num)] = ccoef;
  dump[fmt::format("d{}", test_num)] = dcoef;
  dump[fmt::format("f{}", test_num)] = field;
  dump[fmt::format("sol{}", test_num)] = sol;
  dump[fmt::format("error{}", test_num)] = error;
  dump[fmt::format("absolute_error{}", test_num)] = abs_error;
  dump[fmt::format("max_error{}", test_num)] = max_error;
}

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);
  {
    Options* options = Options::getRoot()->getSection("petsc2nd");
    auto invert = Laplacian::create(options);
    options = Options::getRoot()->getSection("petsc4th");
    auto invert_4th = Laplacian::create(options);

    Options dump;

    // Solving equations of the form d*Delp2(f) + 1/c*Grad_perp(c).Grad_perp(f) + a*f = b for various f, a, c, d
    Field3D f1;
    Field3D a1;
    Field3D c1;
    Field3D d1;
    BoutReal p;
    BoutReal q; //Use to set parameters in constructing trial functions

    using bout::globals::mesh;

    // Only Neumann x-boundary conditions are implemented so far, so test functions should be Neumann in x and periodic in z.
    // Use Field3D's, but solver only works on FieldPerp slices, so only use 1 y-point
    const BoutReal nx = mesh->GlobalNx - 2 * mesh->xstart - 1;
    const BoutReal nz = mesh->GlobalNz;

    /////////////////////////////////////////////////////
    // Test 1: Gaussian x-profiles, 2nd order Krylov
    p = 0.39503274;
    q = 0.20974396;
    f1.allocate();
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          const BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
          const BoutReal z = BoutReal(jz) / nz;
          //make the gradients zero at both x-boundaries
          f1(jx, jy, jz) = 0. + exp(-(100. * pow(x - p, 2) + 1. - cos(2. * PI * (z - q))))
                           - 50.
                                 * (2. * p * exp(-100. * pow(-p, 2)) * x
                                    + (-p * exp(-100. * pow(-p, 2))
                                       - (1 - p) * exp(-100. * pow(1 - p, 2)))
                                          * pow(x, 2))
                                 * exp(-(1. - cos(2. * PI * (z - q))));
          ASSERT0(finite(f1(jx, jy, jz)));
        }
      }
    }
    if (mesh->firstX()) {
      for (int jx = mesh->xstart - 1; jx >= 0; jx--) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            const BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            const BoutReal z = BoutReal(jz) / nz;
            //make the gradients zero at both x-boundaries
            f1(jx, jy, jz) = 0.
                             + exp(-(60. * pow(x - p, 2) + 1. - cos(2. * PI * (z - q))))
                             - 50.
                                   * (2. * p * exp(-60. * pow(-p, 2)) * x
                                      + (-p * exp(-60. * pow(-p, 2))
                                         - (1 - p) * exp(-60. * pow(1 - p, 2)))
                                            * pow(x, 2))
                                   * exp(-(1. - cos(2. * PI * (z - q))));
            ASSERT0(finite(f1(jx, jy, jz)));
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int jx = mesh->xend + 1; jx < mesh->LocalNx; jx++) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            const BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            const BoutReal z = BoutReal(jz) / nz;
            //make the gradients zero at both x-boundaries
            f1(jx, jy, jz) = 0.
                             + exp(-(60. * pow(x - p, 2) + 1. - cos(2. * PI * (z - q))))
                             - 50.
                                   * (2. * p * exp(-60. * pow(-p, 2)) * x
                                      + (-p * exp(-60. * pow(-p, 2))
                                         - (1 - p) * exp(-60. * pow(1 - p, 2)))
                                            * pow(x, 2))
                                   * exp(-(1. - cos(2. * PI * (z - q))));
            ASSERT0(finite(f1(jx, jy, jz)));
          }
        }
      }
    }

    f1.applyBoundary("neumann");

    p = 0.512547;
    q = 0.30908712;
    d1.allocate();
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
          BoutReal z = BoutReal(jz) / nz;
          d1(jx, jy, jz) =
              1. + 0.2 * exp(-50. * pow(x - p, 2) / 4.) * sin(2. * PI * (z - q) * 3.);
        }
      }
    }
    if (mesh->firstX()) {
      for (int jx = mesh->xstart - 1; jx >= 0; jx--) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            d1(jx, jy, jz) =
                1. + 0.2 * exp(-50. * pow(x - p, 2) / 4.) * sin(2. * PI * (z - q) * 3.);
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int jx = mesh->xend + 1; jx < mesh->LocalNx; jx++) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            d1(jx, jy, jz) =
                1. + 0.2 * exp(-50. * pow(x - p, 2) / 4.) * sin(2. * PI * (z - q) * 3.);
          }
        }
      }
    }

    p = 0.18439023;
    q = 0.401089473;
    c1.allocate();
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
          BoutReal z = BoutReal(jz) / nz;
          c1(jx, jy, jz) =
              1. + 0.15 * exp(-50. * pow(x - p, 2) * 2.) * sin(2. * PI * (z - q) * 2.);
        }
      }
    }
    if (mesh->firstX()) {
      for (int jx = mesh->xstart - 1; jx >= 0; jx--) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            c1(jx, jy, jz) =
                1. + 0.15 * exp(-50. * pow(x - p, 2) * 2.) * sin(2. * PI * (z - q) * 2.);
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int jx = mesh->xend + 1; jx < mesh->LocalNx; jx++) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            c1(jx, jy, jz) =
                1. + 0.15 * exp(-50. * pow(x - p, 2) * 2.) * sin(2. * PI * (z - q) * 2.);
          }
        }
      }
    }

    p = 0.612547;
    q = 0.30908712;
    a1.allocate();
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
          BoutReal z = BoutReal(jz) / nz;
          a1(jx, jy, jz) =
              -1. + 0.1 * exp(-50. * pow(x - p, 2) * 2.5) * sin(2. * PI * (z - q) * 7.);
        }
      }
    }
    if (mesh->firstX()) {
      for (int jx = mesh->xstart - 1; jx >= 0; jx--) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            a1(jx, jy, jz) =
                -1. + 0.1 * exp(-50. * pow(x - p, 2) * 2.5) * sin(2. * PI * (z - q) * 7.);
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int jx = mesh->xend + 1; jx < mesh->LocalNx; jx++) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            a1(jx, jy, jz) =
                -1. + 0.1 * exp(-50. * pow(x - p, 2) * 2.5) * sin(2. * PI * (z - q) * 7.);
          }
        }
      }
    }

    checkData(f1);
    checkData(a1);
    checkData(c1);
    checkData(d1);

    mesh->communicate(f1, a1, c1, d1);

    Field3D b1 = d1 * Delp2(f1) + Grad_perp(c1) * Grad_perp(f1) / c1 + a1 * f1;
    apply_flat_boundary(b1);

    int test_num = 0;
    check_laplace(++test_num, "PETSc 2nd order", *invert, INVERT_AC_GRAD, INVERT_AC_GRAD, a1, c1,
                  d1, b1, f1, mesh->ystart, dump);

    /////////////////////////////////////////////////
    // Test 2: Gaussian x-profiles, 4th order Krylov

    check_laplace(++test_num, "PETSc 4th order", *invert_4th, INVERT_AC_GRAD, INVERT_AC_GRAD, a1,
                  c1, d1, b1, f1, mesh->ystart, dump);

    ////////////////////////////////////////////////////////////////////////////////////////
    // Test 3+4: Gaussian x-profiles, z-independent coefficients and compare with SPT method

    const Field2D a3 = DC(a1);
    const Field2D c3 = DC(c1);
    const Field2D d3 = DC(d1);
    Field3D b3 = d3 * Delp2(f1) + Grad_perp(c3) * Grad_perp(f1) / c3 + a3 * f1;
    apply_flat_boundary(b3);

    check_laplace(++test_num, "with coefficients constant in z, PETSc 2nd order", *invert,
                  INVERT_AC_GRAD, INVERT_AC_GRAD, a3, c3, d3, b3, f1, mesh->ystart, dump);

    Options* SPT_options = Options::getRoot()->getSection("SPT");
    auto invert_SPT = Laplacian::create(SPT_options);

    check_laplace(++test_num, "with coefficients constant in z, default solver", *invert_SPT,
                  INVERT_AC_GRAD, INVERT_AC_GRAD | INVERT_DC_GRAD, a3, c3, d3, b3, f1,
                  mesh->ystart, dump);

    //////////////////////////////////////////////
    // Test 5: Cosine x-profiles, 2nd order Krylov
    Field3D f5, a5, c5, d5;

    p = 0.623901;
    q = 0.01209489;
    f5.allocate();
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          const BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
          const BoutReal z = BoutReal(jz) / nz;
          //make the gradients zero at both x-boundaries
          f5(jx, jy, jz) = 0. + exp(-(50. * pow(x - p, 2) + 1. - cos(2. * PI * (z - q))))
                           - 50.
                                 * (2. * p * exp(-50. * pow(-p, 2)) * x
                                    + (-p * exp(-50. * pow(-p, 2))
                                       - (1 - p) * exp(-50. * pow(1 - p, 2)))
                                          * pow(x, 2))
                                 * exp(-(1. - cos(2. * PI * (z - q))));
        }
      }
    }
    if (mesh->firstX()) {
      for (int jx = mesh->xstart - 1; jx >= 0; jx--) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            //make the gradients zero at both x-boundaries
            f5(jx, jy, jz) = 0.
                             + exp(-(50. * pow(x - p, 2) + 1. - cos(2. * PI * (z - q))))
                             - 50.
                                   * (2. * p * exp(-50. * pow(-p, 2)) * x
                                      + (-p * exp(-50. * pow(-p, 2))
                                         - (1 - p) * exp(-50. * pow(1 - p, 2)))
                                            * pow(x, 2))
                                   * exp(-(1. - cos(2. * PI * (z - q))));
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int jx = mesh->xend + 1; jx < mesh->LocalNx; jx++) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            //make the gradients zero at both x-boundaries
            f5(jx, jy, jz) = 0.
                             + exp(-(50. * pow(x - p, 2) + 1. - cos(2. * PI * (z - q))))
                             - 50.
                                   * (2. * p * exp(-50. * pow(-p, 2)) * x
                                      + (-p * exp(-50. * pow(-p, 2))
                                         - (1 - p) * exp(-50. * pow(1 - p, 2)))
                                            * pow(x, 2))
                                   * exp(-(1. - cos(2. * PI * (z - q))));
          }
        }
      }
    }

    p = 0.63298589;
    q = 0.889237890;
    d5.allocate();
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
          BoutReal z = BoutReal(jz) / nz;
          d5(jx, jy, jz) = 1. + p * cos(2. * PI * x) * sin(2. * PI * (z - q) * 3.);
        }
      }
    }
    if (mesh->firstX()) {
      for (int jx = mesh->xstart - 1; jx >= 0; jx--) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            d5(jx, jy, jz) = 1. + p * cos(2. * PI * x) * sin(2. * PI * (z - q) * 3.);
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int jx = mesh->xend + 1; jx < mesh->LocalNx; jx++) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            d5(jx, jy, jz) = 1. + p * cos(2. * PI * x) * sin(2. * PI * (z - q) * 3.);
          }
        }
      }
    }

    p = 0.160983834;
    q = 0.73050121087;
    c5.allocate();
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
          BoutReal z = BoutReal(jz) / nz;
          c5(jx, jy, jz) = 1. + p * cos(2. * PI * x * 5) * sin(2. * PI * (z - q) * 2.);
        }
      }
    }
    if (mesh->firstX()) {
      for (int jx = mesh->xstart - 1; jx >= 0; jx--) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            c5(jx, jy, jz) = 1. + p * cos(2. * PI * x * 5) * sin(2. * PI * (z - q) * 2.);
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int jx = mesh->xend + 1; jx < mesh->LocalNx; jx++) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            c5(jx, jy, jz) = 1. + p * cos(2. * PI * x * 5) * sin(2. * PI * (z - q) * 2.);
          }
        }
      }
    }

    p = 0.5378950;
    q = 0.2805870;
    a5.allocate();
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
          BoutReal z = BoutReal(jz) / nz;
          a5(jx, jy, jz) = -1. + p * cos(2. * PI * x * 2.) * sin(2. * PI * (z - q) * 7.);
        }
      }
    }
    if (mesh->firstX()) {
      for (int jx = mesh->xstart - 1; jx >= 0; jx--) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            a5(jx, jy, jz) =
                -1. + p * cos(2. * PI * x * 2.) * sin(2. * PI * (z - q) * 7.);
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int jx = mesh->xend + 1; jx < mesh->LocalNx; jx++) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            BoutReal x = BoutReal(mesh->getGlobalXIndex(jx) - mesh->xstart) / nx;
            BoutReal z = BoutReal(jz) / nz;
            a5(jx, jy, jz) =
                -1. + p * cos(2. * PI * x * 2.) * sin(2. * PI * (z - q) * 7.);
          }
        }
      }
    }

    f5.applyBoundary("neumann");
    mesh->communicate(f5, a5, c5, d5);

    Field3D b5 = d5 * Delp2(f5) + Grad_perp(c5) * Grad_perp(f5) / c5 + a5 * f5;
    apply_flat_boundary(b5);

    check_laplace(++test_num, "different profiles, PETSc 2nd order", *invert, INVERT_AC_GRAD,
                  INVERT_AC_GRAD, a5, c5, d5, b5, f5, mesh->ystart, dump);

    //////////////////////////////////////////////
    // Test 6: Cosine x-profiles, 4th order Krylov

    check_laplace(++test_num, "different profiles, PETSc 4th order", *invert_4th, INVERT_AC_GRAD,
                  INVERT_AC_GRAD, a5, c5, d5, b5, f5, mesh->ystart, dump);

    //////////////////////////////////////////////////////////////////////////////////////
    // Test 7+8: Cosine x-profiles, z-independent coefficients and compare with SPT method

    const Field2D a7 = DC(a5);
    const Field2D c7 = DC(c5);
    const Field2D d7 = DC(d5);
    Field3D b7 = d7 * Delp2(f5) + Grad_perp(c7) * Grad_perp(f5) / c7 + a7 * f5;
    apply_flat_boundary(b7);

    check_laplace(++test_num, "different profiles, with coefficients constant in z, PETSc 2nd order",
        *invert, INVERT_AC_GRAD, INVERT_AC_GRAD, a7, c7, d7, b7, f5, mesh->ystart, dump);

    check_laplace(++test_num,
                  "different profiles, with coefficients constant in z, default solver",
                  *invert_SPT, INVERT_AC_GRAD, INVERT_AC_GRAD | INVERT_DC_GRAD, a7, c7,
                  d7, b7, f5, mesh->ystart, dump);

    // Write and close the output file
    bout::writeDefaultOutputFile(dump);

    MPI_Barrier(BoutComm::get()); // Wait for all processors to write data
  }

  bout::checkForUnusedOptions();

  BoutFinalise();
  return 0;
}

BoutReal max_error_at_ystart(const Field3D& error) {
  const auto* mesh = error.getMesh();
  BoutReal local_max_error = error(mesh->xstart, mesh->ystart, 0);

  for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      if (local_max_error < error(jx, mesh->ystart, jz)) {
        local_max_error = error(jx, mesh->ystart, jz);
      }
    }
  }

  BoutReal max_error = BoutNaN;

  MPI_Allreduce(&local_max_error, &max_error, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());

  return max_error;
}

void apply_flat_boundary(Field3D& bcoef) {
  const Mesh& mesh = *bcoef.getMesh();
  if (mesh.firstX()) {
    for (int jx = mesh.xstart - 1; jx >= 0; jx--) {
      for (int jy = 0; jy < mesh.LocalNy; jy++) {
        for (int jz = 0; jz < mesh.LocalNz; jz++) {
          bcoef(jx, jy, jz) = bcoef(jx + 1, jy, jz);
        }
      }
    }
  }
  if (mesh.lastX()) {
    for (int jx = mesh.xend + 1; jx < mesh.LocalNx; jx++) {
      for (int jy = 0; jy < mesh.LocalNy; jy++) {
        for (int jz = 0; jz < mesh.LocalNz; jz++) {
          bcoef(jx, jy, jz) = bcoef(jx - 1, jy, jz);
        }
      }
    }
  }
}
