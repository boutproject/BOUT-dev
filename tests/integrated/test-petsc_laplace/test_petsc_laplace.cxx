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

#include "bout/bout.hxx" // NOLINT
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/build_config.hxx"
#include "bout/constants.hxx"
#include "bout/difops.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/field_factory.hxx"
#include "bout/invert_laplace.hxx"
#include "bout/options.hxx"
#include "bout/options_io.hxx"
#include "bout/output.hxx"
#include "bout/traits.hxx"
#include "bout/vecops.hxx"

#include "fmt/format.h"
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <string_view>

namespace {
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

template <class T>
Field3D forward_laplace(const Field3D& field, const T& acoef, const T& ccoef,
                        const T& dcoef) {
  auto bcoef =
      dcoef * Delp2(field) + Grad_perp(ccoef) * Grad_perp(field) / ccoef + acoef * field;
  apply_flat_boundary(bcoef);
  return bcoef;
}

Field3D generate_f1(const Mesh& mesh);

Field3D generate_d1(Mesh& mesh) {
  Options options{{"p_d1", 0.512547}, {"q_d1", 0.30908712}};

  return FieldFactory::get()->create3D(
      "1. + (0.2 * exp(-50. * (x - p_d1)^2 / 4.) * sin((z - 2. * pi * q_d1) * 3.))",
      &options, &mesh);
}

Field3D generate_c1(Mesh& mesh) {
  Options options{{"p_c1", 0.18439023}, {"q_c1", 0.401089473}};

  return FieldFactory::get()->create3D(
      "1. + (0.15 * exp(-50. * (x - p_c1)^2 * 2.) * sin((z - 2. * pi * q_c1) * 2.))",
      &options, &mesh);
}

Field3D generate_a1(Mesh& mesh) {
  Options options{{"p_a1", 0.612547}, {"q_a1", 0.30908712}};

  return FieldFactory::get()->create3D(
      "-1. + (0.1 * exp(-50. * (x - p_a1)^2 * 2.5) * sin((z - 2. * pi * q_a1) * 7.))",
      &options, &mesh);
}

Field3D generate_f5(Mesh& mesh) {
  Options options{{"p_f5", 0.623901}, {"q_f5", 0.01209489}};

  auto result = FieldFactory::get()->create3D(
      "exp(-((50. * (x - p_f5)^2) + 1. - cos((z - 2. * pi * q_f5))))"
      "- 50 * (2. * p_f5 * exp(-50. * (-p_f5)^2) * x"
      "         + (-p_f5 * exp(-50. * (-p_f5)^2) - (1 - p_f5) * exp(-50. * (1 - "
      "p_f5)^2)) * x^2)"
      "     * exp(-(1. - cos(z - 2. * pi * q_f5)))",
      &options, &mesh);
  result.applyBoundary("neumann");
  checkData(result);
  return result;
}

Field3D generate_d5(Mesh& mesh) {
  Options options{{"p_d5", 0.63298589}, {"q_d5", 0.889237890}};

  return FieldFactory::get()->create3D(
      "1. + (p_d5 * cos(2. * pi * x) * sin((z - 2 * pi * q_d5) * 3.))", &options, &mesh);
}

Field3D generate_c5(Mesh& mesh) {
  Options options{{"p_c5", 0.160983834}, {"q_c5", 0.73050121087}};

  return FieldFactory::get()->create3D(
      "1. + (p_c5 * cos(2. * pi * x * 5) * sin((z - 2 * pi * q_c5) * 2.))", &options,
      &mesh);
}

Field3D generate_a5(Mesh& mesh) {
  Options options{{"p_a5", 0.5378950}, {"q_a5", 0.2805870}};

  return FieldFactory::get()->create3D(
      "-1. + (p_a5 * cos(2. * pi * x * 2.) * sin((z - 2. * pi * q_a5) * 7.))", &options,
      &mesh);
}
} // namespace

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);
  {
    // Not be used for 3D metrics
    Options::root()["laplace"].setConditionallyUsed();

    // For 3D metrics, we need to use one of PETSc's preconditioners
    // and mumps seems to be the best? Or at least HYPRE doesn't work
    // on this test with more than one processor.
#ifdef PETSC_HAVE_MUMPS
    constexpr auto petsc_has_mumps = true;
#else
    constexpr auto petsc_has_mumps = false;
#endif

    // Preconditioner to use for 3D metrics
    constexpr auto petsc_pc = petsc_has_mumps ? "mumps" : "lu";

    auto& options_2nd = Options::root()["petsc2nd"];
    if constexpr (bout::build::use_metric_3d) {
      options_2nd["pctype"] = petsc_pc;
      options_2nd["rightprec"].setConditionallyUsed();
      options_2nd["precon"].setConditionallyUsed();
    }

    auto invert = Laplacian::create(&options_2nd);

    auto& options_4th = Options::root()["petsc4th"];
    if constexpr (bout::build::use_metric_3d) {
      options_4th["pctype"] = petsc_pc;
      options_4th["rightprec"].setConditionallyUsed();
      options_4th["precon"].setConditionallyUsed();
    }
    auto invert_4th = Laplacian::create(&options_4th);

    Options dump;

    // Solving equations of the form d*Delp2(f) + 1/c*Grad_perp(c).Grad_perp(f) + a*f = b for various f, a, c, d
    using bout::globals::mesh;

    Field3D loop_x;
    Field3D loop_z;
    loop_x.allocate();
    loop_z.allocate();
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      const BoutReal x = mesh->GlobalX(jx);
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        for (int jz = mesh->zstart; jz <= mesh->zend; jz++) {
          const BoutReal z = mesh->GlobalZ(jz);
          loop_x(jx, jy, jz) = x;
          loop_z(jx, jy, jz) = z;
        }
      }
    }
    dump["loop_x"] = loop_x;
    dump["loop_z"] = loop_z;
    dump["exp_x"] = FieldFactory::get()->create3D("x");
    dump["exp_z"] = FieldFactory::get()->create3D("z");

    // Only Neumann x-boundary conditions are implemented so far, so test functions should be Neumann in x and periodic in z.
    // Use Field3D's, but solver only works on FieldPerp slices, so only use 1 y-point

    /////////////////////////////////////////////////////
    // Test 1: Gaussian x-profiles, 2nd order Krylov
    Field3D f_1 = generate_f1(*mesh);
    Field3D a_1 = generate_a1(*mesh);
    Field3D c_1 = generate_c1(*mesh);
    Field3D d_1 = generate_d1(*mesh);

    mesh->communicate(f_1, a_1, c_1, d_1);

    const Field3D b_1 = forward_laplace(f_1, a_1, c_1, d_1);

    check_laplace(1, "PETSc 2nd order", *invert, INVERT_AC_GRAD, INVERT_AC_GRAD, a_1, c_1,
                  d_1, b_1, f_1, mesh->ystart, dump);

    /////////////////////////////////////////////////
    // Test 2: Gaussian x-profiles, 4th order Krylov

    check_laplace(2, "PETSc 4th order", *invert_4th, INVERT_AC_GRAD, INVERT_AC_GRAD, a_1,
                  c_1, d_1, b_1, f_1, mesh->ystart, dump);

    ////////////////////////////////////////////////////////////////////////////////////////
    // Test 3+4: Gaussian x-profiles, z-independent coefficients

    const Field2D a_3 = DC(a_1);
    const Field2D c_3 = DC(c_1);
    const Field2D d_3 = DC(d_1);
    const Field3D b_3 = forward_laplace(f_1, a_3, c_3, d_3);

    check_laplace(3, "with coefficients constant in z, PETSc 2nd order", *invert,
                  INVERT_AC_GRAD, INVERT_AC_GRAD, a_3, c_3, d_3, b_3, f_1, mesh->ystart,
                  dump);

    //////////////////////////////////////////////
    // Test 5: Cosine x-profiles, 2nd order Krylov
    Field3D f_5 = generate_f5(*mesh);
    Field3D a_5 = generate_a5(*mesh);
    Field3D c_5 = generate_c5(*mesh);
    Field3D d_5 = generate_d5(*mesh);

    mesh->communicate(f_5, a_5, c_5, d_5);

    const Field3D b_5 = forward_laplace(f_5, a_5, c_5, d_5);

    check_laplace(5, "different profiles, PETSc 2nd order", *invert, INVERT_AC_GRAD,
                  INVERT_AC_GRAD, a_5, c_5, d_5, b_5, f_5, mesh->ystart, dump);

    //////////////////////////////////////////////
    // Test 6: Cosine x-profiles, 4th order Krylov

    check_laplace(6, "different profiles, PETSc 4th order", *invert_4th, INVERT_AC_GRAD,
                  INVERT_AC_GRAD, a_5, c_5, d_5, b_5, f_5, mesh->ystart, dump);

    //////////////////////////////////////////////////////////////////////////////////////
    // Test 7+8: Cosine x-profiles, z-independent coefficients and compare with SPT method

    const Field2D a_7 = DC(a_5);
    const Field2D c_7 = DC(c_5);
    const Field2D d_7 = DC(d_5);
    const Field3D b_7 = forward_laplace(f_5, a_7, c_7, d_7);

    check_laplace(7,
                  "different profiles, with coefficients constant in z, PETSc 2nd order",
                  *invert, INVERT_AC_GRAD, INVERT_AC_GRAD, a_7, c_7, d_7, b_7, f_5,
                  mesh->ystart, dump);

    // Write and close the output file
    bout::writeDefaultOutputFile(dump);

    MPI_Barrier(BoutComm::get()); // Wait for all processors to write data
  }

  bout::checkForUnusedOptions();

  BoutFinalise();
  return 0;
}

namespace {
BoutReal max_error_at_ystart(const Field3D& error) {
  const auto* mesh = error.getMesh();
  BoutReal local_max_error = 0.0;

  for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
    for (int jz = mesh->zstart; jz <= mesh->zend; jz++) {
      local_max_error = std::max(local_max_error, error(jx, mesh->ystart, jz));
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
        for (int jz = mesh.zstart; jz <= mesh.zend; jz++) {
          bcoef(jx, jy, jz) = bcoef(jx + 1, jy, jz);
        }
      }
    }
  }
  if (mesh.lastX()) {
    for (int jx = mesh.xend + 1; jx < mesh.LocalNx; jx++) {
      for (int jy = 0; jy < mesh.LocalNy; jy++) {
        for (int jz = mesh.zstart; jz <= mesh.zend; jz++) {
          bcoef(jx, jy, jz) = bcoef(jx - 1, jy, jz);
        }
      }
    }
  }
}

Field3D generate_f1(const Mesh& mesh) {
  constexpr BoutReal p = 0.39503274; // NOLINT
  constexpr BoutReal q = 0.20974396; // NOLINT

  Field3D result;
  result.allocate();
  for (int jx = mesh.xstart; jx <= mesh.xend; jx++) {
    const BoutReal x = mesh.GlobalX(jx);
    for (int jy = 0; jy < mesh.LocalNy; jy++) {
      for (int jz = mesh.zstart; jz <= mesh.zend; jz++) {
        const BoutReal z = mesh.GlobalZ(jz);
        //make the gradients zero at both x-boundaries
        result(jx, jy, jz) =
            0. + exp(-((100. * pow(x - p, 2)) + 1. - cos(2. * PI * (z - q))))
            - 50.
                  * (2. * p * exp(-100. * pow(-p, 2)) * x
                     + (-p * exp(-100. * pow(-p, 2))
                        - (1 - p) * exp(-100. * pow(1 - p, 2)))
                           * pow(x, 2))
                  * exp(-(1. - cos(2. * PI * (z - q))));
      }
    }
  }
  if (mesh.firstX()) {
    for (int jx = mesh.xstart - 1; jx >= 0; jx--) {
      const BoutReal x = mesh.GlobalX(jx);
      for (int jy = 0; jy < mesh.LocalNy; jy++) {
        for (int jz = mesh.zstart; jz <= mesh.zend; jz++) {
          const BoutReal z = mesh.GlobalZ(jz);
          //make the gradients zero at both x-boundaries
          result(jx, jy, jz) =
              0. + exp(-((60. * pow(x - p, 2)) + 1. - cos(2. * PI * (z - q))))
              - 50.
                    * (2. * p * exp(-60. * pow(-p, 2)) * x
                       + (-p * exp(-60. * pow(-p, 2))
                          - (1 - p) * exp(-60. * pow(1 - p, 2)))
                             * pow(x, 2))
                    * exp(-(1. - cos(2. * PI * (z - q))));
        }
      }
    }
  }
  if (mesh.lastX()) {
    for (int jx = mesh.xend + 1; jx < mesh.LocalNx; jx++) {
      const BoutReal x = mesh.GlobalX(jx);
      for (int jy = 0; jy < mesh.LocalNy; jy++) {
        for (int jz = mesh.zstart; jz <= mesh.zend; jz++) {
          const BoutReal z = mesh.GlobalZ(jz);
          //make the gradients zero at both x-boundaries
          result(jx, jy, jz) =
              0. + exp(-((60. * pow(x - p, 2)) + 1. - cos(2. * PI * (z - q))))
              - 50.
                    * (2. * p * exp(-60. * pow(-p, 2)) * x
                       + (-p * exp(-60. * pow(-p, 2))
                          - (1 - p) * exp(-60. * pow(1 - p, 2)))
                             * pow(x, 2))
                    * exp(-(1. - cos(2. * PI * (z - q))));
        }
      }
    }
  }

  checkData(result);
  result.applyBoundary("neumann");
  return result;
}
} // namespace
