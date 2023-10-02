/**************************************************************************
 * Testing Perpendicular Laplacian inversion using PETSc solvers
 *
 **************************************************************************
 * Copyright 2019 J.T. Omotani, B.D. Dudson
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

#include <bout/bout.hxx>
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/options.hxx>

using bout::globals::mesh;

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  auto* coords = mesh->getCoordinates();

  // Solving equations of the form
  // Div(A Grad_perp(f)) + B*f = rhs
  // A*Laplace_perp(f) + Grad_perp(A).Grad_perp(f) + B*f = rhs
  Field2D f;
  Field2D a;
  Field2D b;

  initial_profile("f", f);
  initial_profile("a", a);
  initial_profile("b", b);

  // Apply boundary conditions to f, so the boundary cells match the way boundary
  // conditions will be applied to sol
  f.setBoundary("f");
  f.applyBoundary();

  ////////////////////////////////////////////////////////////////////////////////////////

  Field2D rhs;
  const bool include_y_derivs = Options::root()["laplacexy"]["include_y_derivs"];
  if (include_y_derivs) {
    rhs = a * Laplace_perp(f) + Grad_perp(a) * Grad_perp(f) + b * f;
  } else {
    rhs = a * Delp2(f, CELL_DEFAULT, false) + coords->g11 * DDX(a) * DDX(f) + b * f;
  }

  LaplaceXY laplacexy;
  laplacexy.setCoefs(a, b);

  Field2D solution = laplacexy.solve(rhs, 0.);
  Field2D relative_error = (f - solution) / f;
  Field2D absolute_error = f - solution;
  BoutReal max_error = max(abs(absolute_error), true);

  output.write("Magnitude of maximum absolute error is {}\n", max_error);

  mesh->communicate(solution);
  Field2D rhs_check;
  if (include_y_derivs) {
    rhs_check =
        a * Laplace_perp(solution) + Grad_perp(a) * Grad_perp(solution) + b * solution;
  } else {
    rhs_check = a * Delp2(solution, CELL_DEFAULT, false)
                + coords->g11 * DDX(a) * DDX(solution) + b * solution;
  }

  Options dump;
  dump["a"] = a;
  dump["b"] = b;
  dump["f"] = f;
  dump["sol"] = solution;
  dump["relative_error"] = relative_error;
  dump["absolute_error"] = absolute_error;
  dump["max_error"] = max_error;
  dump["rhs"] = rhs;
  dump["rhs_check"] = rhs_check;
  bout::writeDefaultOutputFile(dump);

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  BoutFinalise();
  return 0;
}
