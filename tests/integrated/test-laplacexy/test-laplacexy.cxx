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

#include <bout.hxx>
#include <bout/constants.hxx>
#include <bout/invert/laplacexy.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <options.hxx>

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  auto coords = mesh->getCoordinates();

  auto& opt = Options::root();
  
  LaplaceXY laplacexy;

  bool include_y_derivs = opt["laplacexy"]["include_y_derivs"];

  // Solving equations of the form
  // Div(A Grad_perp(f)) + B*f = rhs
  // A*Laplace_perp(f) + Grad_perp(A).Grad_perp(f) + B*f = rhs
  Field2D f, a, b, sol;
  Field2D error, absolute_error; //Absolute value of relative error: abs((f - sol)/f)
  BoutReal max_error; //Output of test

  initial_profile("f", f);
  initial_profile("a", a);
  initial_profile("b", b);

  // Apply boundary conditions to f, so the boundary cells match the way boundary
  // conditions will be applied to sol
  f.setBoundary("f");
  f.applyBoundary();

  ////////////////////////////////////////////////////////////////////////////////////////

  Field2D rhs;
  if (include_y_derivs) {
    rhs = a*Laplace_perp(f) + Grad_perp(a)*Grad_perp(f) + b*f;
  } else {
    rhs = a*Delp2(f, CELL_DEFAULT, false) + coords->g11*DDX(a)*DDX(f) + b*f;
  }

  laplacexy.setCoefs(a, b);
  
  sol = laplacexy.solve(rhs, 0.);
  error = (f - sol)/f;
  absolute_error = f - sol;
  max_error = max(abs(absolute_error), true);

  output<<"Magnitude of maximum absolute error is "<<max_error<<endl;

  dump.add(a, "a");
  dump.add(b, "b");
  dump.add(f, "f");
  dump.add(sol, "sol");
  dump.add(error, "error");
  dump.add(absolute_error, "absolute_error");
  dump.add(max_error, "max_error");

  dump.write();
  dump.close();

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  BoutFinalise();
  return 0;
}
