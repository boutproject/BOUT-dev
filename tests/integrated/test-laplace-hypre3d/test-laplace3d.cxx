/**************************************************************************
 * Testing 3d inversion of perpendicular Laplacian
 *
 **************************************************************************
 * Copyright 2019 J.T. Omotani, C. MacMackin
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

#include "bout.hxx"
#include "initialprofiles.hxx"
#include "invert_laplace.hxx"

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  ///////////////////////////////////////////////////////////////////////////////////////
  // Initialise input
  ///////////////////////////////////////////////////////////////////////////////////////
  Field3D f, rhs;
  initial_profile("rhs", rhs);

  // initial profile of f only used to set boundary values
  initial_profile("f", f);

  auto* mesh = f.getMesh();

  Field3D guess = zeroFrom(rhs);

  // Copy boundary values into boundary cells
  for (auto it = mesh->iterateBndryLowerY(); !it.isDone(); it.next()) {
    int x = it.ind;
    int y = mesh->ystart - 1;
    if (x == mesh->xstart) {
      for (int z = mesh->zstart; z <= mesh->zend; z++) {
	guess(x-1, y, z) = 0.5*(f(x-1, y - 1, z) + f(x-1, y, z));
      }
    }
    for (int z = mesh->zstart; z <= mesh->zend; z++) {
      guess(x, y, z) = 0.5*(f(x, y, z) + f(x, y + 1, z));
    }
    if (x == mesh->xend) {
      for (int z = mesh->zstart; z <= mesh->zend; z++) {
	guess(x+1, y, z) = 0.5*(f(x+1, y - 1, z) + f(x+1, y, z));
      }
    }
  }
  for (auto it = mesh->iterateBndryUpperY(); !it.isDone(); it.next()) {
    int x = it.ind;
    int y = mesh->yend + 1;
    if (x == mesh->xstart) {
      for (int z = mesh->zstart; z <= mesh->zend; z++) {
	guess(x-1, y, z) = 0.5*(f(x-1, y - 1, z) + f(x-1, y, z));
      }
    }
    for (int z = mesh->zstart; z <= mesh->zend; z++) {
      guess(x, y, z) = 0.5*(f(x, y - 1, z) + f(x, y, z));
    }
    if (x == mesh->xend) {
      for (int z = mesh->zstart; z <= mesh->zend; z++) {
	guess(x+1, y, z) = 0.5*(f(x+1, y - 1, z) + f(x+1, y, z));
      }
    }
  }
  if (mesh->firstX()) {
    int x = mesh->xstart - 1;
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      for (int z = mesh->zstart; z <= mesh->zend; z++) {
        guess(x, y, z) = 0.5*(f(x, y, z) + f(x + 1, y, z));
      }
    }
  }
  if (mesh->lastX()) {
    int x = mesh->xend + 1;
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      for (int z = mesh->zstart; z <= mesh->zend; z++) {
        guess(x, y, z) = 0.5*(f(x - 1, y, z) + f(x, y, z));
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  // Set up Laplace solver
  ///////////////////////////////////////////////////////////////////////////////////////
  auto laplace_solver = Laplacian::create();
  Field3D A, C1, C2, D;
  initial_profile("A", A);
  initial_profile("C1", C1);
  initial_profile("C2", C2);
  mesh->communicate(C2);
  C2.setBoundary("C2");
  C2.applyParallelBoundary();
  initial_profile("D", D);
  laplace_solver->setCoefA(A);
  laplace_solver->setCoefC1(C1);
  laplace_solver->setCoefC2(C2);
  laplace_solver->setCoefD(D);

  ///////////////////////////////////////////////////////////////////////////////////////
  // Solve
  ///////////////////////////////////////////////////////////////////////////////////////
  f = laplace_solver->solve(rhs, guess);

  ///////////////////////////////////////////////////////////////////////////////////////
  // Calculate error
  ///////////////////////////////////////////////////////////////////////////////////////
  Field3D rhs_check = D*Laplace_perp(f) + Grad_perp(C2)*Grad_perp(f)/C1 + A*f;
  Field3D error = rhs_check - rhs;
  BoutReal error_max = max(abs(error), true);

  SAVE_ONCE(f, rhs, rhs_check, error, error_max);
  bout::globals::dump.write();

  laplace_solver.reset(nullptr);
  BoutFinalise();

  return 0;
}
