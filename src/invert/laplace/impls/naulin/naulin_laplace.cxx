/*!
 * \file cyclic.cxx
 *
 * \brief FFT + Tridiagonal solver in serial or parallel
 *
 * Not particularly optimised: Each y slice is solved sequentially
 *
 * CHANGELOG
 * =========
 *
 * Jan 2014: Brendan Shanahan <bws502@york.ac.uk>
 *         * Added DST option
 *
 **************************************************************************
 * Copyright 2013 B.D.Dudson
 *
 * Contact: Ben Dudson, benjamin.dudson@york.ac.uk
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
 */

#include <boutexception.hxx>
#include <bout/mesh.hxx>
#include <bout/coordinates.hxx>
#include <bout/sys/timer.hxx>
#include <derivs.hxx>
#include <globals.hxx>
#include <output.hxx>

#include "naulin_laplace.hxx"

LaplaceNaulin::LaplaceNaulin(Options *opt)
    : Laplacian(opt), C1coef(1.0), C2coef(0.0), Dcoef(1.0),
      delp2solver(nullptr), naulinsolver_mean_its(0.), ncalls(0) {

  // Get options
  OPTION(opt, rtol, 1.e-7);
  OPTION(opt, atol, 1.e-20);
  OPTION(opt, maxits, 100);
  delp2solver = create(opt->getSection("delp2solver"));
  // Use same flags for FFT solver as for NaulinSolver
  delp2solver->setGlobalFlags(global_flags);
  delp2solver->setInnerBoundaryFlags(inner_boundary_flags);
  delp2solver->setOuterBoundaryFlags(outer_boundary_flags);

  static bool first = true;
  if (first) {
    SAVE_REPEAT(naulinsolver_mean_its);
    first = false;
  }
}

LaplaceNaulin::~LaplaceNaulin() {
  delete delp2solver;
}

const Field3D LaplaceNaulin::solve(const Field3D &rhs, const Field3D &x0) {
  // Rearrange equation so first term is just Delp2(x):
  //   D*Delp2(x) + 1/C1*Grad_perp(C2).Grad_perp(phi) = rhs
  //   -> Delp2(x) + 1/(C1*D)*Grad_perp(C2).Grad_perp(phi) = rhs/D

  Timer timer("invert"); ///< Start timer

  CELL_LOC location = rhs.getLocation();
  ASSERT1(Dcoef.getLocation() == location);
  ASSERT1(C1coef.getLocation() == location);
  ASSERT1(C2coef.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  Mesh *mesh = rhs.getMesh();
  Coordinates *coords = mesh->coordinates();
  Field3D x(x0); // Result

  Field3D rhsOverD = rhs/Dcoef;
  Field3D ddx_c = DDX(C2coef, location, DIFF_C2);
  Field3D ddz_c = DDZ(C2coef, location, DIFF_FFT);
  Field3D oneOverC1coefTimesDcoef = 1./C1coef/Dcoef;

  BoutReal error_rel = 1e20, error_abs=1e20;
  int count = 0;
  while (error_rel>rtol && error_abs>atol) {

    copy_x_boundaries(x, x0, mesh);
    Field3D ddx_x = DDX(x, location, DIFF_C2);
    Field3D ddz_x = DDZ(x, location, DIFF_FFT);
    Field3D b = rhsOverD - (coords->g11*ddx_c*ddx_x + coords->g33*ddz_c*ddz_x + coords->g13*(ddx_c*ddz_x + ddz_c*ddx_x))*oneOverC1coefTimesDcoef;
    // NB need to pass x in case boundary flags require 'x0', even if
    // delp2solver is not iterative and does not use an initial guess
    Field3D xnew = delp2solver->solve(b, x);
    Field3D difference = xnew-x;
    error_abs = max(abs(difference, RGN_NOBNDRY), true, RGN_NOBNDRY);
    error_rel = error_abs / sqrt(mean(SQ(x), true, RGN_NOBNDRY)); // use sqrt(mean(SQ)) to make sure we do not divide by zero at a point
    x = xnew;
    mesh->communicate(x);
    //output<<"NaulinSolver: "<<count<<" "<<error_rel<<" "<<error_abs<<endl;

    count++;
    if (count>maxits)
      throw BoutException("LaplaceNaulin error: Took more than maxits=%i iterations to converge.", maxits);
  }

  ncalls++;
  naulinsolver_mean_its = (naulinsolver_mean_its*BoutReal(ncalls-1) + BoutReal(count))/BoutReal(ncalls);

  return x;
}

void LaplaceNaulin::copy_x_boundaries(Field3D &x, const Field3D &x0, Mesh *mesh) {
  if (mesh->firstX()) {
    for (int i=mesh->xstart-1; i>=0; i--)
      for (int j=mesh->ystart; j<=mesh->yend; j++)
        for (int k=0; k<mesh->LocalNz; k++)
          x(i, j, k) = x0(i, j, k);
  }
  if (mesh->lastX()) {
    for (int i=mesh->xend+1; i<mesh->LocalNx; i++)
      for (int j=mesh->ystart; j<=mesh->yend; j++)
        for (int k=0; k<mesh->LocalNz; k++)
          x(i, j, k) = x0(i, j, k);
  }
}
