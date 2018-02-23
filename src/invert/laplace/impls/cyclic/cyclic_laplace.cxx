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

#include <globals.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <bout/sys/timer.hxx>
#include <bout/constants.hxx>
#include <output.hxx>

#include "cyclic_laplace.hxx"

LaplaceCyclic::LaplaceCyclic(Options *opt)
    : Laplacian(opt), Acoef(0.0), Ccoef(1.0), Dcoef(1.0) {
  // Get options

  OPTION(opt, dst, false);

  if(dst) {
    nmode = mesh->LocalNz-2;
  }else
    nmode = maxmode+1; // Number of Z modes. maxmode set in invert_laplace.cxx from options

  // Allocate arrays

  int nsys = nmode;     // Number of tridiagonal systems to solve simulaneously

  xs = mesh->xstart; // Starting X index
  if(mesh->firstX() && !mesh->periodicX){ // Only want to include guard cells at boundaries (unless periodic in x)
	  xs = 0;
  }
  xe = mesh->xend;   // Last X index
  if(mesh->lastX() && !mesh->periodicX){ // Only want to include guard cells at boundaries (unless periodic in x)
	  xe = mesh->LocalNx-1;
  }
  int n = xe - xs + 1;  // Number of X points on this processor,
                        // including boundaries but not guard cells

  a   = Matrix<dcomplex>(nsys, n);
  b   = Matrix<dcomplex>(nsys, n);
  c   = Matrix<dcomplex>(nsys, n);
  xcmplx = Matrix<dcomplex>(nsys, n);
  bcmplx = Matrix<dcomplex>(nsys, n);

  if(dst)
    k1d = Array<dcomplex>(mesh->LocalNz);         // DST has different k space
  else
    k1d = Array<dcomplex>((mesh->LocalNz)/2 + 1); // ZFFT routine expects input of this length

  // Create a cyclic reduction object, operating on dcomplex values
  cr = new CyclicReduce<dcomplex>(mesh->getXcomm(), n);
  cr->setPeriodic(mesh->periodicX);
}

LaplaceCyclic::~LaplaceCyclic() {
  // Delete tridiagonal solver
  delete cr;
}

const FieldPerp LaplaceCyclic::solve(const FieldPerp &rhs, const FieldPerp &x0) {
  Mesh *mesh = rhs.getMesh();
  FieldPerp x(mesh); // Result
  x.allocate();

  Coordinates *coord = mesh->coordinates();

  int jy = rhs.getIndex();  // Get the Y index
  x.setIndex(jy);

  // Get the width of the boundary

  int inbndry = 2, outbndry=2;
  if(global_flags & INVERT_BOTH_BNDRY_ONE) {
    inbndry = outbndry = 1;
  }
  if(inner_boundary_flags & INVERT_BNDRY_ONE)
    inbndry = 1;
  if(outer_boundary_flags & INVERT_BNDRY_ONE)
    outbndry = 1;

  if(dst) {
    // Loop over X indices, including boundaries but not guard cells. (unless periodic in x)
    for(int ix=xs; ix <= xe; ix++) {
      // Take DST in Z direction and put result in k1d

      if(((ix < inbndry) && (inner_boundary_flags & INVERT_SET) && mesh->firstX()) ||
         ((xe-ix < outbndry) && (outer_boundary_flags & INVERT_SET) && mesh->lastX())) {
        // Use the values in x0 in the boundary
        DST(x0[ix]+1, mesh->LocalNz-2 , std::begin(k1d));
      }else {
        DST(rhs[ix]+1, mesh->LocalNz-2 , std::begin(k1d));
      }

      // Copy into array, transposing so kz is first index
      for(int kz = 0; kz < nmode; kz++)
        bcmplx(kz, ix-xs) = k1d[kz];
    }

    // Get elements of the tridiagonal matrix
    // including boundary conditions
    for(int kz = 0; kz < nmode; kz++) {
      BoutReal zlen = coord->dz*(mesh->LocalNz-3);
      BoutReal kwave=kz*2.0*PI/(2.*zlen); // wave number is 1/[rad]; DST has extra 2.

      tridagMatrix(&a(kz,0), &b(kz,0), &c(kz,0), &bcmplx(kz,0), jy,
                   kz,    // wave number index
                   kwave, // kwave (inverse wave length)
                   global_flags, inner_boundary_flags, outer_boundary_flags, &Acoef,
                   &Ccoef, &Dcoef,
                   false); // Don't include guard cells in arrays
    }

    // Solve tridiagonal systems

    cr->setCoefs(nmode, a, b, c);
    cr->solve(nmode, bcmplx, xcmplx);

    // FFT back to real space
    for(int ix=xs; ix <= xe; ix++) {
      for(int kz = 0; kz < nmode; kz++)
        k1d[kz] = xcmplx(kz, ix-xs);

      for(int kz=nmode;kz<(mesh->LocalNz);kz++)
        k1d[kz] = 0.0; // Filtering out all higher harmonics

      DST_rev(std::begin(k1d), mesh->LocalNz-2, x[ix]+1);

      x(ix, 0) = -x(ix, 2);
      x(ix, mesh->LocalNz-1) = -x(ix, mesh->LocalNz-3);
    }
  }else {
    // Loop over X indices, including boundaries but not guard cells (unless periodic in x)
    for(int ix=xs; ix <= xe; ix++) {
      // Take FFT in Z direction, apply shift, and put result in k1d

    if(((ix < inbndry) && (inner_boundary_flags & INVERT_SET) && mesh->firstX()) ||
       ((xe-ix < outbndry) && (outer_boundary_flags & INVERT_SET) && mesh->lastX())) {
        // Use the values in x0 in the boundary
      rfft(x0[ix], mesh->LocalNz, std::begin(k1d));
    }else {
      rfft(rhs[ix], mesh->LocalNz, std::begin(k1d));
      }

      // Copy into array, transposing so kz is first index
      for(int kz = 0; kz < nmode; kz++)
        bcmplx(kz, ix-xs) = k1d[kz];
    }

    // Get elements of the tridiagonal matrix
    // including boundary conditions
    for(int kz = 0; kz < nmode; kz++) {
      BoutReal kwave=kz*2.0*PI/(coord->zlength()); // wave number is 1/[rad]
      tridagMatrix(&a(kz,0), &b(kz,0), &c(kz,0), &bcmplx(kz,0), jy,
                   kz,    // True for the component constant (DC) in Z
                   kwave, // Z wave number
                   global_flags, inner_boundary_flags, outer_boundary_flags, &Acoef,
                   &Ccoef, &Dcoef,
                   false); // Don't include guard cells in arrays
    }

    // Solve tridiagonal systems

    cr->setCoefs(nmode, a, b, c);
    cr->solve(nmode, bcmplx, xcmplx);

    // FFT back to real space
    for(int ix=xs; ix <= xe; ix++) {
      for(int kz = 0; kz < nmode; kz++)
        k1d[kz] = xcmplx(kz, ix-xs);

      for(int kz=nmode;kz<(mesh->LocalNz)/2 + 1;kz++)
        k1d[kz] = 0.0; // Filtering out all higher harmonics

      irfft(std::begin(k1d), mesh->LocalNz, x[ix]);
    }
  }
  return x;
}
