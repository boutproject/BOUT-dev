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
#include <bout/mesh.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <bout/sys/timer.hxx>
#include <bout/constants.hxx>
#include <output.hxx>

#include "cyclic_laplace.hxx"

LaplaceCyclic::LaplaceCyclic(Options *opt, const CELL_LOC loc, Mesh *mesh_in)
    : Laplacian(opt, loc, mesh_in), Acoef(0.0), C1coef(1.0), C2coef(1.0), Dcoef(1.0) {
  Acoef.setLocation(location);
  C1coef.setLocation(location);
  C2coef.setLocation(location);
  Dcoef.setLocation(location);

  // Get options

  OPTION(opt, dst, false);

  if(dst) {
    nmode = localmesh->LocalNz-2;
  }else
    nmode = maxmode+1; // Number of Z modes. maxmode set in invert_laplace.cxx from options

  // Note nmode == nsys of cyclic_reduction

  // Allocate arrays

  xs = localmesh->xstart; // Starting X index
  if(localmesh->firstX() && !localmesh->periodicX){ // Only want to include guard cells at boundaries (unless periodic in x)
	  xs = 0;
  }
  xe = localmesh->xend;   // Last X index
  if(localmesh->lastX() && !localmesh->periodicX){ // Only want to include guard cells at boundaries (unless periodic in x)
	  xe = localmesh->LocalNx-1;
  }
  int n = xe - xs + 1;  // Number of X points on this processor,
                        // including boundaries but not guard cells

  a.reallocate(nmode, n);
  b.reallocate(nmode, n);
  c.reallocate(nmode, n);
  xcmplx.reallocate(nmode, n);
  bcmplx.reallocate(nmode, n);

  // Create a cyclic reduction object, operating on dcomplex values
  cr = new CyclicReduce<dcomplex>(localmesh->getXcomm(), n);
  cr->setPeriodic(localmesh->periodicX);
}

LaplaceCyclic::~LaplaceCyclic() {
  // Delete tridiagonal solver
  delete cr;
}

FieldPerp LaplaceCyclic::solve(const FieldPerp& rhs, const FieldPerp& x0) {
  ASSERT1(localmesh == rhs.getMesh() && localmesh == x0.getMesh());
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  FieldPerp x{emptyFrom(rhs)}; // Result

  int jy = rhs.getIndex();  // Get the Y index
  x.setIndex(jy);

  // Get the width of the boundary

  // If the flags to assign that only one guard cell should be used is set
  int inbndry = localmesh->xstart, outbndry=localmesh->xstart;
  if((global_flags & INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2))  {
    inbndry = outbndry = 1;
  }
  if(inner_boundary_flags & INVERT_BNDRY_ONE)
    inbndry = 1;
  if(outer_boundary_flags & INVERT_BNDRY_ONE)
    outbndry = 1;

  if(dst) {
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d =
          Array<dcomplex>(localmesh->LocalNz); // ZFFT routine expects input of this length

      // Loop over X indices, including boundaries but not guard cells. (unless periodic
      // in x)
      BOUT_OMP(for)
      for (int ix = xs; ix <= xe; ix++) {
        // Take DST in Z direction and put result in k1d

        if (((ix < inbndry) && (inner_boundary_flags & INVERT_SET) && localmesh->firstX()) ||
            ((localmesh->LocalNx - ix - 1 < outbndry) && (outer_boundary_flags & INVERT_SET) &&
             localmesh->lastX())) {
          // Use the values in x0 in the boundary
          DST(x0[ix] + 1, localmesh->LocalNz - 2, std::begin(k1d));
        } else {
          DST(rhs[ix] + 1, localmesh->LocalNz - 2, std::begin(k1d));
        }

        // Copy into array, transposing so kz is first index
        for (int kz = 0; kz < nmode; kz++)
          bcmplx(kz, ix - xs) = k1d[kz];
      }

      // Get elements of the tridiagonal matrix
      // including boundary conditions
      BOUT_OMP(for nowait)
      for (int kz = 0; kz < nmode; kz++) {
        BoutReal zlen = coords->dz * (localmesh->LocalNz - 3);
        BoutReal kwave =
            kz * 2.0 * PI / (2. * zlen); // wave number is 1/[rad]; DST has extra 2.

        tridagMatrix(&a(kz, 0), &b(kz, 0), &c(kz, 0), &bcmplx(kz, 0), jy,
                     kz,    // wave number index
                     kwave, // kwave (inverse wave length)
                     global_flags, inner_boundary_flags, outer_boundary_flags, &Acoef,
                     &C1coef, &C2coef, &Dcoef,
                     false); // Don't include guard cells in arrays
      }
    }

    // Solve tridiagonal systems
    cr->setCoefs(a, b, c);
    cr->solve(bcmplx, xcmplx);

    // FFT back to real space
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d =
          Array<dcomplex>(localmesh->LocalNz); // ZFFT routine expects input of this length

      BOUT_OMP(for nowait)
      for (int ix = xs; ix <= xe; ix++) {
        for (int kz = 0; kz < nmode; kz++)
          k1d[kz] = xcmplx(kz, ix - xs);

        for (int kz = nmode; kz < (localmesh->LocalNz); kz++)
          k1d[kz] = 0.0; // Filtering out all higher harmonics

        DST_rev(std::begin(k1d), localmesh->LocalNz - 2, x[ix] + 1);

        x(ix, 0) = -x(ix, 2);
        x(ix, localmesh->LocalNz - 1) = -x(ix, localmesh->LocalNz - 3);
      }
    }
  }else {
    BOUT_OMP(parallel)
    {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>((localmesh->LocalNz) / 2 +
                                 1); // ZFFT routine expects input of this length

      // Loop over X indices, including boundaries but not guard cells (unless periodic in
      // x)
      BOUT_OMP(for)
      for (int ix = xs; ix <= xe; ix++) {
        // Take FFT in Z direction, apply shift, and put result in k1d

        if (((ix < inbndry) && (inner_boundary_flags & INVERT_SET) && localmesh->firstX()) ||
            ((localmesh->LocalNx - ix - 1 < outbndry) && (outer_boundary_flags & INVERT_SET) &&
             localmesh->lastX())) {
          // Use the values in x0 in the boundary
          rfft(x0[ix], localmesh->LocalNz, std::begin(k1d));
        } else {
          rfft(rhs[ix], localmesh->LocalNz, std::begin(k1d));
        }

        // Copy into array, transposing so kz is first index
        for (int kz = 0; kz < nmode; kz++)
          bcmplx(kz, ix - xs) = k1d[kz];
      }

      // Get elements of the tridiagonal matrix
      // including boundary conditions
      BOUT_OMP(for nowait)
      for (int kz = 0; kz < nmode; kz++) {
        BoutReal kwave = kz * 2.0 * PI / (coords->zlength()); // wave number is 1/[rad]
        tridagMatrix(&a(kz, 0), &b(kz, 0), &c(kz, 0), &bcmplx(kz, 0), jy,
                     kz,    // True for the component constant (DC) in Z
                     kwave, // Z wave number
                     global_flags, inner_boundary_flags, outer_boundary_flags, &Acoef,
                     &C1coef, &C2coef, &Dcoef,
                     false); // Don't include guard cells in arrays
      }
    }

    // Solve tridiagonal systems
    cr->setCoefs(a, b, c);
    cr->solve(bcmplx, xcmplx);

    // FFT back to real space
    BOUT_OMP(parallel)
    {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>((localmesh->LocalNz) / 2 +
                                 1); // ZFFT routine expects input of this length

      bool zero_DC = false;
      if(global_flags & INVERT_ZERO_DC) {
        // No DC component
        zero_DC = true;
      }

      BOUT_OMP(for nowait)
      for (int ix = xs; ix <= xe; ix++) {
        if (zero_DC) {
          k1d[0] = 0.;
        }

        for (int kz = zero_DC; kz < nmode; kz++)
          k1d[kz] = xcmplx(kz, ix - xs);

        for (int kz = nmode; kz < (localmesh->LocalNz) / 2 + 1; kz++)
          k1d[kz] = 0.0; // Filtering out all higher harmonics

        irfft(std::begin(k1d), localmesh->LocalNz, x[ix]);
      }
    }
  }
  return x;
}

Field3D LaplaceCyclic::solve(const Field3D& rhs, const Field3D& x0) {
  TRACE("LaplaceCyclic::solve(Field3D, Field3D)");

  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  ASSERT1(localmesh == rhs.getMesh() && localmesh == x0.getMesh());

  Timer timer("invert");

  Field3D x{emptyFrom(rhs)}; // Result

  // Get the width of the boundary

  // If the flags to assign that only one guard cell should be used is set
  int inbndry = localmesh->xstart, outbndry = localmesh->xstart;
  if ((global_flags & INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2)) {
    inbndry = outbndry = 1;
  }
  if (inner_boundary_flags & INVERT_BNDRY_ONE)
    inbndry = 1;
  if (outer_boundary_flags & INVERT_BNDRY_ONE)
    outbndry = 1;

  int nx = xe - xs + 1; // Number of X points on this processor

  // Get range of Y indices
  int ys = localmesh->ystart, ye = localmesh->yend;

  if (localmesh->hasBndryLowerY()) {
    if (include_yguards)
      ys = 0; // Mesh contains a lower boundary and we are solving in the guard cells

    ys += extra_yguards_lower;
  }
  if (localmesh->hasBndryUpperY()) {
    if (include_yguards)
      ye = localmesh->LocalNy -
           1; // Contains upper boundary and we are solving in the guard cells

    ye -= extra_yguards_upper;
  }

  const int ny = (ye - ys + 1); // Number of Y points
  const int nsys = nmode * ny;  // Number of systems of equations to solve
  const int nxny = nx * ny;     // Number of points in X-Y

  auto a3D = Matrix<dcomplex>(nsys, nx);
  auto b3D = Matrix<dcomplex>(nsys, nx);
  auto c3D = Matrix<dcomplex>(nsys, nx);

  auto xcmplx3D = Matrix<dcomplex>(nsys, nx);
  auto bcmplx3D = Matrix<dcomplex>(nsys, nx);

  if (dst) {
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d =
          Array<dcomplex>(localmesh->LocalNz); // ZFFT routine expects input of this length

      // Loop over X and Y indices, including boundaries but not guard cells.
      // (unless periodic in x)
      BOUT_OMP(for)
      for (int ind = 0; ind < nxny; ++ind) {
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        // Take DST in Z direction and put result in k1d

        if (((ix < inbndry) && (inner_boundary_flags & INVERT_SET) && localmesh->firstX()) ||
            ((localmesh->LocalNx - ix - 1 < outbndry) && (outer_boundary_flags & INVERT_SET) &&
             localmesh->lastX())) {
          // Use the values in x0 in the boundary
          DST(x0(ix, iy) + 1, localmesh->LocalNz - 2, std::begin(k1d));
        } else {
          DST(rhs(ix, iy) + 1, localmesh->LocalNz - 2, std::begin(k1d));
        }

        // Copy into array, transposing so kz is first index
        for (int kz = 0; kz < nmode; kz++) {
          bcmplx3D((iy - ys) * nmode + kz, ix - xs) = k1d[kz];
        }
      }

      // Get elements of the tridiagonal matrix
      // including boundary conditions
      BOUT_OMP(for nowait)
      for (int ind = 0; ind < nsys; ind++) {
        // ind = (iy - ys) * nmode + kz
        int iy = ys + ind / nmode;
        int kz = ind % nmode;

        BoutReal zlen = coords->dz * (localmesh->LocalNz - 3);
        BoutReal kwave =
            kz * 2.0 * PI / (2. * zlen); // wave number is 1/[rad]; DST has extra 2.

        tridagMatrix(&a3D(ind, 0), &b3D(ind, 0), &c3D(ind, 0), &bcmplx3D(ind, 0), iy,
                     kz,    // wave number index
                     kwave, // kwave (inverse wave length)
                     global_flags, inner_boundary_flags, outer_boundary_flags, &Acoef,
                     &C1coef, &C2coef, &Dcoef,
                     false); // Don't include guard cells in arrays
      }
    }

    // Solve tridiagonal systems
    cr->setCoefs(a3D, b3D, c3D);
    cr->solve(bcmplx3D, xcmplx3D);

    // FFT back to real space
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d =
          Array<dcomplex>(localmesh->LocalNz); // ZFFT routine expects input of this length

      BOUT_OMP(for nowait)
      for (int ind = 0; ind < nxny; ++ind) { // Loop over X and Y
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        for (int kz = 0; kz < nmode; kz++) {
          k1d[kz] = xcmplx3D((iy - ys) * nmode + kz, ix - xs);
        }

        for (int kz = nmode; kz < localmesh->LocalNz; kz++)
          k1d[kz] = 0.0; // Filtering out all higher harmonics

        DST_rev(std::begin(k1d), localmesh->LocalNz - 2, &x(ix, iy, 1));

        x(ix, iy, 0) = -x(ix, iy, 2);
        x(ix, iy, localmesh->LocalNz - 1) = -x(ix, iy, localmesh->LocalNz - 3);
      }
    }
  } else {
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>(localmesh->LocalNz / 2 +
                                 1); // ZFFT routine expects input of this length

      // Loop over X and Y indices, including boundaries but not guard cells
      // (unless periodic in x)

      BOUT_OMP(for)
      for (int ind = 0; ind < nxny; ++ind) {
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        // Take FFT in Z direction, apply shift, and put result in k1d

        if (((ix < inbndry) && (inner_boundary_flags & INVERT_SET) && localmesh->firstX()) ||
            ((localmesh->LocalNx - ix - 1 < outbndry) && (outer_boundary_flags & INVERT_SET) &&
             localmesh->lastX())) {
          // Use the values in x0 in the boundary
          rfft(x0(ix, iy), localmesh->LocalNz, std::begin(k1d));
        } else {
          rfft(rhs(ix, iy), localmesh->LocalNz, std::begin(k1d));
        }

        // Copy into array, transposing so kz is first index
        for (int kz = 0; kz < nmode; kz++)
          bcmplx3D((iy - ys) * nmode + kz, ix - xs) = k1d[kz];
      }

      // Get elements of the tridiagonal matrix
      // including boundary conditions
      BOUT_OMP(for nowait)
      for (int ind = 0; ind < nsys; ind++) {
        // ind = (iy - ys) * nmode + kz
        int iy = ys + ind / nmode;
        int kz = ind % nmode;

        BoutReal kwave = kz * 2.0 * PI / (coords->zlength()); // wave number is 1/[rad]
        tridagMatrix(&a3D(ind, 0), &b3D(ind, 0), &c3D(ind, 0), &bcmplx3D(ind, 0), iy,
                     kz,    // True for the component constant (DC) in Z
                     kwave, // Z wave number
                     global_flags, inner_boundary_flags, outer_boundary_flags, &Acoef,
                     &C1coef, &C2coef, &Dcoef,
                     false); // Don't include guard cells in arrays
      }
    }

    // Solve tridiagonal systems
    cr->setCoefs(a3D, b3D, c3D);
    cr->solve(bcmplx3D, xcmplx3D);

    // FFT back to real space
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>((localmesh->LocalNz) / 2 +
                                 1); // ZFFT routine expects input of this length

      bool zero_DC = false;
      if(global_flags & INVERT_ZERO_DC) {
        // No DC component
        zero_DC = true;
      }

      BOUT_OMP(for nowait)
      for (int ind = 0; ind < nxny; ++ind) { // Loop over X and Y
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        if (zero_DC) {
          k1d[0] = 0.;
        }

        for (int kz = zero_DC; kz < nmode; kz++)
          k1d[kz] = xcmplx3D((iy - ys) * nmode + kz, ix - xs);

        for (int kz = nmode; kz < localmesh->LocalNz / 2 + 1; kz++)
          k1d[kz] = 0.0; // Filtering out all higher harmonics

        irfft(std::begin(k1d), localmesh->LocalNz, x(ix, iy));
      }
    }
  }
  return x;
}
