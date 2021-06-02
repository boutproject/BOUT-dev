/**************************************************************************
 * Perpendicular Laplacian inversion. Parallel code using FFTs in z
 * and parallel cyclic reduction in x.
 *
 **************************************************************************
 * Copyright 2021 Joseph Parker
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
 * This solver incorporates code from the repo:
 * https://github.com/jihoonakang/parallel_tdma_cpp
 * by Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and
 * Technology Information
 *
 * If using this solver, please cite
 * @misc{kang2019ptdma,
 *  title  = {Parallel tri-diagonal matrix solver using cyclic reduction (CR), parallel CR
 *(PCR), and Thomas+PCR hybrid algorithm}, author = {Kang, Ji-Hoon}, url    =
 *https://github.com/jihoonakang/parallel_tdma_cpp}, year   = {2019}
 * }
 *
 **************************************************************************/

#include "pcr.hxx"
#include "globals.hxx"

#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout/openmpwrap.hxx>
#include <bout/sys/timer.hxx>
#include <boutexception.hxx>
#include <cmath>
#include <fft.hxx>
#include <lapack_routines.hxx>
#include <utils.hxx>

#include "boutcomm.hxx"
#include <output.hxx>

#include <bout/scorepwrapper.hxx>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

using namespace std;

bout::ArgumentHelper<LaplacePCR>::ArgumentHelper(Options& options)
    : bout::ArgumentHelper<Laplacian>(options),
      dst(options["dst"].doc("Use DST instead of FFT").withDefault(false)) {}

LaplacePCR::LaplacePCR(Options* opt, CELL_LOC loc, Mesh* mesh_in)
    : Laplacian(opt, loc, mesh_in), Acoef(0.0, localmesh), C1coef(1.0, localmesh),
      C2coef(1.0, localmesh), Dcoef(1.0, localmesh), nmode(maxmode + 1),
      ncx(localmesh->LocalNx), ny(localmesh->LocalNy), avec(ny, nmode, ncx),
      bvec(ny, nmode, ncx), cvec(ny, nmode, ncx) {

  Acoef.setLocation(location);
  C1coef.setLocation(location);
  C2coef.setLocation(location);
  Dcoef.setLocation(location);

  // Number of X procs must be a power of 2
  const int nxpe = localmesh->getNXPE();
  if (!is_pow2(nxpe)) {
    throw BoutException("LaplacePCR error: NXPE must be a power of 2");
  }

  // Number of x points must be a power of 2
  if (!is_pow2(localmesh->GlobalNxNoBoundaries)) {
    throw BoutException("LaplacePCR error: GlobalNx must be a power of 2");
  }

  Acoef.setLocation(location);
  C1coef.setLocation(location);
  C2coef.setLocation(location);
  Dcoef.setLocation(location);

  // Get options

  OPTION(opt, dst, false);

  if (dst) {
    nmode = localmesh->LocalNz - 2;
  } else {
    nmode = maxmode + 1; // Number of Z modes. maxmode set in
                         // invert_laplace.cxx from options
  }

  // Note nmode == nsys of cyclic_reduction

  // Allocate arrays

  xs = localmesh->xstart; // Starting X index
  if (localmesh->firstX()
      && !localmesh->periodicX) { // Only want to include guard cells at boundaries
                                  // (unless periodic in x)
    xs = 0;
  }
  xe = localmesh->xend; // Last X index
  if (localmesh->lastX()
      && !localmesh->periodicX) { // Only want to include guard cells at boundaries
                                  // (unless periodic in x)
    xe = localmesh->LocalNx - 1;
  }
  int n = xe - xs + 1; // Number of X points on this processor,
                       // including boundaries but not guard cells

  a.reallocate(nmode, n);
  b.reallocate(nmode, n);
  c.reallocate(nmode, n);
  xcmplx.reallocate(nmode, n);
  bcmplx.reallocate(nmode, n);

  xproc = localmesh->getXProcIndex(); // Local rank in x proc space
  const int yproc = localmesh->getYProcIndex();
  nprocs = localmesh->getNXPE();                    // Number of processors in x
  myrank = yproc * nprocs + xproc;                  // Global rank for communication
  n_mpi = localmesh->GlobalNxNoBoundaries / nprocs; // Number of internal x
                                                    // grid points per proc
}

FieldPerp LaplacePCR::solve(const FieldPerp& rhs, const FieldPerp& x0) {
  ASSERT1(localmesh == rhs.getMesh() && localmesh == x0.getMesh());
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  FieldPerp x{emptyFrom(rhs)}; // Result

  int jy = rhs.getIndex(); // Get the Y index
  x.setIndex(jy);

  // Get the width of the boundary

  // If the flags to assign that only one guard cell should be used is set
  int inbndry = localmesh->xstart;
  int outbndry = localmesh->xstart;
  if (((global_flags & INVERT_BOTH_BNDRY_ONE) != 0) || (localmesh->xstart < 2)) {
    inbndry = outbndry = 1;
  }
  if ((inner_boundary_flags & INVERT_BNDRY_ONE) != 0) {
    inbndry = 1;
  }
  if ((outer_boundary_flags & INVERT_BNDRY_ONE) != 0) {
    outbndry = 1;
  }

  if (dst) {
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>(
          localmesh->LocalNz); // ZFFT routine expects input of this length

      // Loop over X indices, including boundaries but not guard cells. (unless periodic
      // in x)
      BOUT_OMP(for)
      for (int ix = xs; ix <= xe; ix++) {
        // Take DST in Z direction and put result in k1d

        if (((ix < inbndry) && ((inner_boundary_flags & INVERT_SET) != 0)
             && localmesh->firstX())
            || ((localmesh->LocalNx - ix - 1 < outbndry)
                && ((outer_boundary_flags & INVERT_SET) != 0) && localmesh->lastX())) {
          // Use the values in x0 in the boundary
          DST(x0[ix] + 1, localmesh->LocalNz - 2, std::begin(k1d));
        } else {
          DST(rhs[ix] + 1, localmesh->LocalNz - 2, std::begin(k1d));
        }

        // Copy into array, transposing so kz is first index
        for (int kz = 0; kz < nmode; kz++) {
          bcmplx(kz, ix - xs) = k1d[kz];
        }
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
    cr_pcr_solver(a, b, c, bcmplx, xcmplx);

    // FFT back to real space
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>(
          localmesh->LocalNz); // ZFFT routine expects input of this length

      BOUT_OMP(for nowait)
      for (int ix = xs; ix <= xe; ix++) {
        for (int kz = 0; kz < nmode; kz++) {
          k1d[kz] = xcmplx(kz, ix - xs);
        }

        for (int kz = nmode; kz < (localmesh->LocalNz); kz++) {
          k1d[kz] = 0.0; // Filtering out all higher harmonics
        }

        DST_rev(std::begin(k1d), localmesh->LocalNz - 2, x[ix] + 1);

        x(ix, 0) = -x(ix, 2);
        x(ix, localmesh->LocalNz - 1) = -x(ix, localmesh->LocalNz - 3);
      }
    }
  } else {
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>((localmesh->LocalNz) / 2
                                 + 1); // ZFFT routine expects input of this length

      // Loop over X indices, including boundaries but not guard cells (unless periodic in
      // x)
      BOUT_OMP(for)
      for (int ix = xs; ix <= xe; ix++) {
        // Take FFT in Z direction, apply shift, and put result in k1d

        if (((ix < inbndry) && ((inner_boundary_flags & INVERT_SET) != 0)
             && localmesh->firstX())
            || ((localmesh->LocalNx - ix - 1 < outbndry)
                && ((outer_boundary_flags & INVERT_SET) != 0) && localmesh->lastX())) {
          // Use the values in x0 in the boundary
          rfft(x0[ix], localmesh->LocalNz, std::begin(k1d));
        } else {
          rfft(rhs[ix], localmesh->LocalNz, std::begin(k1d));
        }

        // Copy into array, transposing so kz is first index
        for (int kz = 0; kz < nmode; kz++) {
          bcmplx(kz, ix - xs) = k1d[kz];
        }
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
    cr_pcr_solver(a, b, c, bcmplx, xcmplx);

    // FFT back to real space
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>((localmesh->LocalNz) / 2
                                 + 1); // ZFFT routine expects input of this length

      const bool zero_DC = (global_flags & INVERT_ZERO_DC) != 0;

      BOUT_OMP(for nowait)
      for (int ix = xs; ix <= xe; ix++) {
        if (zero_DC) {
          k1d[0] = 0.;
        }

        for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
          k1d[kz] = xcmplx(kz, ix - xs);
        }

        for (int kz = nmode; kz < (localmesh->LocalNz) / 2 + 1; kz++) {
          k1d[kz] = 0.0; // Filtering out all higher harmonics
        }

        irfft(std::begin(k1d), localmesh->LocalNz, x[ix]);
      }
    }
  }

  checkData(x);

  return x;
}

Field3D LaplacePCR::solve(const Field3D& rhs, const Field3D& x0) {
  TRACE("LaplacePCR::solve(Field3D, Field3D)");

  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  ASSERT1(localmesh == rhs.getMesh() && localmesh == x0.getMesh());

  Timer timer("invert");

  Field3D x{emptyFrom(rhs)}; // Result

  // Get the width of the boundary

  // If the flags to assign that only one guard cell should be used is set
  int inbndry = localmesh->xstart;
  int outbndry = localmesh->xstart;
  if (((global_flags & INVERT_BOTH_BNDRY_ONE) != 0) || (localmesh->xstart < 2)) {
    inbndry = outbndry = 1;
  }
  if ((inner_boundary_flags & INVERT_BNDRY_ONE) != 0) {
    inbndry = 1;
  }
  if ((outer_boundary_flags & INVERT_BNDRY_ONE) != 0) {
    outbndry = 1;
  }

  int nx = xe - xs + 1; // Number of X points on this processor

  // Get range of Y indices
  int ys = localmesh->ystart;
  int ye = localmesh->yend;

  if (localmesh->hasBndryLowerY()) {
    if (include_yguards) {
      ys = 0; // Mesh contains a lower boundary and we are solving in the guard cells
    }

    ys += extra_yguards_lower;
  }
  if (localmesh->hasBndryUpperY()) {
    if (include_yguards) {
      ye = localmesh->LocalNy - 1; // Contains upper boundary and we are
                                   // solving in the guard cells:
    }
    ye -= extra_yguards_upper;
  }

  const int ny = (ye - ys + 1); // Number of Y points
  nsys = nmode * ny;            // Number of systems of equations to solve
  const int nxny = nx * ny;     // Number of points in X-Y

  auto a3D = Matrix<dcomplex>(nsys, nx);
  auto b3D = Matrix<dcomplex>(nsys, nx);
  auto c3D = Matrix<dcomplex>(nsys, nx);

  auto xcmplx3D = Matrix<dcomplex>(nsys, nx);
  auto bcmplx3D = Matrix<dcomplex>(nsys, nx);

  if (dst) {
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>(
          localmesh->LocalNz); // ZFFT routine expects input of this length

      // Loop over X and Y indices, including boundaries but not guard cells.
      // (unless periodic in x)
      BOUT_OMP(for)
      for (int ind = 0; ind < nxny; ++ind) {
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        // Take DST in Z direction and put result in k1d

        if (((ix < inbndry) && ((inner_boundary_flags & INVERT_SET) != 0)
             && localmesh->firstX())
            || ((localmesh->LocalNx - ix - 1 < outbndry)
                && ((outer_boundary_flags & INVERT_SET) != 0) && localmesh->lastX())) {
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
    cr_pcr_solver(a3D, b3D, c3D, bcmplx3D, xcmplx3D);

    // FFT back to real space
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>(
          localmesh->LocalNz); // ZFFT routine expects input of this length

      BOUT_OMP(for nowait)
      for (int ind = 0; ind < nxny; ++ind) { // Loop over X and Y
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        for (int kz = 0; kz < nmode; kz++) {
          k1d[kz] = xcmplx3D((iy - ys) * nmode + kz, ix - xs);
        }

        for (int kz = nmode; kz < localmesh->LocalNz; kz++) {
          k1d[kz] = 0.0; // Filtering out all higher harmonics
        }

        DST_rev(std::begin(k1d), localmesh->LocalNz - 2, &x(ix, iy, 1));

        x(ix, iy, 0) = -x(ix, iy, 2);
        x(ix, iy, localmesh->LocalNz - 1) = -x(ix, iy, localmesh->LocalNz - 3);
      }
    }
  } else {
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>(localmesh->LocalNz / 2
                                 + 1); // ZFFT routine expects input of this length

      // Loop over X and Y indices, including boundaries but not guard cells
      // (unless periodic in x)

      BOUT_OMP(for)
      for (int ind = 0; ind < nxny; ++ind) {
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        // Take FFT in Z direction, apply shift, and put result in k1d

        if (((ix < inbndry) && ((inner_boundary_flags & INVERT_SET) != 0)
             && localmesh->firstX())
            || ((localmesh->LocalNx - ix - 1 < outbndry)
                && ((outer_boundary_flags & INVERT_SET) != 0) && localmesh->lastX())) {
          // Use the values in x0 in the boundary
          rfft(x0(ix, iy), localmesh->LocalNz, std::begin(k1d));
        } else {
          rfft(rhs(ix, iy), localmesh->LocalNz, std::begin(k1d));
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
    cr_pcr_solver(a3D, b3D, c3D, bcmplx3D, xcmplx3D);

    // FFT back to real space
    BOUT_OMP(parallel) {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>((localmesh->LocalNz) / 2
                                 + 1); // ZFFT routine expects input of this length

      const bool zero_DC = (global_flags & INVERT_ZERO_DC) != 0;

      BOUT_OMP(for nowait)
      for (int ind = 0; ind < nxny; ++ind) { // Loop over X and Y
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        if (zero_DC) {
          k1d[0] = 0.;
        }

        for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
          k1d[kz] = xcmplx3D((iy - ys) * nmode + kz, ix - xs);
        }

        for (int kz = nmode; kz < localmesh->LocalNz / 2 + 1; kz++) {
          k1d[kz] = 0.0; // Filtering out all higher harmonics
        }

        irfft(std::begin(k1d), localmesh->LocalNz, x(ix, iy));
      }
    }
  }

  checkData(x);

  return x;
}

/**
 * @brief   CR-PCR solver: cr_forward_multiple + pcr_forward_single + cr_backward_multiple
 * @param   a_mpi (input) Lower off-diagonal coeff., which is assigned to local private
 * pointer a
 * @param   b_mpi (input) Diagonal coeff., which is assigned to local private pointer b
 * @param   c_mpi (input) Upper off-diagonal coeff.,, which is assigned to local private
 * pointer c
 * @param   r_mpi (input) RHS vector, which is assigned to local private pointer r
 * @param   x_mpi (output) Solution vector, which is assigned to local private pointer x
 */
void LaplacePCR ::cr_pcr_solver(Matrix<dcomplex>& a_mpi, Matrix<dcomplex>& b_mpi,
                                Matrix<dcomplex>& c_mpi, Matrix<dcomplex>& r_mpi,
                                Matrix<dcomplex>& x_mpi) {

  const int xstart = localmesh->xstart;
  const int xend = localmesh->xend;
  const int nx = xend - xstart + 1; // number of interior points

  // Handle boundary points so that the PCR algorithm works with arrays of
  // the same size on each rank.
  // Note that this modifies the coefficients of b and r in the first and last
  // interior rows. We can continue to use b_mpi and r_mpi arrays directly as
  // their original values are no longer required.
  eliminate_boundary_rows(a_mpi, b_mpi, c_mpi, r_mpi);

  // nsys = nmode * ny;  // Number of systems of equations to solve
  aa.reallocate(nsys, nx + 2);
  bb.reallocate(nsys, nx + 2);
  cc.reallocate(nsys, nx + 2);
  r.reallocate(nsys, nx + 2);
  x.reallocate(nsys, nx + 2);

  for (int kz = 0; kz < nsys; kz++) {
    aa(kz, 0) = 0;
    bb(kz, 0) = 1;
    cc(kz, 0) = 0;
    r(kz, 0) = 0;
    x(kz, 0) = 0;
    for (int ix = 0; ix < nx; ix++) {
      // The offset xstart - xs ensures that this copies interior points.
      // If a proc has boundary points, these are included in *_mpi, but we
      // don't want to copy them.
      // xs = xstart if a proc has no boundary points
      // xs = 0 if a proc has boundary points
      aa(kz, ix + 1) = a_mpi(kz, ix + xstart - xs);
      bb(kz, ix + 1) = b_mpi(kz, ix + xstart - xs);
      cc(kz, ix + 1) = c_mpi(kz, ix + xstart - xs);
      r(kz, ix + 1) = r_mpi(kz, ix + xstart - xs);
      x(kz, ix + 1) = x_mpi(kz, ix + xstart - xs);
    }
    aa(kz, nx + 1) = 0;
    bb(kz, nx + 1) = 1;
    cc(kz, nx + 1) = 0;
    r(kz, nx + 1) = 0;
    x(kz, nx + 1) = 0;
  }

  // Perform parallel cyclic reduction
  cr_forward_multiple_row(aa, bb, cc, r);
  pcr_forward_single_row(aa, bb, cc, r, x); // Including 2x2 solver
  cr_backward_multiple_row(aa, bb, cc, r, x);
  // End of PCR

  // Copy solution back to bout format - this is correct on interior rows, but
  // not boundary rows
  for (int kz = 0; kz < nsys; kz++) {
    for (int ix = 0; ix < nx; ix++) {
      x_mpi(kz, ix + xstart - xs) = x(kz, ix + 1);
    }
  }

  // Ensure solution is also correct on boundary rows
  apply_boundary_conditions(a_mpi, b_mpi, c_mpi, r_mpi, x_mpi);
}

/**
 * Eliminate boundary rows - perform row elimination to uncouple the first and
 * last interior rows from their respective boundary rows. This is necessary
 * to ensure we pass a square system of interior rows to the PCR library.
 */
void LaplacePCR ::eliminate_boundary_rows(const Matrix<dcomplex>& a, Matrix<dcomplex>& b,
                                          const Matrix<dcomplex>& c,
                                          Matrix<dcomplex>& r) {

  if (localmesh->firstX()) {
    // x index is first interior row
    const int xstart = localmesh->xstart;
    for (int kz = 0; kz < nsys; kz++) {
      b(kz, xstart) =
          b(kz, xstart) - c(kz, xstart - 1) * a(kz, xstart) / b(kz, xstart - 1);
      r(kz, xstart) =
          r(kz, xstart) - r(kz, xstart - 1) * a(kz, xstart) / b(kz, xstart - 1);
      // Row elimination would set a to zero, but value is unused:
      // a(kz,xstart) = 0.0;
    }
  }
  if (localmesh->lastX()) {
    int n = xe - xs + 1; // actual length of array
    int xind = n - localmesh->xstart - 1;
    for (int kz = 0; kz < nsys; kz++) {
      // x index is last interior row
      b(kz, xind) = b(kz, xind) - c(kz, xind) * a(kz, xind + 1) / b(kz, xind + 1);
      r(kz, xind) = r(kz, xind) - c(kz, xind) * r(kz, xind + 1) / b(kz, xind + 1);
      // Row elimination would set c to zero, but value is unused:
      // c(kz,xind) = 0.0;
    }
  }
}

/**
 * Apply the boundary conditions on the first and last X processors
 */
void LaplacePCR ::apply_boundary_conditions(const Matrix<dcomplex>& a,
                                            const Matrix<dcomplex>& b,
                                            const Matrix<dcomplex>& c,
                                            const Matrix<dcomplex>& r,
                                            Matrix<dcomplex>& x) {

  if (localmesh->firstX()) {
    for (int kz = 0; kz < nsys; kz++) {
      for (int ix = localmesh->xstart - 1; ix >= 0; ix--) {
        x(kz, ix) = (r(kz, ix) - c(kz, ix) * x(kz, ix + 1)) / b(kz, ix);
      }
    }
  }
  if (localmesh->lastX()) {
    int n = xe - xs + 1; // actual length of array
    for (int kz = 0; kz < nsys; kz++) {
      for (int ix = n - localmesh->xstart; ix < n; ix++) {
        x(kz, ix) = (r(kz, ix) - a(kz, ix) * x(kz, ix - 1)) / b(kz, ix);
      }
    }
  }
}

/**
 * @brief   Forward elimination of CR until a single row per MPI process remains.
 * @details After a single row per MPI process remains, PCR or CR between a single row is
 * performed.
 */
void LaplacePCR ::cr_forward_multiple_row(Matrix<dcomplex>& a, Matrix<dcomplex>& b,
                                          Matrix<dcomplex>& c,
                                          Matrix<dcomplex>& r) const {
  MPI_Comm comm = BoutComm::get();
  Array<dcomplex> alpha(nsys);
  Array<dcomplex> gamma(nsys);
  Array<dcomplex> sbuf(4 * nsys);
  Array<dcomplex> rbuf(4 * nsys);

  MPI_Status status;
  MPI_Status status1;
  Array<MPI_Request> request(2);

  /// Variable nlevel is used to indicates when single row remains.
  const int nlevel = log2(n_mpi);
  int dist_row = 1;
  int dist2_row = 2;

  for (int l = 0; l < nlevel; l++) {
    const int start = dist2_row;
    /// Data exchange is performed using MPI send/recv for each succesive reduction
    if (xproc < nprocs - 1) {
      MPI_Irecv(&rbuf[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank + 1, 0, comm, &request[0]);
    }
    if (xproc > 0) {
      for (int kz = 0; kz < nsys; kz++) {
        sbuf[0 + 4 * kz] = a(kz, dist_row);
        sbuf[1 + 4 * kz] = b(kz, dist_row);
        sbuf[2 + 4 * kz] = c(kz, dist_row);
        sbuf[3 + 4 * kz] = r(kz, dist_row);
      }
      MPI_Isend(&sbuf[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank - 1, 0, comm, &request[1]);
    }
    if (xproc < nprocs - 1) {
      MPI_Wait(&request[0], &status1);
      for (int kz = 0; kz < nsys; kz++) {
        a(kz, n_mpi + 1) = rbuf[0 + 4 * kz];
        b(kz, n_mpi + 1) = rbuf[1 + 4 * kz];
        c(kz, n_mpi + 1) = rbuf[2 + 4 * kz];
        r(kz, n_mpi + 1) = rbuf[3 + 4 * kz];
      }
    }

    /// Odd rows of remained rows are reduced to even rows of remained rows in each
    /// reduction step. Index in of global last row is out of range, but we treat it as a
    /// = c = r = 0 and b = 1 in main function.
    for (int i = start; i <= n_mpi; i += dist2_row) {
      const int ip = i - dist_row;
      const int in = min(i + dist_row, n_mpi + 1);
      for (int kz = 0; kz < nsys; kz++) {
        alpha[kz] = -a(kz, i) / b(kz, ip);
        gamma[kz] = -c(kz, i) / b(kz, in);

        b(kz, i) += (alpha[kz] * c(kz, ip) + gamma[kz] * a(kz, in));
        a(kz, i) = alpha[kz] * a(kz, ip);
        c(kz, i) = gamma[kz] * c(kz, in);
        r(kz, i) += (alpha[kz] * r(kz, ip) + gamma[kz] * r(kz, in));
      }
    }
    /// As reduction continues, the indices of required coefficients doubles.
    dist2_row *= 2;
    dist_row *= 2;

    if (xproc > 0) {
      MPI_Wait(&request[1], &status);
    }
  }
}

/**
 * @brief   Backward substitution of CR after single-row solution per MPI process is
 * obtained.
 */
void LaplacePCR ::cr_backward_multiple_row(Matrix<dcomplex>& a, Matrix<dcomplex>& b,
                                           Matrix<dcomplex>& c, Matrix<dcomplex>& r,
                                           Matrix<dcomplex>& x) const {
  MPI_Comm comm = BoutComm::get();

  MPI_Status status;
  MPI_Request request[2];
  auto recvvec = Array<dcomplex>(nsys);
  auto sendvec = Array<dcomplex>(nsys);

  const int nlevel = log2(n_mpi);
  int dist_row = n_mpi / 2;

  /// Each rank requires a solution on last row of previous rank.
  if (xproc > 0) {
    MPI_Irecv(&recvvec[0], nsys, MPI_DOUBLE_COMPLEX, myrank - 1, 100, comm, request);
  }
  if (xproc < nprocs - 1) {
    for (int kz = 0; kz < nsys; kz++) {
      sendvec[kz] = x(kz, n_mpi);
    }
    MPI_Isend(&sendvec[0], nsys, MPI_DOUBLE_COMPLEX, myrank + 1, 100, comm, request + 1);
  }
  if (xproc > 0) {
    MPI_Wait(request, &status);
    for (int kz = 0; kz < nsys; kz++) {
      x(kz, 0) = recvvec[kz];
    }
  }
  for (int l = nlevel - 1; l >= 0; l--) {
    const int dist2_row = dist_row * 2;
    for (int i = n_mpi - dist_row; i >= 0; i -= dist2_row) {
      const int ip = i - dist_row;
      const int in = i + dist_row;
      for (int kz = 0; kz < nsys; kz++) {
        x(kz, i) = r(kz, i) - c(kz, i) * x(kz, in) - a(kz, i) * x(kz, ip);
        x(kz, i) = x(kz, i) / b(kz, i);
      }
    }
    dist_row = dist_row / 2;
  }
  if (xproc < nprocs - 1) {
    MPI_Wait(request + 1, &status);
  }
}

/**
 * @brief   PCR between a single row per MPI process and 2x2 matrix solver between i and
 * i+nprocs/2 rows.
 */
void LaplacePCR ::pcr_forward_single_row(Matrix<dcomplex>& a, Matrix<dcomplex>& b,
                                         Matrix<dcomplex>& c, Matrix<dcomplex>& r,
                                         Matrix<dcomplex>& x) const {

  Array<dcomplex> alpha(nsys);
  Array<dcomplex> gamma(nsys);
  Array<dcomplex> sbuf(4 * nsys);
  Array<dcomplex> rbuf0(4 * nsys);
  Array<dcomplex> rbuf1(4 * nsys);

  MPI_Status status;
  Array<MPI_Request> request(4);
  MPI_Comm comm = BoutComm::get();

  const int nlevel = log2(nprocs);
  const int nhprocs = nprocs / 2;
  int dist_rank = 1;
  int dist2_rank = 2;

  /// Parallel cyclic reduction continues until 2x2 matrix are made between a pair of
  /// rank, (myrank, myrank+nhprocs).
  for (int l = 0; l < nlevel - 1; l++) {

    /// Rank is newly calculated in each level to find communication pair.
    /// Nprocs is also newly calculated as myrank is changed.
    const int myrank_level = xproc / dist_rank;
    const int nprocs_level = nprocs / dist_rank;

    /// All rows exchange data for reduction and perform reduction successively.
    /// Coefficients are updated for every rows.
    for (int kz = 0; kz < nsys; kz++) {
      sbuf[0 + 4 * kz] = a(kz, n_mpi);
      sbuf[1 + 4 * kz] = b(kz, n_mpi);
      sbuf[2 + 4 * kz] = c(kz, n_mpi);
      sbuf[3 + 4 * kz] = r(kz, n_mpi);
    }

    if ((myrank_level + 1) % 2 == 0) {
      if (xproc + dist_rank < nprocs) {
        MPI_Irecv(&rbuf1[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank + dist_rank, 202, comm,
                  &request[0]);
        MPI_Isend(&sbuf[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank + dist_rank, 203, comm,
                  &request[1]);
      }
      if (xproc - dist_rank >= 0) {
        MPI_Irecv(&rbuf0[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank - dist_rank, 200, comm,
                  &request[2]);
        MPI_Isend(&sbuf[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank - dist_rank, 201, comm,
                  &request[3]);
      }
      if (xproc + dist_rank < nprocs) {
        MPI_Wait(&request[0], &status);
        for (int kz = 0; kz < nsys; kz++) {
          a(kz, n_mpi + 1) = rbuf1[0 + 4 * kz];
          b(kz, n_mpi + 1) = rbuf1[1 + 4 * kz];
          c(kz, n_mpi + 1) = rbuf1[2 + 4 * kz];
          r(kz, n_mpi + 1) = rbuf1[3 + 4 * kz];
        }
        MPI_Wait(&request[1], &status);
      }
      if (xproc - dist_rank >= 0) {
        MPI_Wait(&request[2], &status);
        for (int kz = 0; kz < nsys; kz++) {
          a(kz, 0) = rbuf0[0 + 4 * kz];
          b(kz, 0) = rbuf0[1 + 4 * kz];
          c(kz, 0) = rbuf0[2 + 4 * kz];
          r(kz, 0) = rbuf0[3 + 4 * kz];
        }
        MPI_Wait(&request[3], &status);
      }
    } else if ((myrank_level + 1) % 2 == 1) {
      if (xproc + dist_rank < nprocs) {
        MPI_Irecv(&rbuf1[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank + dist_rank, 201, comm,
                  &request[0]);
        MPI_Isend(&sbuf[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank + dist_rank, 200, comm,
                  &request[1]);
      }
      if (xproc - dist_rank >= 0) {
        MPI_Irecv(&rbuf0[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank - dist_rank, 203, comm,
                  &request[2]);
        MPI_Isend(&sbuf[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank - dist_rank, 202, comm,
                  &request[3]);
      }
      if (xproc + dist_rank < nprocs) {
        MPI_Wait(&request[0], &status);
        for (int kz = 0; kz < nsys; kz++) {
          a(kz, n_mpi + 1) = rbuf1[0 + 4 * kz];
          b(kz, n_mpi + 1) = rbuf1[1 + 4 * kz];
          c(kz, n_mpi + 1) = rbuf1[2 + 4 * kz];
          r(kz, n_mpi + 1) = rbuf1[3 + 4 * kz];
        }
        MPI_Wait(&request[1], &status);
      }
      if (xproc - dist_rank >= 0) {
        MPI_Wait(&request[2], &status);
        for (int kz = 0; kz < nsys; kz++) {
          a(kz, 0) = rbuf0[0 + 4 * kz];
          b(kz, 0) = rbuf0[1 + 4 * kz];
          c(kz, 0) = rbuf0[2 + 4 * kz];
          r(kz, 0) = rbuf0[3 + 4 * kz];
        }
        MPI_Wait(&request[3], &status);
      }
    }

    const int i = n_mpi;
    const int ip = 0;
    const int in = i + 1;
    if (myrank_level == 0) {
      for (int kz = 0; kz < nsys; kz++) {
        alpha[kz] = 0.0;
      }
    } else {
      for (int kz = 0; kz < nsys; kz++) {
        alpha[kz] = -a(kz, i) / b(kz, ip);
      }
    }
    if (myrank_level == nprocs_level - 1) {
      for (int kz = 0; kz < nsys; kz++) {
        gamma[kz] = 0.0;
      }
    } else {
      for (int kz = 0; kz < nsys; kz++) {
        gamma[kz] = -c(kz, i) / b(kz, in);
      }
    }

    for (int kz = 0; kz < nsys; kz++) {
      b(kz, i) += (alpha[kz] * c(kz, ip) + gamma[kz] * a(kz, in));
      a(kz, i) = alpha[kz] * a(kz, ip);
      c(kz, i) = gamma[kz] * c(kz, in);
      r(kz, i) += (alpha[kz] * r(kz, ip) + gamma[kz] * r(kz, in));
    }

    dist_rank *= 2;
    dist2_rank *= 2;
  }

  /// Solving 2x2 matrix. All pair of ranks, myrank and myrank+nhprocs, solves it
  /// simultaneously.
  for (int kz = 0; kz < nsys; kz++) {
    sbuf[0 + 4 * kz] = a(kz, n_mpi);
    sbuf[1 + 4 * kz] = b(kz, n_mpi);
    sbuf[2 + 4 * kz] = c(kz, n_mpi);
    sbuf[3 + 4 * kz] = r(kz, n_mpi);
  }
  if (xproc < nhprocs) {
    MPI_Irecv(&rbuf1[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank + nhprocs, 300, comm,
              &request[0]);
    MPI_Isend(&sbuf[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank + nhprocs, 301, comm,
              &request[1]);

    MPI_Wait(&request[0], &status);
    for (int kz = 0; kz < nsys; kz++) {
      a(kz, n_mpi + 1) = rbuf1[0 + 4 * kz];
      b(kz, n_mpi + 1) = rbuf1[1 + 4 * kz];
      c(kz, n_mpi + 1) = rbuf1[2 + 4 * kz];
      r(kz, n_mpi + 1) = rbuf1[3 + 4 * kz];
    }

    const int i = n_mpi;
    const int in = n_mpi + 1;

    for (int kz = 0; kz < nsys; kz++) {
      const dcomplex det = b(kz, i) * b(kz, in) - c(kz, i) * a(kz, in);
      x(kz, i) = (r(kz, i) * b(kz, in) - r(kz, in) * c(kz, i)) / det;
      x(kz, in) = (r(kz, in) * b(kz, i) - r(kz, i) * a(kz, in)) / det;
    }
    MPI_Wait(&request[1], &status);

  } else if (xproc >= nhprocs) {

    // nhprocs=0 if and only if NXPE=1. This check skips communication and
    // allows the serial case to work:
    if (nhprocs > 0) {
      MPI_Irecv(&rbuf0[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank - nhprocs, 301, comm,
                &request[2]);
      MPI_Isend(&sbuf[0], 4 * nsys, MPI_DOUBLE_COMPLEX, myrank - nhprocs, 300, comm,
                &request[3]);

      MPI_Wait(&request[2], &status);
      for (int kz = 0; kz < nsys; kz++) {
        a(kz, 0) = rbuf0[0 + 4 * kz];
        b(kz, 0) = rbuf0[1 + 4 * kz];
        c(kz, 0) = rbuf0[2 + 4 * kz];
        r(kz, 0) = rbuf0[3 + 4 * kz];
      }
    }

    const int ip = 0;
    const int i = n_mpi;

    for (int kz = 0; kz < nsys; kz++) {
      const dcomplex det = b(kz, ip) * b(kz, i) - c(kz, ip) * a(kz, i);
      x(kz, ip) = (r(kz, ip) * b(kz, i) - r(kz, i) * c(kz, ip)) / det;
      x(kz, i) = (r(kz, i) * b(kz, ip) - r(kz, ip) * a(kz, i)) / det;
    }

    if (nhprocs > 0) {
      MPI_Wait(&request[3], &status);
    }
  }
}

/**
 * @brief   Solution check
 * @param   *a_ver Coefficients of a with original values
 * @param   *b_ver Coefficients of b with original values
 * @param   *c_ver Coefficients of c with original values
 * @param   *r_ver RHS vector with original values
 * @param   *x_sol Solution vector
 */
void LaplacePCR ::verify_solution(const Matrix<dcomplex>& a_ver,
                                  const Matrix<dcomplex>& b_ver,
                                  const Matrix<dcomplex>& c_ver,
                                  const Matrix<dcomplex>& r_ver,
                                  const Matrix<dcomplex>& x_sol) const {
  output.write("Verify solution\n");
  const int nx = xe - xs + 1; // Number of X points on this processor,
                              // including boundaries but not guard cells
  Matrix<dcomplex> y_ver(nsys, nx + 2);
  Matrix<dcomplex> error(nsys, nx + 2);

  MPI_Status status;
  Array<MPI_Request> request(4);
  Array<dcomplex> sbufup(nsys);
  Array<dcomplex> sbufdown(nsys);
  Array<dcomplex> rbufup(nsys);
  Array<dcomplex> rbufdown(nsys);

  // nsys = nmode * ny;  // Number of systems of equations to solve
  Matrix<dcomplex> x_ver(nsys, nx + 2);

  for (int kz = 0; kz < nsys; kz++) {
    for (int ix = 0; ix < nx; ix++) {
      x_ver(kz, ix + 1) = x_sol(kz, ix);
    }
  }

  if (xproc > 0) {
    MPI_Irecv(&rbufdown[0], nsys, MPI_DOUBLE_COMPLEX, myrank - 1, 901, MPI_COMM_WORLD,
              &request[1]);
    for (int kz = 0; kz < nsys; kz++) {
      sbufdown[kz] = x_ver(kz, 1);
    }
    MPI_Isend(&sbufdown[0], nsys, MPI_DOUBLE_COMPLEX, myrank - 1, 900, MPI_COMM_WORLD,
              &request[0]);
  }
  if (xproc < nprocs - 1) {
    MPI_Irecv(&rbufup[0], nsys, MPI_DOUBLE_COMPLEX, myrank + 1, 900, MPI_COMM_WORLD,
              &request[3]);
    for (int kz = 0; kz < nsys; kz++) {
      sbufup[kz] = x_ver(kz, nx);
    }
    MPI_Isend(&sbufup[0], nsys, MPI_DOUBLE_COMPLEX, myrank + 1, 901, MPI_COMM_WORLD,
              &request[2]);
  }

  if (xproc > 0) {
    MPI_Wait(&request[0], &status);
    MPI_Wait(&request[1], &status);
    for (int kz = 0; kz < nsys; kz++) {
      x_ver(kz, 0) = rbufdown[kz];
    }
  }
  if (xproc < nprocs - 1) {
    MPI_Wait(&request[2], &status);
    MPI_Wait(&request[3], &status);
    for (int kz = 0; kz < nsys; kz++) {
      x_ver(kz, nx + 1) = rbufup[kz];
    }
  }

  BoutReal max_error = 0.0;
  for (int kz = 0; kz < nsys; kz++) {
    for (int i = 0; i < nx; i++) {
      y_ver(kz, i) = a_ver(kz, i) * x_ver(kz, i) + b_ver(kz, i) * x_ver(kz, i + 1)
                     + c_ver(kz, i) * x_ver(kz, i + 2);
      error(kz, i) = y_ver(kz, i) - r_ver(kz, i);
      max_error = std::max(max_error, std::abs(error(kz, i)));
      output.write("abs error {}, r={}, y={}, kz {}, i {},  a={}, b={}, c={}, x-= {}, "
                   "x={}, x+ = {}\n",
                   error(kz, i).real(), r_ver(kz, i).real(), y_ver(kz, i).real(), kz, i,
                   a_ver(kz, i).real(), b_ver(kz, i).real(), c_ver(kz, i).real(),
                   x_ver(kz, i).real(), x_ver(kz, i + 1).real(), x_ver(kz, i + 2).real());
    }
  }
  output.write("max abs error {}\n", max_error);
}
