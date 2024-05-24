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

#include "bout/build_defines.hxx"

#if not BOUT_USE_METRIC_3D

#include "cyclic_laplace.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/fft.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/output.hxx>
#include <bout/sys/timer.hxx>
#include <bout/utils.hxx>

#include <vector>

LaplaceCyclic::LaplaceCyclic(Options* opt, const CELL_LOC loc, Mesh* mesh_in,
                             Solver* UNUSED(solver))
    : Laplacian(opt, loc, mesh_in), Acoef(0.0), C1coef(1.0), C2coef(1.0), Dcoef(1.0) {

  Acoef.setLocation(location);
  C1coef.setLocation(location);
  C2coef.setLocation(location);
  Dcoef.setLocation(location);

  // Get options

  dst = (*opt)["dst"]
            .doc("Use Discrete Sine Transform in Z to enforce Dirichlet boundaries in Z")
            .withDefault<bool>(false);

  if (dst) {
    nmode = localmesh->LocalNz - 2;
  } else {
    // Number of Z modes. maxmode set in invert_laplace.cxx from options
    nmode = maxmode + 1;
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

  int jy = rhs.getIndex(); // Get the Y index
  x.setIndex(jy);

  // Get the width of the boundary

  // If the flags to assign that only one guard cell should be used is set
  int inbndry = localmesh->xstart, outbndry = localmesh->xstart;
  if (isGlobalFlagSet(INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2)) {
    inbndry = outbndry = 1;
  }
  if (isInnerBoundaryFlagSet(INVERT_BNDRY_ONE)) {
    inbndry = 1;
  }
  if (isOuterBoundaryFlagSet(INVERT_BNDRY_ONE)) {
    outbndry = 1;
  }

  if (dst) {
    BOUT_OMP_PERF(parallel)
    {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>(
          localmesh->LocalNz); // ZFFT routine expects input of this length

      // Loop over X indices, including boundaries but not guard cells. (unless periodic
      // in x)
      BOUT_OMP_PERF(for)
      for (int ix = xs; ix <= xe; ix++) {
        // Take DST in Z direction and put result in k1d

        if (((ix < inbndry) && isInnerBoundaryFlagSetOnFirstX(INVERT_SET))
            || ((localmesh->LocalNx - ix - 1 < outbndry)
                && isOuterBoundaryFlagSetOnLastX(INVERT_SET))) {
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
      BoutReal zlen = getUniform(coords->dz) * (localmesh->LocalNz - 3);
      BOUT_OMP_PERF(for nowait)
      for (int kz = 0; kz < nmode; kz++) {
        // wave number is 1/[rad]; DST has extra 2.
        BoutReal kwave = kz * 2.0 * PI / (2. * zlen);

        tridagMatrix(&a(kz, 0), &b(kz, 0), &c(kz, 0), &bcmplx(kz, 0), jy,
                     kz,    // wave number index
                     kwave, // kwave (inverse wave length)
                     &Acoef, &C1coef, &C2coef, &Dcoef,
                     false,  // Don't include guard cells in arrays
                     false); // Z domain not periodic
      }
    }

    // Solve tridiagonal systems
    cr->setCoefs(a, b, c);
    cr->solve(bcmplx, xcmplx);

    // FFT back to real space
    BOUT_OMP_PERF(parallel)
    {
      /// Create a local thread-scope working array

      // ZFFT routine expects input of this length
      auto k1d = Array<dcomplex>(localmesh->LocalNz);

      BOUT_OMP_PERF(for nowait)
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
    const BoutReal zlength = getUniform(coords->zlength());
    BOUT_OMP_PERF(parallel)
    {
      /// Create a local thread-scope working array
      // ZFFT routine expects input of this length
      auto k1d = Array<dcomplex>((localmesh->LocalNz) / 2 + 1);

      // Loop over X indices, including boundaries but not guard
      // cells (unless periodic in x)
      BOUT_OMP_PERF(for)
      for (int ix = xs; ix <= xe; ix++) {
        // Take FFT in Z direction, apply shift, and put result in k1d

        if (((ix < inbndry) && isInnerBoundaryFlagSetOnFirstX(INVERT_SET))
            || ((localmesh->LocalNx - ix - 1 < outbndry)
                && isOuterBoundaryFlagSetOnLastX(INVERT_SET))) {
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
      BOUT_OMP_PERF(for nowait)
      for (int kz = 0; kz < nmode; kz++) {
        BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
        tridagMatrix(&a(kz, 0), &b(kz, 0), &c(kz, 0), &bcmplx(kz, 0), jy,
                     kz,    // True for the component constant (DC) in Z
                     kwave, // Z wave number
                     &Acoef, &C1coef, &C2coef, &Dcoef,
                     false); // Don't include guard cells in arrays
      }
    }

    // Solve tridiagonal systems
    cr->setCoefs(a, b, c);
    cr->solve(bcmplx, xcmplx);

    if (localmesh->periodicX) {
      // Subtract X average of kz=0 mode
      BoutReal local[2] = {
          0.0,                               // index 0 = sum of coefficients
          static_cast<BoutReal>(xe - xs + 1) // number of grid cells
      };
      for (int ix = xs; ix <= xe; ix++) {
        local[0] += xcmplx(0, ix - xs).real();
      }
      BoutReal global[2];
      MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, localmesh->getXcomm());
      BoutReal avg = global[0] / global[1];
      for (int ix = xs; ix <= xe; ix++) {
        xcmplx(0, ix - xs) -= avg;
      }
    }

    // FFT back to real space
    BOUT_OMP_PERF(parallel)
    {
      /// Create a local thread-scope working array
      // ZFFT routine expects input of this length
      auto k1d = Array<dcomplex>((localmesh->LocalNz) / 2 + 1);

      const bool zero_DC = isGlobalFlagSet(INVERT_ZERO_DC);

      BOUT_OMP_PERF(for nowait)
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
  if (isGlobalFlagSet(INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2)) {
    inbndry = outbndry = 1;
  }
  if (isInnerBoundaryFlagSet(INVERT_BNDRY_ONE)) {
    inbndry = 1;
  }
  if (isOuterBoundaryFlagSet(INVERT_BNDRY_ONE)) {
    outbndry = 1;
  }

  int nx = xe - xs + 1; // Number of X points on this processor

  // Get range of Y indices
  int ys = localmesh->ystart, ye = localmesh->yend;

  if (localmesh->hasBndryLowerY()) {
    if (include_yguards) {
      ys = 0; // Mesh contains a lower boundary and we are solving in the guard cells
    }

    ys += extra_yguards_lower;
  }
  if (localmesh->hasBndryUpperY()) {
    if (include_yguards) {
      // Contains upper boundary and we are solving in the guard cells
      ye = localmesh->LocalNy - 1;
    }
    ye -= extra_yguards_upper;
  }

  const int ny = (ye - ys + 1); // Number of Y points
  const int nsys = nmode * ny;  // Number of systems of equations to solve
  const int nxny = nx * ny;     // Number of points in X-Y

  // This is just to silence static analysis
  ASSERT0(ny > 0);

  auto a3D = Matrix<dcomplex>(nsys, nx);
  auto b3D = Matrix<dcomplex>(nsys, nx);
  auto c3D = Matrix<dcomplex>(nsys, nx);

  auto xcmplx3D = Matrix<dcomplex>(nsys, nx);
  auto bcmplx3D = Matrix<dcomplex>(nsys, nx);

  if (dst) {
    BOUT_OMP_PERF(parallel)
    {
      /// Create a local thread-scope working array
      // ZFFT routine expects input of this length
      auto k1d = Array<dcomplex>(localmesh->LocalNz);

      // Loop over X and Y indices, including boundaries but not guard cells.
      // (unless periodic in x)
      BOUT_OMP_PERF(for)
      for (int ind = 0; ind < nxny; ++ind) {
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        // Take DST in Z direction and put result in k1d

        if (((ix < inbndry) && isInnerBoundaryFlagSetOnFirstX(INVERT_SET))
            || ((localmesh->LocalNx - ix - 1 < outbndry)
                && isOuterBoundaryFlagSetOnLastX(INVERT_SET))) {
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
      const BoutReal zlen = getUniform(coords->dz) * (localmesh->LocalNz - 3);
      BOUT_OMP_PERF(for nowait)
      for (int ind = 0; ind < nsys; ind++) {
        // ind = (iy - ys) * nmode + kz
        int iy = ys + ind / nmode;
        int kz = ind % nmode;

        // wave number is 1/[rad]; DST has extra 2.
        BoutReal kwave = kz * 2.0 * PI / (2. * zlen);

        tridagMatrix(&a3D(ind, 0), &b3D(ind, 0), &c3D(ind, 0), &bcmplx3D(ind, 0), iy,
                     kz,    // wave number index
                     kwave, // kwave (inverse wave length)
                     &Acoef, &C1coef, &C2coef, &Dcoef,
                     false,  // Don't include guard cells in arrays
                     false); // Z domain not periodic
      }
    }

    // Solve tridiagonal systems
    cr->setCoefs(a3D, b3D, c3D);
    cr->solve(bcmplx3D, xcmplx3D);

    // FFT back to real space
    BOUT_OMP_PERF(parallel)
    {
      /// Create a local thread-scope working array
      // ZFFT routine expects input of length LocalNz
      auto k1d = Array<dcomplex>(localmesh->LocalNz);

      BOUT_OMP_PERF(for nowait)
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
    const BoutReal zlength = getUniform(coords->zlength());
    BOUT_OMP_PERF(parallel)
    {
      /// Create a local thread-scope working array
      // ZFFT routine expects input of this length
      auto k1d = Array<dcomplex>(localmesh->LocalNz / 2 + 1);

      // Loop over X and Y indices, including boundaries but not guard cells
      // (unless periodic in x)

      BOUT_OMP_PERF(for)
      for (int ind = 0; ind < nxny; ++ind) {
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
        int ix = xs + ind / ny;
        int iy = ys + ind % ny;

        // Take FFT in Z direction, apply shift, and put result in k1d

        if (((ix < inbndry) && isInnerBoundaryFlagSetOnFirstX(INVERT_SET))
            || ((localmesh->LocalNx - ix - 1 < outbndry)
                && isOuterBoundaryFlagSetOnLastX(INVERT_SET))) {
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
      BOUT_OMP_PERF(for nowait)
      for (int ind = 0; ind < nsys; ind++) {
        // ind = (iy - ys) * nmode + kz
        int iy = ys + ind / nmode;
        int kz = ind % nmode;

        BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
        tridagMatrix(&a3D(ind, 0), &b3D(ind, 0), &c3D(ind, 0), &bcmplx3D(ind, 0), iy,
                     kz,    // True for the component constant (DC) in Z
                     kwave, // Z wave number
                     &Acoef, &C1coef, &C2coef, &Dcoef,
                     false); // Don't include guard cells in arrays
      }
    }

    // Solve tridiagonal systems
    cr->setCoefs(a3D, b3D, c3D);
    cr->solve(bcmplx3D, xcmplx3D);

    if (localmesh->periodicX) {
      // Subtract X average of kz=0 mode
      std::vector<BoutReal> local(ny + 1, 0.0);
      for (int y = 0; y < ny; y++) {
        for (int ix = xs; ix <= xe; ix++) {
          local[y] += xcmplx3D(y * nmode, ix - xs).real();
        }
      }
      local[ny] = static_cast<BoutReal>(xe - xs + 1);

      // Global reduce
      std::vector<BoutReal> global(ny + 1, 0.0);
      MPI_Allreduce(local.data(), global.data(), ny + 1, MPI_DOUBLE, MPI_SUM,
                    localmesh->getXcomm());
      // Subtract average from kz=0 modes
      for (int y = 0; y < ny; y++) {
        BoutReal avg = global[y] / global[ny];
        for (int ix = xs; ix <= xe; ix++) {
          xcmplx3D(y * nmode, ix - xs) -= avg;
        }
      }
    }

    // FFT back to real space
    BOUT_OMP_PERF(parallel)
    {
      /// Create a local thread-scope working array
      auto k1d = Array<dcomplex>((localmesh->LocalNz) / 2
                                 + 1); // ZFFT routine expects input of this length

      const bool zero_DC = isGlobalFlagSet(INVERT_ZERO_DC);

      BOUT_OMP_PERF(for nowait)
      for (int ind = 0; ind < nxny; ++ind) { // Loop over X and Y
        // ind = (ix - xs)*(ye - ys + 1) + (iy - ys)
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

void LaplaceCyclic ::verify_solution(const Matrix<dcomplex>& a_ver,
                                     const Matrix<dcomplex>& b_ver,
                                     const Matrix<dcomplex>& c_ver,
                                     const Matrix<dcomplex>& r_ver,
                                     const Matrix<dcomplex>& x_sol, const int nsys) {
  output.write("Verify solution\n");
  const int nx = xe - xs + 1; // Number of X points on this processor,
                              // including boundaries but not guard cells
  const int xproc = localmesh->getXProcIndex();
  const int yproc = localmesh->getYProcIndex();
  const int nprocs = localmesh->getNXPE();
  const int myrank = yproc * nprocs + xproc;
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

#endif // BOUT_USE_METRIC_3D
