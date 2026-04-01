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

#include "../../common_transform.hxx"

#include "bout/array.hxx"
#include "bout/dcomplex.hxx"
#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/fft.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/output.hxx>
#include <bout/sys/timer.hxx>
#include <bout/utils.hxx>

#include <array>
#include <vector>

LaplaceCyclic::LaplaceCyclic(Options* opt, const CELL_LOC loc, Mesh* mesh_in,
                             Solver* UNUSED(solver))
    : Laplacian(opt, loc, mesh_in), Acoef(0.0), C1coef(1.0), C2coef(1.0), Dcoef(1.0) {

  bout::fft::assertZSerial(*localmesh, "`cyclic` inversion");

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

  // Get the width of the boundary

  // If the flags to assign that only one guard cell should be used is set
  int inbndry = localmesh->xstart;
  int outbndry = localmesh->xstart;
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
    const DSTTransform transform(
        *localmesh, nmode, xs, xe, 0, 0, localmesh->zstart, localmesh->zend, inbndry,
        outbndry, isInnerBoundaryFlagSetOnFirstX(INVERT_SET),
        isOuterBoundaryFlagSetOnLastX(INVERT_SET), isGlobalFlagSet(INVERT_ZERO_DC));

    auto matrices = transform.forward(*this, rhs, x0, Acoef, C1coef, C2coef, Dcoef);

    // Solve tridiagonal systems
    cr->setCoefs(a, b, c);
    cr->solve(bcmplx, xcmplx);

    return transform.backward(rhs, xcmplx);
  }

  const FFTTransform transform(
      *localmesh, nmode, xs, xe, 0, 0, localmesh->zstart, localmesh->zend, inbndry,
      outbndry, isInnerBoundaryFlagSetOnFirstX(INVERT_SET),
      isOuterBoundaryFlagSetOnLastX(INVERT_SET), isGlobalFlagSet(INVERT_ZERO_DC));

  auto matrices = transform.forward(*this, rhs, x0, Acoef, C1coef, C2coef, Dcoef);

  // Solve tridiagonal systems
  cr->setCoefs(matrices.a, matrices.b, matrices.c);
  cr->solve(matrices.bcmplx, xcmplx);

  if (localmesh->periodicX) {
    // Subtract X average of kz=0 mode
    std::array<BoutReal, 2> local = {
        0.0,                               // index 0 = sum of coefficients
        static_cast<BoutReal>(xe - xs + 1) // number of grid cells
    };
    for (int ix = xs; ix <= xe; ix++) {
      local[0] += xcmplx(0, ix - xs).real();
    }
    std::array<BoutReal, 2> global{};
    MPI_Allreduce(local.data(), global.data(), 2, MPI_DOUBLE, MPI_SUM,
                  localmesh->getXcomm());
    const BoutReal avg = global[0] / global[1];
    for (int ix = xs; ix <= xe; ix++) {
      xcmplx(0, ix - xs) -= avg;
    }
  }

  return transform.backward(rhs, xcmplx);
}

Field3D LaplaceCyclic::solve(const Field3D& rhs, const Field3D& x0) {

  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  ASSERT1(localmesh == rhs.getMesh() && localmesh == x0.getMesh());

  Timer timer("invert");

  Field3D x{emptyFrom(rhs)}; // Result

  // Get the width of the boundary

  // If the flags to assign that only one guard cell should be used is set
  int inbndry = localmesh->xstart;
  int outbndry = localmesh->xstart;
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
      // Contains upper boundary and we are solving in the guard cells
      ye = localmesh->LocalNy - 1;
    }
    ye -= extra_yguards_upper;
  }

  const int ny = (ye - ys + 1); // Number of Y points
  const int nsys = nmode * ny;  // Number of systems of equations to solve

  // This is just to silence static analysis
  ASSERT0(ny > 0);

  auto xcmplx3D = Matrix<dcomplex>(nsys, nx);

  if (dst) {
    const DSTTransform transform(
        *localmesh, nmode, xs, xe, ys, ye, localmesh->zstart, localmesh->zend, inbndry,
        outbndry, isInnerBoundaryFlagSetOnFirstX(INVERT_SET),
        isOuterBoundaryFlagSetOnLastX(INVERT_SET), isGlobalFlagSet(INVERT_ZERO_DC));

    auto matrices = transform.forward(*this, rhs, x0, Acoef, C1coef, C2coef, Dcoef);

    // Solve tridiagonal systems
    cr->setCoefs(matrices.a, matrices.b, matrices.c);
    cr->solve(matrices.bcmplx, xcmplx3D);

    return transform.backward(rhs, xcmplx3D);
  }
  const FFTTransform transform(
      *localmesh, nmode, xs, xe, ys, ye, localmesh->zstart, localmesh->zend, inbndry,
      outbndry, isInnerBoundaryFlagSetOnFirstX(INVERT_SET),
      isOuterBoundaryFlagSetOnLastX(INVERT_SET), isGlobalFlagSet(INVERT_ZERO_DC));

  auto matrices = transform.forward(*this, rhs, x0, Acoef, C1coef, C2coef, Dcoef);

  // Solve tridiagonal systems
  cr->setCoefs(matrices.a, matrices.b, matrices.c);
  cr->solve(matrices.bcmplx, xcmplx3D);

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

  return transform.backward(rhs, xcmplx3D);
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
    MPI_Irecv(&rbufdown[0], nsys, MPI_DOUBLE_COMPLEX, myrank - 1, 901, BoutComm::get(),
              &request[1]);
    for (int kz = 0; kz < nsys; kz++) {
      sbufdown[kz] = x_ver(kz, 1);
    }
    MPI_Isend(&sbufdown[0], nsys, MPI_DOUBLE_COMPLEX, myrank - 1, 900, BoutComm::get(),
              &request[0]);
  }
  if (xproc < nprocs - 1) {
    MPI_Irecv(&rbufup[0], nsys, MPI_DOUBLE_COMPLEX, myrank + 1, 900, BoutComm::get(),
              &request[3]);
    for (int kz = 0; kz < nsys; kz++) {
      sbufup[kz] = x_ver(kz, nx);
    }
    MPI_Isend(&sbufup[0], nsys, MPI_DOUBLE_COMPLEX, myrank + 1, 901, BoutComm::get(),
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
