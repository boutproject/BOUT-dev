/**************************************************************************
 * Perpendicular Laplacian inversion. Serial code using FFT
 * and tridiagonal solver.
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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

#include "globals.hxx"
#include "serial_tri.hxx"

#include <bout/mesh.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <lapack_routines.hxx>
#include <bout/constants.hxx>
#include <bout/openmpwrap.hxx>
#include <cmath>

#include <output.hxx>

LaplaceSerialTri::LaplaceSerialTri(Options *opt, CELL_LOC loc, Mesh *mesh_in)
    : Laplacian(opt, loc, mesh_in), A(0.0), C(1.0), D(1.0) {
  A.setLocation(location);
  C.setLocation(location);
  D.setLocation(location);

  if(!localmesh->firstX() || !localmesh->lastX()) {
    throw BoutException("LaplaceSerialTri only works for localmesh->NXPE = 1");
  }
}

FieldPerp LaplaceSerialTri::solve(const FieldPerp& b) { return solve(b, b); }

/*!
 * Solve Ax=b for x given b
 *
 * This function will
 *      1. Take the fourier transform of the y-slice given in the input
 *      2. For each fourier mode
 *          a) Set up the tridiagonal matrix
 *          b) Call the solver which inverts the matrix Ax_mode = b_mode
 *      3. Collect all the modes in a 2D array
 *      4. Back transform the y-slice
 *
 * Input:
 * \param[in] b     A 2D variable that will be fourier decomposed, each fourier
 *                  mode of this variable is going to be the right hand side of
 *                  the equation Ax = b
 * \param[in] x0    Variable used to set BC (if the right flags are set, see
 *                  the user manual)
 *
 * \return          The inverted variable.
 */
FieldPerp LaplaceSerialTri::solve(const FieldPerp& b, const FieldPerp& x0) {
  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  FieldPerp x{emptyFrom(b)};

  int jy = b.getIndex();

  int ncz = localmesh->LocalNz; // No of z pnts
  int ncx = localmesh->LocalNx; // No of x pnts

  BoutReal kwaveFactor = 2.0 * PI / coords->zlength();

  // Setting the width of the boundary.
  // NOTE: The default is a width of 2 guard cells
  int inbndry = localmesh->xstart, outbndry=localmesh->xstart;

  // If the flags to assign that only one guard cell should be used is set
  if((global_flags & INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2))  {
    inbndry = outbndry = 1;
  }
  if (inner_boundary_flags & INVERT_BNDRY_ONE)
    inbndry = 1;
  if (outer_boundary_flags & INVERT_BNDRY_ONE)
    outbndry = 1;

  /* Allocation fo
   * bk   = The fourier transformed of b, where b is one of the inputs in
   *        LaplaceSerialTri::solve()
   * bk1d = The 1d array of bk
   * xk   = The fourier transformed of x, where x the output of
   *        LaplaceSerialTri::solve()
   * xk1d = The 1d array of xk
   */
  auto bk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto bk1d = Array<dcomplex>(ncx);
  auto xk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto xk1d = Array<dcomplex>(ncx);

  // Initialise xk to 0 as we only visit 0<= kz <= maxmode in solve
  for (int ix = 0; ix < ncx; ix++) {
    for (int kz = maxmode + 1; kz < ncz / 2 + 1; kz++) {
      xk(ix, kz) = 0.0;
    }
  }

  /* Coefficents in the tridiagonal solver matrix
   * Following the notation in "Numerical recipes"
   * avec is the lower diagonal of the matrix
   * bvec is the diagonal of the matrix
   * cvec is the upper diagonal of the matrix
   * NOTE: Do not confuse avec, bvec and cvec with the A, C, and D coefficients
   *       above
   */
  auto avec = Array<dcomplex>(ncx);
  auto bvec = Array<dcomplex>(ncx);
  auto cvec = Array<dcomplex>(ncx);

  BOUT_OMP(parallel for)
  for (int ix = 0; ix < ncx; ix++) {
    /* This for loop will set the bk (initialized by the constructor)
     * bk is the z fourier modes of b in z
     * If the INVERT_SET flag is set (meaning that x0 will be used to set the
     * bounadry values),
     */
    if (((ix < inbndry) && (inner_boundary_flags & INVERT_SET)) ||
        ((ncx - 1 - ix < outbndry) && (outer_boundary_flags & INVERT_SET))) {
      // Use the values in x0 in the boundary

      // x0 is the input
      // bk is the output
      rfft(x0[ix], ncz, &bk(ix, 0));

    } else {
      // b is the input
      // bk is the output
      rfft(b[ix], ncz, &bk(ix, 0));
    }
  }

  /* Solve differential equation in x for each fourier mode
   * Note that only the non-degenerate fourier modes are being used (i.e. the
   * offset and all the modes up to the Nyquist frequency)
   */
  for (int kz = 0; kz <= maxmode; kz++) {

    // set bk1d
    for (int ix = 0; ix < ncx; ix++) {
      // Get bk of the current fourier mode
      bk1d[ix] = bk(ix, kz);
    }

    /* Set the matrix A used in the inversion of Ax=b
     * by calling tridagCoef and setting the BC
     *
     * Note that A, C and D in
     *
     * D*Laplace_perp(x) + (1/C)Grad_perp(C)*Grad_perp(x) + Ax = B
     *
     * has nothing to do with
     * avec - the lower diagonal of the tridiagonal matrix
     * bvec - the main diagonal
     * cvec - the upper diagonal
    */
    tridagMatrix(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(bk1d),
                 jy,
                 // wave number index
                 kz,
                 // wave number (different from kz only if we are taking a part
                 // of the z-domain [and not from 0 to 2*pi])
                 kz * kwaveFactor, global_flags, inner_boundary_flags,
                 outer_boundary_flags, &A, &C, &D);

    ///////// PERFORM INVERSION /////////
    if (!localmesh->periodicX) {
      // Call tridiagonal solver
      tridag(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(bk1d),
             std::begin(xk1d), ncx);

    } else {
      // Periodic in X, so cyclic tridiagonal

      int xs = localmesh->xstart;
      cyclic_tridag(&avec[xs], &bvec[xs], &cvec[xs], &bk1d[xs], &xk1d[xs], ncx - 2 * xs);

      // Copy boundary regions
      for (int ix = 0; ix < xs; ix++) {
        xk1d[ix] = xk1d[ncx - 2 * xs + ix];
        xk1d[ncx - xs + ix] = xk1d[xs + ix];
      }
    }

    // If the global flag is set to INVERT_KX_ZERO
    if ((global_flags & INVERT_KX_ZERO) && (kz == 0)) {
      dcomplex offset(0.0);
      for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
        offset += xk1d[ix];
      }
      offset /= static_cast<BoutReal>(localmesh->xend - localmesh->xstart + 1);
      for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
        xk1d[ix] -= offset;
      }
    }

    // Store the solution xk for the current fourier mode in a 2D array
    for (int ix = 0; ix < ncx; ix++) {
      xk(ix, kz) = xk1d[ix];
    }
  }

  // Done inversion, transform back
  for (int ix = 0; ix < ncx; ix++) {

    if(global_flags & INVERT_ZERO_DC)
      xk(ix, 0) = 0.0;

    irfft(&xk(ix, 0), ncz, x[ix]);

#if CHECK > 2
    for(int kz=0;kz<ncz;kz++)
      if(!finite(x(ix,kz)))
        throw BoutException("Non-finite at {:d}, {:d}, {:d}", ix, jy, kz);
#endif
  }

  checkData(x);

  return x; // Result of the inversion
}
