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

#include <boutexception.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <lapack_routines.hxx>
#include <bout/constants.hxx>
#include <cmath>

#include <output.hxx>

LaplaceSerialTri::LaplaceSerialTri(Options *opt) : Laplacian(opt), A(0.0), C(1.0), D(1.0) {

  if(!mesh->firstX() || !mesh->lastX()) {
    throw BoutException("LaplaceSerialTri only works for mesh->NXPE = 1");
  }

  // Allocate memory

  int ncz = mesh->LocalNz;

  bk = matrix<dcomplex>(mesh->LocalNx, ncz/2 + 1);
  bk1d = new dcomplex[mesh->LocalNx];

  //Initialise bk to 0 as we only visit 0<= kz <= maxmode in solve
  for(int kz=maxmode+1; kz < ncz/2 + 1; kz++){
    for (int ix=0; ix<mesh->LocalNx; ix++){
      bk[ix][kz] = 0.0;
    }
  }

  xk = matrix<dcomplex>(mesh->LocalNx, ncz/2 + 1);
  xk1d = new dcomplex[mesh->LocalNx];

  //Initialise xk to 0 as we only visit 0<= kz <= maxmode in solve
  for(int kz=maxmode+1; kz < ncz/2 + 1; kz++){
    for (int ix=0; ix<mesh->LocalNx; ix++){
      xk[ix][kz] = 0.0;
    }
  }

  avec = new dcomplex[mesh->LocalNx];
  bvec = new dcomplex[mesh->LocalNx];
  cvec = new dcomplex[mesh->LocalNx];
}

LaplaceSerialTri::~LaplaceSerialTri() {
  free_matrix(bk);
  delete[] bk1d;
  free_matrix(xk);
  delete[] xk1d;

  delete[] avec;
  delete[] bvec;
  delete[] cvec;
}

const FieldPerp LaplaceSerialTri::solve(const FieldPerp &b) {
  return solve(b,b);   // Call the solver below
}

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
 * \param[out] x    The inverted variable.
 */
const FieldPerp LaplaceSerialTri::solve(const FieldPerp &b, const FieldPerp &x0) {
  Mesh * mesh = b.getMesh();
  FieldPerp x(mesh);
  x.allocate();

  Coordinates *coord = mesh->coordinates();

  int jy = b.getIndex();
  x.setIndex(jy);

  int ncz = mesh->LocalNz; // No of z pnts (counts from 1 to easily convert to kz)
  int ncx = mesh->LocalNx-1; // No of x pnts (counts from 0)

  // Setting the width of the boundary.
  // NOTE: The default is a width of 2 guard cells
  int inbndry = 2, outbndry=2;

  // If the flags to assign that only one guard cell should be used is set
  if(global_flags & INVERT_BOTH_BNDRY_ONE) {
    inbndry = outbndry = 1;
  }
  if(inner_boundary_flags & INVERT_BNDRY_ONE)
    inbndry = 1;
  if(outer_boundary_flags & INVERT_BNDRY_ONE)
    outbndry = 1;

  #pragma omp parallel for
  for(int ix=0;ix<mesh->LocalNx;ix++) {
    /* This for loop will set the bk (initialized by the constructor)
     * bk is the z fourier modes of b in z
     * If the INVERT_SET flag is set (meaning that x0 will be used to set the
     * bounadry values),
     */
    if(((ix < inbndry) && (inner_boundary_flags & INVERT_SET)) ||
       ((ncx-ix < outbndry) && (outer_boundary_flags & INVERT_SET))) {
      // Use the values in x0 in the boundary

      // x0 is the input
      // bk is the output
      rfft(x0[ix], ncz, bk[ix]);

    }else {
      // b is the input
      // bk is the output
      rfft(b[ix], ncz, bk[ix]);
    }
  }

  /* Solve differential equation in x for each fourier mode
   * Note that only the non-degenerate fourier modes are being used (i.e. the
   * offset and all the modes up to the Nyquist frequency)
   */
  for(int kz=0;kz<=maxmode;kz++) {

    // set bk1d
    for(int ix=0;ix<=ncx;ix++)
      // Get bk of the current fourier mode
      bk1d[ix] = bk[ix][kz];

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
    tridagMatrix(avec, bvec, cvec, bk1d, jy,
                 // wave number index
                 kz,
                 // wave number (different from kz only if we are taking a part
                 // of the z-domain [and not from 0 to 2*pi])
                 kz*2.0*PI/coord->zlength(),
                 global_flags, inner_boundary_flags, outer_boundary_flags,
                 &A, &C, &D);

    ///////// PERFORM INVERSION /////////
    if(!mesh->periodicX) {
      // Call tridiagonal solver
      tridag(avec, bvec, cvec, bk1d, xk1d, mesh->LocalNx);

    } else {
      // Periodic in X, so cyclic tridiagonal
      cyclic_tridag(avec+2, bvec+2, cvec+2, bk1d+2, xk1d+2, mesh->LocalNx-4);

      // Copy boundary regions
      for(int ix=0;ix<2;ix++) {
        xk1d[ix] = xk1d[mesh->LocalNx-4+ix];
        xk1d[mesh->LocalNx-2+ix] = xk1d[2+ix];
      }
    }

    // If the global flag is set to INVERT_KX_ZERO
    if((global_flags & INVERT_KX_ZERO) && (kz == 0)) {
      dcomplex offset(0.0);
      for(int ix=0;ix<=ncx;ix++)
        offset += bk1d[ix];
      offset /= static_cast<BoutReal>(ncx + 1);
      for(int ix=0;ix<=ncx;ix++)
        bk1d[ix] -= offset;
    }

    // Store the solution xk for the current fourier mode in a 2D array
    for (int ix=0; ix<=ncx; ix++){
      xk[ix][kz]=xk1d[ix];
    }
  }

  // Done inversion, transform back
  for(int ix=0; ix<=ncx; ix++){

    if(global_flags & INVERT_ZERO_DC)
      xk[ix][0] = 0.0;

    irfft(xk[ix], ncz, x[ix]);

#if CHECK > 2
    for(int kz=0;kz<ncz;kz++)
      if(!finite(x(ix,kz)))
        throw BoutException("Non-finite at %d, %d, %d", ix, jy, kz);
#endif
  }

  return x; // Result of the inversion
}
