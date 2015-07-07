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

  int ncz = mesh->ngz-1;

  bk = cmatrix(mesh->ngx, ncz/2 + 1);
  bk1d = new dcomplex[mesh->ngx];
  
  xk = cmatrix(mesh->ngx, ncz/2 + 1);
  xk1d = new dcomplex[mesh->ngx];
  
  avec = new dcomplex[mesh->ngx];
  bvec = new dcomplex[mesh->ngx];
  cvec = new dcomplex[mesh->ngx];
}

LaplaceSerialTri::~LaplaceSerialTri() {
  free_cmatrix(bk);
  delete[] bk1d;
  free_cmatrix(xk);
  delete[] xk1d;
  
  delete[] avec;
  delete[] bvec;
  delete[] cvec;
}

const FieldPerp LaplaceSerialTri::solve(const FieldPerp &b) {
  return solve(b,b);   // Call the solver below
}

const FieldPerp LaplaceSerialTri::solve(const FieldPerp &b, const FieldPerp &x0) {
  FieldPerp x;
  x.allocate();

  int jy = b.getIndex();
  x.setIndex(jy);
  
  int ncz = mesh->ngz-1;
  int ncx = mesh->ngx-1;
  
  // Used when checking if one should set the boundary
  int inbndry = 2, outbndry=2;
  
  if(global_flags & INVERT_BOTH_BNDRY_ONE) {
    inbndry = outbndry = 1;
  }
  if(inner_boundary_flags & INVERT_BNDRY_ONE)
    inbndry = 1;
  if(outer_boundary_flags & INVERT_BNDRY_ONE)
    outbndry = 1;
  
  #pragma omp parallel for
  for(int ix=0;ix<mesh->ngx;ix++) {
    // for fixed ix,jy set a complex vector rho(z)
    
    if(((ix < inbndry) && (inner_boundary_flags & INVERT_SET)) ||
       ((ncx-ix < outbndry) && (outer_boundary_flags & INVERT_SET))) {
      // Use the values in x0 in the boundary

      // x0 and mesh->zShift are the inputs
      // bk is the output
      ZFFT(x0[ix], mesh->zShift[ix][jy], bk[ix]);
      
    }else
      // b and mesh->zShift are the inputs
      // bk is the output
      ZFFT(b[ix], mesh->zShift[ix][jy], bk[ix]);
  }

  // Solve differential equation in x for each fourier mode of bk (up to half of max)
  for(int iz=0;iz<=ncz/2;iz++) {

    // set bk1d
    BoutReal flt;
    if (iz>maxmode) flt=0.0; else flt=1.0;
      
    for(int ix=0;ix<=ncx;ix++)
      // Get bk of the current fourier mode
      bk1d[ix] = bk[ix][iz] * flt;
    
    // Sets the lower diagonal avec, the diagonal bvec and the upper
    // diagonal cvec in the tridiagonal matrix (see "Numerical recipes" for 
    // notation) with the correct numbers (including the metrics) by calling
    // tridagCoef
    // tridagMatrix also sets the boundary (by setting avec, bvec, cvec and bk
    // for the boundary points) according to the given flags.    
    tridagMatrix(avec, bvec, cvec, bk1d, jy, 
		 iz == 0, // Called "int kz" in the function. Used to check if DC.
		 iz*2.0*PI/mesh->zlength, // kwave
		 global_flags, inner_boundary_flags, outer_boundary_flags,
		 &A, &C, &D);

    ///////// PERFORM INVERSION /////////
    if(!mesh->periodicX) {
      // Call tridiagonal solver
      tridag(avec, bvec, cvec, bk1d, xk1d, mesh->ngx);
      
    } else {
      // Periodic in X, so cyclic tridiagonal
      cyclic_tridag(avec+2, bvec+2, cvec+2, bk1d+2, xk1d+2, mesh->ngx-4);
	
      // Copy boundary regions
      for(int ix=0;ix<2;ix++) {
        xk1d[ix] = xk1d[mesh->ngx-4+ix];
        xk1d[mesh->ngx-2+ix] = xk1d[2+ix];
      }
    }
      
    if((global_flags & INVERT_KX_ZERO) && (iz == 0)) {
      dcomplex offset(0.0);
      for(int ix=0;ix<=ncx;ix++)
        offset += bk1d[ix];
      offset /= (BoutReal) (ncx+1);
      for(int ix=0;ix<=ncx;ix++)
        bk1d[ix] -= offset;
    }
      
    // Set the answer xk for the current fourier mode
    for (int ix=0; ix<=ncx; ix++){
      xk[ix][iz]=xk1d[ix];
    }
  }

  // Done inversion, transform back
  for(int ix=0; ix<=ncx; ix++){
    
    if(global_flags & INVERT_ZERO_DC)
      xk[ix][0] = 0.0;
    
    ZFFT_rev(xk[ix], mesh->zShift[ix][jy], x[ix]);
    
    x[ix][mesh->ngz-1] = x[ix][0]; // enforce periodicity
    
    for(int kz=0;kz<mesh->ngz;kz++)
      if(!finite(x[ix][kz]))
        throw BoutException("Non-finite at %d, %d, %d", ix, jy, kz);
  }

  return x;
}
