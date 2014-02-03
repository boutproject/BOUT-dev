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

LaplaceCyclic::LaplaceCyclic(Options *opt) : Laplacian(opt), A(0.0), C(1.0), D(1.0) {
  // Get options
  
  OPTION(opt, dst, false);

  if(dst) {
    nmode = mesh->ngz-3;
  }else
    nmode = maxmode+1; // Number of Z modes. maxmode set in invert_laplace.cxx from options

  // Allocate arrays
  
  int nsys = nmode;     // Number of tridiagonal systems to solve simulaneously
  
  xs = mesh->xstart; // Starting X index
  if(mesh->firstX())
    xs = 0;
  xe = mesh->xend;  // Last X index
  if(mesh->lastX())
    xe = mesh->ngx-1;

  int n = xe - xs + 1;  // Number of X points on this processor, 
                        // including boundaries but not guard cells
  
  a   = matrix<dcomplex>(nsys, n);
  b   = matrix<dcomplex>(nsys, n);
  c   = matrix<dcomplex>(nsys, n);
  xcmplx = matrix<dcomplex>(nsys, n);
  bcmplx = matrix<dcomplex>(nsys, n);
  
  if(dst)
    k1d = new dcomplex[mesh->ngz-1];         // DST has different k space
  else
    k1d = new dcomplex[(mesh->ngz-1)/2 + 1]; // ZFFT routine expects input of this length
  
  // Create a cyclic reduction object, operating on dcomplex values
  cr = new CyclicReduce<dcomplex>(mesh->getXcomm(), n);
  cr->setPeriodic(mesh->periodicX);
}

LaplaceCyclic::~LaplaceCyclic() {
  // Free arrays
  
  free_matrix(a);
  free_matrix(b);
  free_matrix(c);
  free_matrix(xcmplx);
  free_matrix(bcmplx);
  
  delete[] k1d;
  
  // Delete tridiagonal solver
  delete cr;
}

const FieldPerp LaplaceCyclic::solve(const FieldPerp &rhs, const FieldPerp &x0) {
  FieldPerp x;  // Result
  x.allocate();

  int jy = rhs.getIndex();  // Get the Y index
  x.setIndex(jy);
  
  // Get the width of the boundary
  
  int inbndry = 2, outbndry=2;
  if(flags & INVERT_BNDRY_ONE) {
    inbndry = outbndry = 1;
  }
  if(flags & INVERT_BNDRY_IN_ONE)
    inbndry = 1;
  if(flags & INVERT_BNDRY_OUT_ONE)
    outbndry = 1;

  if(dst) {
    // Loop over X indices, including boundaries but not guard cells
    for(int ix=xs; ix <= xe; ix++) {
      // Take DST in Z direction and put result in k1d  

      if(((ix < inbndry) && (flags & INVERT_IN_SET) && mesh->firstX()) ||
         ((xe-ix < outbndry) && (flags & INVERT_OUT_SET) && mesh->lastX())) {
        // Use the values in x0 in the boundary
        DST(x0[ix]+1, mesh->ngz-3 , k1d);
      }else {
        //	  int length = sizeof(*rhs[ix])/sizeof(FieldPerp);
        DST(rhs[ix]+1, mesh->ngz-3 , k1d);
      }
    
      // Copy into array, transposing so kz is first index
      for(int kz = 0; kz < nmode; kz++)
        bcmplx[kz][ix-xs] = k1d[kz];
    }
  
    // Get elements of the tridiagonal matrix
    // including boundary conditions
    for(int kz = 0; kz < nmode; kz++) {
	
      BoutReal zlen = mesh->dz*(mesh->ngz-4);
	
      BoutReal kwave=kz*2.0*PI/(2.*zlen); // wave number is 1/[rad]; DST has extra 2.
    
      tridagMatrix(a[kz], b[kz], c[kz],
                   bcmplx[kz], 
                   jy, 
                   kz == 0, // True for the component constant (DC) in Z
                   kwave,   // Z wave number
                   flags, 
                   &A, &C, &D,
                   false);  // Don't include guard cells in arrays
    }
 
    // Solve tridiagonal systems

    cr->setCoefs(nmode, a, b, c);
    cr->solve(nmode, bcmplx, xcmplx);
  
    // FFT back to real space
    for(int ix=xs; ix <= xe; ix++) {
      for(int kz = 0; kz < nmode; kz++)
        k1d[kz] = xcmplx[kz][ix-xs];
	
      for(int kz=nmode;kz<(mesh->ngz-1);kz++)  
        k1d[kz] = 0.0; // Filtering out all higher harmonics

      DST_rev(k1d, mesh->ngz-3, x[ix]+1);
	
      x[ix][0] = -x[ix][2];
      x[ix][mesh->ngz-2] = -x[ix][mesh->ngz-4];
      x[ix][mesh->ngz-1] = 0.0; // probably unnecessary
    }
  }else {
    // Loop over X indices, including boundaries but not guard cells
    for(int ix=xs; ix <= xe; ix++) {
      // Take FFT in Z direction, apply shift, and put result in k1d
    
      if(((ix < inbndry) && (flags & INVERT_IN_SET) && mesh->firstX()) ||
         ((xe-ix < outbndry) && (flags & INVERT_OUT_SET) && mesh->lastX())) {
        // Use the values in x0 in the boundary
        ZFFT(x0[ix], mesh->zShift(ix, jy), k1d);
      }else {
        ZFFT(rhs[ix], mesh->zShift(ix, jy), k1d);
      }
    
      // Copy into array, transposing so kz is first index
      for(int kz = 0; kz < nmode; kz++)
        bcmplx[kz][ix-xs] = k1d[kz];
    }
  
    // Get elements of the tridiagonal matrix
    // including boundary conditions
    for(int kz = 0; kz < nmode; kz++) {
      BoutReal kwave=kz*2.0*PI/(mesh->zlength); // wave number is 1/[rad]
    
      tridagMatrix(a[kz], b[kz], c[kz],
                   bcmplx[kz], 
                   jy, 
                   kz == 0, // True for the component constant (DC) in Z
                   kwave,   // Z wave number
                   flags, 
                   &A, &C, &D,
                   false);  // Don't include guard cells in arrays
    }
  
    // Solve tridiagonal systems
  
    cr->setCoefs(nmode, a, b, c);
    cr->solve(nmode, bcmplx, xcmplx);
  
    // FFT back to real space
    for(int ix=xs; ix <= xe; ix++) {
      for(int kz = 0; kz < nmode; kz++)
        k1d[kz] = xcmplx[kz][ix-xs];
	
      for(int kz=nmode;kz<(mesh->ngz-1)/2 + 1;kz++)
        k1d[kz] = 0.0; // Filtering out all higher harmonics

      ZFFT_rev(k1d, mesh->zShift(ix, jy), x[ix]);
    
      x[ix][mesh->ngz-1] = x[ix][0]; // probably unnecessary
    }
  }
  return x;
}
