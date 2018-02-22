/**************************************************************************
 * Perpendicular Laplacian inversion
 * 
 * This code uses the Parallel Diagonally Dominant algorithm. This is very efficient
 * (constant number of communications), but achieves this by neglecting "small"
 * corrections. For ELM simulations these seem to be non-negligable, hence:
 *
 * CHECK IF THIS ALGORITHM PRODUCES REASONABLE RESULTS FOR YOUR PROBLEM
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

#include <bout/globals.hxx>
#include <bout/boutexception.hxx>
#include <bout/utils.hxx>
#include <bout/fft.hxx>
#include <bout/lapack_routines.hxx>

#include "pdd.hxx"

const FieldPerp LaplacePDD::solve(const FieldPerp &b) {
  static PDD_data data;
  static bool allocated = false;
  if(!allocated) {
    data.bk = NULL;
    allocated = true;
  }

  FieldPerp x(b.getMesh());
  x.allocate();
  
  start(b, data);
  next(data);
  finish(data, x);
  
  return x;
}

const Field3D LaplacePDD::solve(const Field3D &b) {
  Mesh *mesh = b.getMesh();
  Field3D x(mesh);
  x.allocate();
  FieldPerp xperp(mesh);
  xperp.allocate();
  
  int ys = mesh->ystart, ye = mesh->yend;
  if(mesh->hasBndryLowerY())
    ys = 0; // Mesh contains a lower boundary
  if(mesh->hasBndryUpperY())
    ye = mesh->LocalNy-1; // Contains upper boundary
  
  if(low_mem) {
    // Solve one slice at a time
    for(int jy=ys; jy <= ye; jy++) {
      x = solve(sliceXZ(b, jy));
    }
  }else {
    // Overlap multiple inversions
    
    static PDD_data *data = NULL;
    
    if(data == NULL) {
      data = new PDD_data[ye - ys + 1];
      data -= ys; // Re-number indices to start at jstart
      for(int jy=ys;jy<=ye;jy++)
        data[jy].bk = NULL; // Mark as unallocated for PDD routine
    }
    /// PDD algorithm communicates twice, so done in 3 stages
      
    for(int jy=ys; jy <= ye; jy++)
      start(sliceXZ(b,jy), data[jy]);
    
    for(int jy=ys; jy <= ye; jy++)
      next(data[jy]);
    
    for(int jy=ys; jy <= ye; jy++) {
      finish(data[jy], xperp);
      x = xperp;
    }
  }

  x.setLocation(b.getLocation()); 

  return x;
}

/// Laplacian inversion using Parallel Diagonal Dominant (PDD) method
/*!
 *
 * July 2008: Adapted from serial version to run in parallel (split in X) for tridiagonal system
 * i.e. no 4th order inversion yet.
 *
 * \note This code stores intermediate results and takes significantly more memory than
 * the serial version. This can be balanced against communication time i.e. faster communications
 * can allow less memory use.
 *
 * @param[in] data  Internal data used for multiple calls in parallel mode
 * @param[in] stage Which stage of the inversion, used to overlap calculation and communications.
 */
void LaplacePDD::start(const FieldPerp &b, PDD_data &data) {
  int ix, kz;
  Mesh *mesh = b.getMesh();

  int ncz = mesh->LocalNz;

  data.jy = b.getIndex();

  if(mesh->firstX() && mesh->lastX())
    throw BoutException("Error: PDD method only works for NXPE > 1\n");
  
  if(mesh->periodicX) {
      throw BoutException("LaplacePDD does not work with periodicity in the x direction (mesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
    }

  if(data.bk == NULL) {
    // Need to allocate working memory
    
    // RHS vector
    data.bk = matrix<dcomplex>(maxmode + 1, mesh->LocalNx);
    
    // Matrix to be solved
    data.avec = matrix<dcomplex>(maxmode + 1, mesh->LocalNx);
    data.bvec = matrix<dcomplex>(maxmode + 1, mesh->LocalNx);
    data.cvec = matrix<dcomplex>(maxmode + 1, mesh->LocalNx);
    
    // Working vectors
    data.v = matrix<dcomplex>(maxmode + 1, mesh->LocalNx);
    data.w = matrix<dcomplex>(maxmode + 1, mesh->LocalNx);

    // Result
    data.xk = matrix<dcomplex>(maxmode + 1, mesh->LocalNx);

    // Communication buffers. Space for 2 complex values for each kz
    data.snd = new BoutReal[4*(maxmode+1)];
    data.rcv = new BoutReal[4*(maxmode+1)];

    data.y2i = new dcomplex[maxmode + 1];
  }

  /// Take FFTs of data
  static dcomplex *bk1d = NULL; ///< 1D in Z for taking FFTs

  if(bk1d == NULL)
    bk1d = new dcomplex[ncz/2 + 1];

  for(ix=0; ix < mesh->LocalNx; ix++) {
    rfft(b[ix], ncz, bk1d);
    for(kz = 0; kz <= maxmode; kz++)
      data.bk[kz][ix] = bk1d[kz];
  }

  /// Create the matrices to be inverted (one for each z point)

  /// Set matrix elements
  tridagMatrix(data.avec, data.bvec, data.cvec, data.bk, data.jy, global_flags,
               inner_boundary_flags, outer_boundary_flags, &Acoef, &Ccoef, &Dcoef);

  for(kz = 0; kz <= maxmode; kz++) {
    // Start PDD algorithm

    // Solve for xtilde, v and w (step 2)

    static dcomplex *e = NULL;
    if(e == NULL) {
      e = new dcomplex[mesh->LocalNx];
      for(ix=0;ix<mesh->LocalNx;ix++)
	e[ix] = 0.0;
    }

    dcomplex v0, x0; // Values to be sent to processor i-1

    if(mesh->firstX()) {
      // Domain includes inner boundary
      tridag(data.avec[kz], data.bvec[kz], data.cvec[kz], 
	     data.bk[kz], data.xk[kz], mesh->xend+1);
      
      // Add C (row m-1) from next processor
      
      e[mesh->xend] = data.cvec[kz][mesh->xend];
      tridag(data.avec[kz], data.bvec[kz], data.cvec[kz], 
	     e, data.w[kz], mesh->xend+1);

    }else if(mesh->lastX()) {
      // Domain includes outer boundary
      tridag(data.avec[kz]+mesh->xstart, 
	     data.bvec[kz]+mesh->xstart, 
	     data.cvec[kz]+mesh->xstart, 
	     data.bk[kz]+mesh->xstart, 
	     data.xk[kz]+mesh->xstart, 
	     mesh->xend - mesh->xend + 1);
      
      // Add A (row 0) from previous processor
      e[0] = data.avec[kz][mesh->xstart];
      tridag(data.avec[kz]+mesh->xstart, 
	     data.bvec[kz]+mesh->xstart, 
	     data.cvec[kz]+mesh->xstart, 
	     e, data.v[kz]+mesh->xstart,
	     mesh->xend+1);
      
      x0 = data.xk[kz][mesh->xstart];
      v0 = data.v[kz][mesh->xstart];

    }else {
      // No boundaries
      tridag(data.avec[kz]+mesh->xstart,
	     data.bvec[kz]+mesh->xstart,
	     data.cvec[kz]+mesh->xstart, 
	     data.bk[kz]+mesh->xstart, 
	     data.xk[kz]+mesh->xstart, 
	     mesh->xend - mesh->xstart + 1);

      // Add A (row 0) from previous processor
      e[0] = data.avec[kz][mesh->xstart];
      tridag(data.avec[kz]+mesh->xstart,
	     data.bvec[kz]+mesh->xstart,
	     data.cvec[kz]+mesh->xstart, 
	     e+mesh->xstart,
	     data.v[kz]+mesh->xstart,
	     mesh->xend - mesh->xstart + 1);
      e[0] = 0.0;
      
      // Add C (row m-1) from next processor
      e[mesh->xend] = data.cvec[kz][mesh->xend];
      tridag(data.avec[kz]+mesh->xstart,
	     data.bvec[kz]+mesh->xstart,
	     data.cvec[kz]+mesh->xstart, 
	     e+mesh->xstart,
	     data.v[kz]+mesh->xstart, 
	     mesh->xend - mesh->xstart + 1);
      e[mesh->xend] = 0.0;
    }
    
    // Put values into communication buffers
    data.snd[4*kz]   = x0.real();
    data.snd[4*kz+1] = x0.imag();
    data.snd[4*kz+2] = v0.real();
    data.snd[4*kz+3] = v0.imag();
  }
  
  // Stage 3: Communicate x0, v0 from node i to i-1
  
  if(!mesh->lastX()) {
    // All except the last processor expect to receive data
    // Post async receive
    data.recv_handle = mesh->irecvXOut(data.rcv, 4*(maxmode+1), PDD_COMM_XV);
  }

  if(!mesh->firstX()) {
    // Send the data
    
    mesh->sendXIn(data.snd, 4*(maxmode+1), PDD_COMM_XV);
  }
}


/// Middle part of the PDD algorithm
void LaplacePDD::next(PDD_data &data) {
  // Wait for x0 and v0 to arrive from processor i+1
  
  if(!mesh->lastX()) {
    mesh->wait(data.recv_handle);

    /*! Now solving on all except the last processor
     * 
     * |    1       w^(i)_(m-1) | | y_{2i}   | = | x^(i)_{m-1} |
     * | v^(i+1)_0       1      | | y_{2i+1} |   | x^(i+1)_0   |
     *
     * Only interested in the value of y_2i however
     */
    
    for(int kz = 0; kz <= maxmode; kz++) {
      dcomplex v0, x0;
      
      // Get x and v0 from processor
      x0 = dcomplex(data.rcv[4*kz], data.rcv[4*kz+1]);
      v0 = dcomplex(data.rcv[4*kz+2], data.rcv[4*kz+3]);
      
      data.y2i[kz] = (data.xk[kz][mesh->xend] - data.w[kz][mesh->xend]*x0) / (1. - data.w[kz][mesh->xend]*v0);
      
    }
  }
  
  if(!mesh->firstX()) {
    // All except pe=0 receive values from i-1. Posting async receive
    data.recv_handle = mesh->irecvXIn(data.rcv, 2*(maxmode+1), PDD_COMM_Y);
  }
  
  if(!mesh->lastX()) {
    // Send value to the (i+1)th processor
    
    for(int kz = 0; kz <= maxmode; kz++) {
      data.snd[2*kz]   = data.y2i[kz].real();
      data.snd[2*kz+1] = data.y2i[kz].imag();
    }
    
    mesh->sendXOut(data.snd, 2*(maxmode+1), PDD_COMM_Y);
  }
}

/// Last part of the PDD algorithm
void LaplacePDD::finish(PDD_data &data, FieldPerp &x) {
  int ix, kz;

  x.allocate();
  x.setIndex(data.jy);
  
  if(!mesh->lastX()) {
    for(kz = 0; kz <= maxmode; kz++) {
      for(ix=0; ix < mesh->LocalNx; ix++)
	data.xk[kz][ix] -= data.w[kz][ix] * data.y2i[kz];
    }
  }

  if(!mesh->firstX()) {
    mesh->wait(data.recv_handle);
  
    for(kz = 0; kz <= maxmode; kz++) {
      dcomplex y2m = dcomplex(data.rcv[2*kz], data.rcv[2*kz+1]);
      
      for(ix=0; ix < mesh->LocalNx; ix++)
	data.xk[kz][ix] -= data.v[kz][ix] * y2m;
    }
  }
  
  // Have result in Fourier space. Convert back to BoutReal space

  static dcomplex *xk1d = NULL; ///< 1D in Z for taking FFTs

  int ncz = mesh->LocalNz;

  if(xk1d == NULL) {
    xk1d = new dcomplex[ncz/2 + 1];
    for(kz=0;kz<=ncz/2;kz++)
      xk1d[kz] = 0.0;
  }

  for(ix=0; ix<mesh->LocalNx; ix++){
    
    for(kz = 0; kz <= maxmode; kz++) {
      xk1d[kz] = data.xk[kz][ix];
    }

    if(global_flags & INVERT_ZERO_DC)
      xk1d[0] = 0.0;

    irfft(xk1d, ncz, x[ix]);
    
  }
}
