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

#include <bout/constants.hxx>
#include <boutexception.hxx>
#include <fft.hxx>
#include <globals.hxx>
#include <lapack_routines.hxx>
#include <utils.hxx>

#include "pdd.hxx"

FieldPerp LaplacePDD::solve(const FieldPerp& b) {
  ASSERT1(localmesh == b.getMesh());
  ASSERT1(b.getLocation() == location);

  PDD_data data;

  FieldPerp x{emptyFrom(b)};
  
  start(b, data);
  next(data);
  finish(data, x);
  
  return x;
}

Field3D LaplacePDD::solve(const Field3D& b) {
  ASSERT1(localmesh == b.getMesh());
  ASSERT1(b.getLocation() == location);

  Field3D x{emptyFrom(b)};
  FieldPerp xperp(localmesh);
  xperp.allocate();
  
  int ys = localmesh->ystart, ye = localmesh->yend;
  if(localmesh->hasBndryLowerY())
    ys = 0; // Mesh contains a lower boundary
  if(localmesh->hasBndryUpperY())
    ye = localmesh->LocalNy-1; // Contains upper boundary
  
  if(low_mem) {
    // Solve one slice at a time
    for(int jy=ys; jy <= ye; jy++) {
      x = solve(sliceXZ(b, jy));
    }
  }else {
    // Overlap multiple inversions

    static PDD_data *data = nullptr;

    if (data == nullptr) {
      data = new PDD_data[ye - ys + 1];
      data -= ys; // Re-number indices to start at jstart
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
///
/// July 2008: Adapted from serial version to run in parallel (split
/// in X) for tridiagonal system i.e. no 4th order inversion yet.
///
/// \note This code stores intermediate results and takes
/// significantly more memory than the serial version. This can be
/// balanced against communication time i.e. faster communications can
/// allow less memory use.
///
/// @param[in]    b  RHS values (Ax = b)
/// @param[in] data  Internal data used for multiple calls in parallel mode
void LaplacePDD::start(const FieldPerp &b, PDD_data &data) {
  ASSERT1(localmesh == b.getMesh());
  ASSERT1(b.getLocation() == location);

  int ix, kz;

  int ncz = localmesh->LocalNz;

  data.jy = b.getIndex();

  if(localmesh->firstX() && localmesh->lastX())
    throw BoutException("Error: PDD method only works for NXPE > 1\n");
  
  if(localmesh->periodicX) {
      throw BoutException("LaplacePDD does not work with periodicity in the x direction (localmesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
    }

    if (data.bk.empty()) {
      // Need to allocate working memory

      // RHS vector
      data.bk.reallocate(maxmode + 1, localmesh->LocalNx);

      // Matrix to be solved
      data.avec.reallocate(maxmode + 1, localmesh->LocalNx);
      data.bvec.reallocate(maxmode + 1, localmesh->LocalNx);
      data.cvec.reallocate(maxmode + 1, localmesh->LocalNx);

      // Working vectors
      data.v.reallocate(maxmode + 1, localmesh->LocalNx);
      data.w.reallocate(maxmode + 1, localmesh->LocalNx);

      // Result
      data.xk.reallocate(maxmode + 1, localmesh->LocalNx);

      // Communication buffers. Space for 2 complex values for each kz
      data.snd.reallocate(4 * (maxmode + 1));
      data.rcv.reallocate(4 * (maxmode + 1));

      data.y2i.reallocate(maxmode + 1);
  }

  /// Take FFTs of data
  Array<dcomplex> bk1d(ncz / 2 + 1); ///< 1D in Z for taking FFTs

  for(ix=0; ix < localmesh->LocalNx; ix++) {
    rfft(b[ix], ncz, std::begin(bk1d));
    for(kz = 0; kz <= maxmode; kz++)
      data.bk(kz, ix) = bk1d[kz];
  }

  /// Create the matrices to be inverted (one for each z point)

  ASSERT1(isConst(coords->zlength()));
  BoutReal kwaveFactor = 2.0 * PI / coords->zlength()(0,0);

  /// Set matrix elements
  for (int kz = 0; kz <= maxmode; kz++) {
    tridagMatrix(&data.avec(kz, 0), &data.bvec(kz, 0), &data.cvec(kz, 0), &data.bk(kz, 0),
                 kz, kz * kwaveFactor, data.jy, global_flags, inner_boundary_flags,
                 outer_boundary_flags, &Acoef, &Ccoef, &Dcoef);
  }

  Array<dcomplex> e(localmesh->LocalNx);
  for (ix = 0; ix < localmesh->LocalNx; ix++)
    e[ix] = 0.0; // Do we need this?

  for(kz = 0; kz <= maxmode; kz++) {
    // Start PDD algorithm

    // Solve for xtilde, v and w (step 2)

    dcomplex v0, x0; // Values to be sent to processor i-1

    if(localmesh->firstX()) {
      // Domain includes inner boundary
      tridag(&data.avec(kz, 0), &data.bvec(kz, 0), &data.cvec(kz, 0), &data.bk(kz, 0),
             &data.xk(kz, 0), localmesh->xend + 1);

      // Add C (row m-1) from next processor

      e[localmesh->xend] = data.cvec(kz, localmesh->xend);
      tridag(&data.avec(kz, 0), &data.bvec(kz, 0), &data.cvec(kz, 0), std::begin(e),
             &data.w(kz, 0), localmesh->xend + 1);

    }else if(localmesh->lastX()) {
      // Domain includes outer boundary
      tridag(&data.avec(kz, localmesh->xstart), &data.bvec(kz, localmesh->xstart),
             &data.cvec(kz, localmesh->xstart), &data.bk(kz, localmesh->xstart),
             &data.xk(kz, localmesh->xstart), localmesh->xend - localmesh->xend + 1);

      // Add A (row 0) from previous processor
      e[0] = data.avec(kz, localmesh->xstart);
      tridag(&data.avec(kz, localmesh->xstart), &data.bvec(kz, localmesh->xstart),
             &data.cvec(kz, localmesh->xstart), std::begin(e), &data.v(kz, localmesh->xstart),
             localmesh->xend + 1);

      x0 = data.xk(kz, localmesh->xstart);
      v0 = data.v(kz, localmesh->xstart);

    }else {
      // No boundaries
      tridag(&data.avec(kz, localmesh->xstart), &data.bvec(kz, localmesh->xstart),
             &data.cvec(kz, localmesh->xstart), &data.bk(kz, localmesh->xstart),
             &data.xk(kz, localmesh->xstart), localmesh->xend - localmesh->xstart + 1);

      // Add A (row 0) from previous processor
      e[0] = data.avec(kz, localmesh->xstart);
      tridag(&data.avec(kz, localmesh->xstart), &data.bvec(kz, localmesh->xstart),
             &data.cvec(kz, localmesh->xstart), &e[localmesh->xstart], &data.v(kz, localmesh->xstart),
             localmesh->xend - localmesh->xstart + 1);
      e[0] = 0.0;
      
      // Add C (row m-1) from next processor
      e[localmesh->xend] = data.cvec(kz, localmesh->xend);
      tridag(&data.avec(kz, localmesh->xstart), &data.bvec(kz, localmesh->xstart),
             &data.cvec(kz, localmesh->xstart), &e[localmesh->xstart], &data.v(kz, localmesh->xstart),
             localmesh->xend - localmesh->xstart + 1);
      e[localmesh->xend] = 0.0;
    }
    
    // Put values into communication buffers
    data.snd[4*kz]   = x0.real();
    data.snd[4*kz+1] = x0.imag();
    data.snd[4*kz+2] = v0.real();
    data.snd[4*kz+3] = v0.imag();
  }
  
  // Stage 3: Communicate x0, v0 from node i to i-1
  
  if(!localmesh->lastX()) {
    // All except the last processor expect to receive data
    // Post async receive
    data.recv_handle =
        localmesh->irecvXOut(std::begin(data.rcv), 4 * (maxmode + 1), PDD_COMM_XV);
  }

  if(!localmesh->firstX()) {
    // Send the data

    localmesh->sendXIn(std::begin(data.snd), 4 * (maxmode + 1), PDD_COMM_XV);
  }
}


/// Middle part of the PDD algorithm
void LaplacePDD::next(PDD_data &data) {
  // Wait for x0 and v0 to arrive from processor i+1
  
  if(!localmesh->lastX()) {
    localmesh->wait(data.recv_handle);

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

      data.y2i[kz] = (data.xk(kz, localmesh->xend) - data.w(kz, localmesh->xend) * x0) /
                     (1. - data.w(kz, localmesh->xend) * v0);
    }
  }
  
  if(!localmesh->firstX()) {
    // All except pe=0 receive values from i-1. Posting async receive
    data.recv_handle =
        localmesh->irecvXIn(std::begin(data.rcv), 2 * (maxmode + 1), PDD_COMM_Y);
  }
  
  if(!localmesh->lastX()) {
    // Send value to the (i+1)th processor
    
    for(int kz = 0; kz <= maxmode; kz++) {
      data.snd[2*kz]   = data.y2i[kz].real();
      data.snd[2*kz+1] = data.y2i[kz].imag();
    }

    localmesh->sendXOut(std::begin(data.snd), 2 * (maxmode + 1), PDD_COMM_Y);
  }
}

/// Last part of the PDD algorithm
void LaplacePDD::finish(PDD_data &data, FieldPerp &x) {
  ASSERT1(x.getLocation() == location);

  int ix, kz;

  x.allocate();
  x.setIndex(data.jy);
  
  if(!localmesh->lastX()) {
    for(kz = 0; kz <= maxmode; kz++) {
      for(ix=0; ix < localmesh->LocalNx; ix++)
        data.xk(kz, ix) -= data.w(kz, ix) * data.y2i[kz];
    }
  }

  if(!localmesh->firstX()) {
    localmesh->wait(data.recv_handle);
  
    for(kz = 0; kz <= maxmode; kz++) {
      dcomplex y2m = dcomplex(data.rcv[2*kz], data.rcv[2*kz+1]);
      
      for(ix=0; ix < localmesh->LocalNx; ix++)
        data.xk(kz, ix) -= data.v(kz, ix) * y2m;
    }
  }
  
  // Have result in Fourier space. Convert back to BoutReal space
  int ncz = localmesh->LocalNz;

  Array<dcomplex> xk1d(ncz / 2 + 1); ///< 1D in Z for taking FFTs
  for (kz = maxmode; kz <= ncz / 2; kz++)
    xk1d[kz] = 0.0;

  for(ix=0; ix<localmesh->LocalNx; ix++){
    
    for(kz = 0; kz <= maxmode; kz++) {
      xk1d[kz] = data.xk(kz, ix);
    }

    if(global_flags & INVERT_ZERO_DC)
      xk1d[0] = 0.0;

    irfft(std::begin(xk1d), ncz, x[ix]);
  }
}
