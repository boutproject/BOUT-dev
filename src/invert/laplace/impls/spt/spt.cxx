/*!
 * \file spt.cxx
 *
 * \brief Simple Parallel Tridiagonal solver
 *
 * Changelog
 * ---------
 * 
 * 2014-06  Ben Dudson <benjamin.dudson@york.ac.uk>
 *     * Removed static variables in functions, changing to class members.
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
 */

#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout/openmpwrap.hxx>
#include <bout/sys/timer.hxx>
#include <boutexception.hxx>
#include <fft.hxx>
#include <globals.hxx>
#include <utils.hxx>

#include "spt.hxx"

LaplaceSPT::LaplaceSPT(Options *opt, const CELL_LOC loc, Mesh *mesh_in)
    : Laplacian(opt, loc, mesh_in), Acoef(0.0), Ccoef(1.0), Dcoef(1.0) {
  Acoef.setLocation(location);
  Ccoef.setLocation(location);
  Dcoef.setLocation(location);

  if(localmesh->periodicX) {
      throw BoutException("LaplaceSPT does not work with periodicity in the x direction (localmesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
    }
	
  // Get start and end indices
  ys = localmesh->ystart;
  ye = localmesh->yend;
  if(localmesh->hasBndryLowerY() && include_yguards)
    ys = 0; // Mesh contains a lower boundary
  if(localmesh->hasBndryUpperY() && include_yguards)
    ye = localmesh->LocalNy-1; // Contains upper boundary
  
  alldata = new SPT_data[ye - ys + 1];
  alldata -= ys; // Re-number indices to start at ys
  for(int jy=ys;jy<=ye;jy++) {
    alldata[jy].comm_tag = SPT_DATA + jy; // Give each one a different tag
  }

  // Temporary array for taking FFTs
  int ncz = localmesh->LocalNz;
  dc1d.reallocate(ncz / 2 + 1);
}

LaplaceSPT::~LaplaceSPT() {
  alldata += ys; // Return to index from 0
  delete[] alldata;
}

FieldPerp LaplaceSPT::solve(const FieldPerp& b) { return solve(b, b); }

FieldPerp LaplaceSPT::solve(const FieldPerp& b, const FieldPerp& x0) {
  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  FieldPerp x{emptyFrom(b)};
  
  if( (inner_boundary_flags & INVERT_SET) || (outer_boundary_flags & INVERT_SET) ) {
    FieldPerp bs = copy(b);
    
    int xbndry = localmesh->xstart;
    // If the flags to assign that only one guard cell should be used is set
    if((global_flags & INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2))
      xbndry = 1;
    if((inner_boundary_flags & INVERT_SET) && localmesh->firstX()) {
      // Copy x0 inner boundary into bs
      for(int ix=0;ix<xbndry;ix++)
        for(int iz=0;iz<localmesh->LocalNz;iz++)
          bs[ix][iz] = x0[ix][iz];
    }
    if((outer_boundary_flags & INVERT_SET) && localmesh->lastX()) {
      // Copy x0 outer boundary into bs
      for(int ix=localmesh->LocalNx-1;ix>=localmesh->LocalNx-xbndry;ix--)
        for(int iz=0;iz<localmesh->LocalNz;iz++)
          bs[ix][iz] = x0[ix][iz];
    }
    start(bs, slicedata);
  }else
    start(b, slicedata);
  finish(slicedata, x);

  checkData(x);
  
  return x;
}

/// Extracts perpendicular slices from 3D fields and inverts separately
/*!
 * In parallel (localmesh->NXPE > 1) this tries to overlap computation and communication.
 * This is done at the expense of more memory useage. Setting low_mem
 * in the config file uses less memory, and less communication overlap
 */
Field3D LaplaceSPT::solve(const Field3D& b) {

  ASSERT1(b.getLocation() == location);
  ASSERT1(localmesh == b.getMesh());

  Timer timer("invert");
  Field3D x{emptyFrom(b)};
  
  for(int jy=ys; jy <= ye; jy++) {
    // And start another one going
    start(sliceXZ(b, jy), alldata[jy]);
    
    // Move each calculation along one processor
    for(int jy2=ys; jy2 < jy; jy2++) 
      next(alldata[jy2]);
  }
  
  bool running = true;
  do {
    // Move each calculation along until the last one is finished
    for(int jy=ys; jy <= ye; jy++)
      running = next(alldata[jy]) == 0;
  }while(running);

  FieldPerp xperp(localmesh);
  xperp.setLocation(location);
  xperp.allocate();
  
  // All calculations finished. Get result
  for(int jy=ys; jy <= ye; jy++) {
    finish(alldata[jy], xperp);
    x = xperp;
  }
  
  return x;
}

Field3D LaplaceSPT::solve(const Field3D& b, const Field3D& x0) {
  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());

  if(  ((inner_boundary_flags & INVERT_SET) && localmesh->firstX()) ||
       ((outer_boundary_flags & INVERT_SET) && localmesh->lastX()) ) {
    Field3D bs = copy(b);
    
    int xbndry = localmesh->xstart;
    // If the flags to assign that only one guard cell should be used is set
    if((global_flags & INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2))
      xbndry = 1;
    
    if((inner_boundary_flags & INVERT_SET) && localmesh->firstX()) {
      // Copy x0 inner boundary into bs
      for(int ix=0;ix<xbndry;ix++)
        for(int iy=0;iy<localmesh->LocalNy;iy++)
          for(int iz=0;iz<localmesh->LocalNz;iz++)
            bs(ix,iy,iz) = x0(ix,iy,iz);
    }
    if((outer_boundary_flags & INVERT_SET) && localmesh->lastX()) {
      // Copy x0 outer boundary into bs
      for(int ix=localmesh->LocalNx-1;ix>=localmesh->LocalNx-xbndry;ix--)
        for(int iy=0;iy<localmesh->LocalNy;iy++)
          for(int iz=0;iz<localmesh->LocalNz;iz++)
            bs(ix,iy,iz) = x0(ix,iy,iz);
    }
    return solve(bs);
  }
  
  return solve(b);
}

/// This is the first half of the Thomas algorithm for parallel calculations
/*!
 * Two complex quantities have to be propagated between processors: bet and u[-1].
 * This routine takes bet and um from the last processor (if start == false),
 * and returns the values to be passed to the next processor in the same variables.
 *
 * @param[in]  a    Vector of matrix coefficients (Left of diagonal)
 * @param[in]  b    Vector of matrix coefficients (Diagonal)
 * @param[in]  c    Vector of matrix coefficients (Right of diagonal)
 * @param[in]  r    RHS vector
 * @param[in]  u    Result vector (Au = r)
 * @param[in]  n    Size of the matrix
 * @param[out] gam  Intermediate values used for backsolve stage
 * @param[inout] bet
 * @param[inout] um
 * @param[in] start
 */
void LaplaceSPT::tridagForward(dcomplex *a, dcomplex *b, dcomplex *c,
                                dcomplex *r, dcomplex *u, int n,
                                dcomplex *gam,
                                dcomplex &bet, dcomplex &um, bool start) {
  int j;
  
  if(start) {
    bet = b[0];
    u[0] = r[0] / bet;
  }else {
    gam[0] = c[-1] / bet; // NOTE: ASSUMES C NOT CHANGING
    bet = b[0] - a[0]*gam[0];
    u[0] = (r[0]-a[0]*um)/bet;
  }
  
  for(j=1;j<n;j++) {
    gam[j] = c[j-1]/bet;
    bet = b[j]-a[j]*gam[j];
    if(bet == 0.0)
      throw BoutException("Tridag: Zero pivot\n");
    
    u[j] = (r[j]-a[j]*u[j-1])/bet;
  }

  um = u[n-1];
}

/// Second (backsolve) part of the Thomas algorithm
/*!
 * @param[inout] u    Result to be solved (Au = r)
 * @param[in]    n    Size of the problem
 * @param[in]    gam  Intermediate values produced by the forward part
 * @param[inout] gp   gam from the processor localmesh->PE_XIND + 1, and returned to localmesh->PE_XIND - 1
 * @param[inout] up   u from processor localmesh->PE_XIND + 1, and returned to localmesh->PE_XIND - 1
 */
void LaplaceSPT::tridagBack(dcomplex *u, int n,
                             dcomplex *gam, dcomplex &gp, dcomplex &up) {
  int j;

  u[n-1] = u[n-1] - gp*up;

  for(j=n-2;j>=0;j--) {
    u[j] = u[j]-gam[j+1]*u[j+1];
  }
  gp = gam[0];
  up = u[0];
}

/// Simple parallelisation of the Thomas tridiagonal solver algorithm
/// (serial code)
///
/// This is a reference code which performs the same operations as the
/// serial code.  To invert a single XZ slice (FieldPerp object), data
/// must pass from the innermost processor (localmesh->PE_XIND = 0) to the
/// outermost (localmesh->PE_XIND = localmesh->NXPE-1) and back again.
///
/// Some parallelism is achieved by running several inversions
/// simultaneously, so while processor #1 is inverting Y=0, processor
/// #0 is starting on Y=1. This works ok as long as the number of
/// slices to be inverted is greater than the number of X processors
/// (MYSUB > localmesh->NXPE).  If MYSUB < localmesh->NXPE then not all
/// processors can be busy at once, and so efficiency will fall
/// sharply.
///
/// @param[in]    b      RHS values (Ax = b)
/// @param[out]   data   Structure containing data needed for second half of inversion
int LaplaceSPT::start(const FieldPerp &b, SPT_data &data) {
  if(localmesh->firstX() && localmesh->lastX())
    throw BoutException("Error: SPT method only works for localmesh->NXPE > 1\n");

  ASSERT1(b.getLocation() == location);

  data.jy = b.getIndex();

  int mm = localmesh->LocalNz/2 + 1;
  data.allocate(mm, localmesh->LocalNx); // Make sure data is allocated. Already allocated -> does nothing
  
  /// Take FFTs of data

  int ncz = localmesh->LocalNz;
  
  for(int ix=0; ix < localmesh->LocalNx; ix++) {
    rfft(b[ix], ncz, std::begin(dc1d));
    for(int kz = 0; kz <= maxmode; kz++)
      data.bk(kz, ix) = dc1d[kz];
  }

  BoutReal kwaveFactor = 2.0 * PI / coords->zlength();

  /// Set matrix elements
  for (int kz = 0; kz <= maxmode; kz++) {
    tridagMatrix(&data.avec(kz, 0), &data.bvec(kz, 0), &data.cvec(kz, 0), &data.bk(kz, 0),
                 data.jy, kz, kz * kwaveFactor, global_flags, inner_boundary_flags,
                 outer_boundary_flags, &Acoef, &Ccoef, &Dcoef);
  }

  data.proc = 0; //< Starts at processor 0
  data.dir = 1;
  
  if(localmesh->firstX()) {
    BOUT_OMP(parallel for)
    for(int kz = 0; kz <= maxmode; kz++) {
      dcomplex bet, u0;
      // Start tridiagonal solve
      tridagForward(&data.avec(kz, 0), &data.bvec(kz, 0), &data.cvec(kz, 0),
                    &data.bk(kz, 0), &data.xk(kz, 0), localmesh->xend + 1, &data.gam(kz, 0),
                    bet, u0, true);
      // Load intermediate values into buffers
      data.buffer[4*kz]     = bet.real();
      data.buffer[4*kz + 1] = bet.imag();
      data.buffer[4*kz + 2] = u0.real();
      data.buffer[4*kz + 3] = u0.imag();
    }
    
    // Send data
    localmesh->sendXOut(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);

  }else if(localmesh->PE_XIND == 1) {
    // Post a receive
    data.recv_handle =
        localmesh->irecvXIn(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
  }
  
  data.proc++; // Now moved onto the next processor
  if(localmesh->NXPE == 2)	
    data.dir = -1; // Special case. Otherwise reversal handled in spt_continue
  
  return 0;
}

/// Shifts the parallelised Thomas algorithm along one processor.
/*!
  Returns non-zero when the calculation is complete.

  @param[inout] data  Structure which keeps track of the calculation
*/
int LaplaceSPT::next(SPT_data &data) {
  if(data.proc < 0) // Already finished
    return 1;
  
  if(localmesh->PE_XIND == data.proc) {
    /// This processor's turn to do inversion

    // Wait for data to arrive
    localmesh->wait(data.recv_handle);

    if(localmesh->lastX()) {
      // Last processor, turn-around
      
      BOUT_OMP(parallel for)
      for(int kz = 0; kz <= maxmode; kz++) {
        dcomplex bet, u0;
        dcomplex gp, up;
	bet = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	u0 = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);
        tridagForward(&data.avec(kz, localmesh->xstart), &data.bvec(kz, localmesh->xstart),
                      &data.cvec(kz, localmesh->xstart), &data.bk(kz, localmesh->xstart),
                      &data.xk(kz, localmesh->xstart), localmesh->xend + 1,
                      &data.gam(kz, localmesh->xstart), bet, u0);

        // Back-substitute
	gp = 0.0;
	up = 0.0;
        tridagBack(&data.xk(kz, localmesh->xstart), localmesh->LocalNx - localmesh->xstart,
                   &data.gam(kz, localmesh->xstart), gp, up);
        data.buffer[4*kz]     = gp.real();
	data.buffer[4*kz + 1] = gp.imag();
	data.buffer[4*kz + 2] = up.real();
	data.buffer[4*kz + 3] = up.imag();
      }

    }else if(data.dir > 0) {
      // In the middle of X, forward direction

      BOUT_OMP(parallel for)
      for(int kz = 0; kz <= maxmode; kz++) {
	dcomplex bet, u0;
	bet = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	u0 = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);
        tridagForward(&data.avec(kz, localmesh->xstart), &data.bvec(kz, localmesh->xstart),
                      &data.cvec(kz, localmesh->xstart), &data.bk(kz, localmesh->xstart),
                      &data.xk(kz, localmesh->xstart), localmesh->xend - localmesh->xstart + 1,
                      &data.gam(kz, localmesh->xstart), bet, u0);
        // Load intermediate values into buffers
	data.buffer[4*kz]     = bet.real();
	data.buffer[4*kz + 1] = bet.imag();
	data.buffer[4*kz + 2] = u0.real();
	data.buffer[4*kz + 3] = u0.imag();
      }
      
    }else if(localmesh->firstX()) {
      // Back to the start

BOUT_OMP(parallel for)
      for(int kz = 0; kz <= maxmode; kz++) {
	dcomplex gp, up;
	gp = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	up = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);

        tridagBack(&data.xk(kz, 0), localmesh->xend + 1, &data.gam(kz, 0), gp, up);
      }

    }else {
      // Middle of X, back-substitution stage

      BOUT_OMP(parallel for)
      for(int kz = 0; kz <= maxmode; kz++) {
	dcomplex gp = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	dcomplex up = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);

        tridagBack(&data.xk(kz, localmesh->xstart), localmesh->xend - localmesh->xstart + 1,
                   &data.gam(kz, localmesh->xstart), gp, up);

        data.buffer[4*kz]     = gp.real();
	data.buffer[4*kz + 1] = gp.imag();
	data.buffer[4*kz + 2] = up.real();
	data.buffer[4*kz + 3] = up.imag();
      }
    }

    if(localmesh->PE_XIND != 0) { // If not finished yet
      /// Send data
      
      if(data.dir > 0) {
        localmesh->sendXOut(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
      }else
        localmesh->sendXIn(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
    }

  }else if(localmesh->PE_XIND == data.proc + data.dir) {
    // This processor is next, post receive
    
    if(data.dir > 0) {
      data.recv_handle =
          localmesh->irecvXIn(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
    }else
      data.recv_handle =
          localmesh->irecvXOut(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
  }
  
  data.proc += data.dir;
  
  if(data.proc == localmesh->NXPE-1)
    data.dir = -1; // Reverses direction at the end

  return 0;
}

/// Finishes the parallelised Thomas algorithm
///
/// @param[inout] data   Structure keeping track of calculation
/// @param[out]   x      The result
void LaplaceSPT::finish(SPT_data &data, FieldPerp &x) {
  int ncx = localmesh->LocalNx-1;
  int ncz = localmesh->LocalNz;

  ASSERT1(x.getLocation() == location);

  x.allocate();
  x.setIndex(data.jy);

  // Make sure calculation has finished
  while(next(data) == 0) {}

  // Have result in Fourier space. Convert back to real space
  
  for(int ix=0; ix<=ncx; ix++){
    
    for(int kz = 0; kz<= maxmode; kz++) {
      dc1d[kz] = data.xk(kz, ix);
    }
    for(int kz = maxmode + 1; kz <= ncz/2; kz++)
      dc1d[kz] = 0.0;

    if(global_flags & INVERT_ZERO_DC)
      dc1d[0] = 0.0;

    irfft(std::begin(dc1d), ncz, x[ix]);
  }

  if(!localmesh->firstX()) {
    // Set left boundary to zero (Prevent unassigned values in corners)
    for(int ix=0; ix<localmesh->xstart; ix++){
      for(int kz=0;kz<localmesh->LocalNz;kz++)
	x(ix,kz) = 0.0;
    }
  }
  if(!localmesh->lastX()) {
    // Same for right boundary
    for(int ix=localmesh->xend+1; ix<localmesh->LocalNx; ix++){
      for(int kz=0;kz<localmesh->LocalNz;kz++)
	x(ix,kz) = 0.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// SPT_data helper class

void LaplaceSPT::SPT_data::allocate(int mm, int nx) {
  bk.reallocate(mm, nx);
  xk.reallocate(mm, nx);

  gam.reallocate(mm, nx);

  // Matrix to be solved
  avec.reallocate(mm, nx);
  bvec.reallocate(mm, nx);
  cvec.reallocate(mm, nx);

  buffer.reallocate(4 * mm);
}

