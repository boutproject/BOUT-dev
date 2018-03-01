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
#include <bout/openmpwrap.hxx>
#include <bout/sys/timer.hxx>
#include <boutexception.hxx>
#include <fft.hxx>
#include <bout/sys/timer.hxx>
#include <bout/scorepwrapper.hxx>
#include <bout/openmpwrap.hxx>
#include <globals.hxx>
#include <utils.hxx>

#include "spt.hxx"

LaplaceSPT::LaplaceSPT(Options *opt)
    : Laplacian(opt), Acoef(0.0), Ccoef(1.0), Dcoef(1.0) {

  if(mesh->periodicX) {
      throw BoutException("LaplaceSPT does not work with periodicity in the x direction (mesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
    }
	
  // Get start and end indices
  ys = mesh->ystart;
  ye = mesh->yend;
  if(mesh->hasBndryLowerY() && include_yguards)
    ys = 0; // Mesh contains a lower boundary
  if(mesh->hasBndryUpperY() && include_yguards)
    ye = mesh->LocalNy-1; // Contains upper boundary
  
  alldata = new SPT_data[ye - ys + 1];
  alldata -= ys; // Re-number indices to start at ys
  for(int jy=ys;jy<=ye;jy++) {
    alldata[jy].comm_tag = SPT_DATA + jy; // Give each one a different tag
  }

  // Temporary array for taking FFTs
  int ncz = mesh->LocalNz;
  dc1d = Array<dcomplex>(ncz / 2 + 1);
}

LaplaceSPT::~LaplaceSPT() {
  alldata += ys; // Return to index from 0
  delete[] alldata;
}

const FieldPerp LaplaceSPT::solve(const FieldPerp &b) {
  SCOREP0();
  return solve(b,b);
}

const FieldPerp LaplaceSPT::solve(const FieldPerp &b, const FieldPerp &x0) {
  SCOREP0();
  Mesh *mesh = b.getMesh();
  FieldPerp x(mesh);
  x.allocate();
  
  if( (inner_boundary_flags & INVERT_SET) || (outer_boundary_flags & INVERT_SET) ) {
    FieldPerp bs = copy(b);
    
    int xbndry = 2;
    if(global_flags & INVERT_BOTH_BNDRY_ONE)
      xbndry = 1;
    if((inner_boundary_flags & INVERT_SET) && mesh->firstX()) {
      // Copy x0 inner boundary into bs
      for(int ix=0;ix<xbndry;ix++)
        for(int iz=0;iz<mesh->LocalNz;iz++)
          bs[ix][iz] = x0[ix][iz];
    }
    if((outer_boundary_flags & INVERT_SET) && mesh->lastX()) {
      // Copy x0 outer boundary into bs
      for(int ix=mesh->LocalNx-1;ix>=mesh->LocalNx-xbndry;ix--)
        for(int iz=0;iz<mesh->LocalNz;iz++)
          bs[ix][iz] = x0[ix][iz];
    }
    start(bs, slicedata);
  }else
    start(b, slicedata);
  finish(slicedata, x);
  
  return x;
}

/// Extracts perpendicular slices from 3D fields and inverts separately
/*!
 * In parallel (mesh->NXPE > 1) this tries to overlap computation and communication.
 * This is done at the expense of more memory useage. Setting low_mem
 * in the config file uses less memory, and less communication overlap
 */
const Field3D LaplaceSPT::solve(const Field3D &b) {
  SCOREP0();
  Timer timer("invert");
  Mesh *mesh = b.getMesh();
  Field3D x(mesh);
  x.allocate();
  
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

  FieldPerp xperp(mesh);
  xperp.allocate();
  
  // All calculations finished. Get result
  for(int jy=ys; jy <= ye; jy++) {
    finish(alldata[jy], xperp);
    x = xperp;
  }
  
  x.setLocation(b.getLocation());
  
  return x;
}

const Field3D LaplaceSPT::solve(const Field3D &b, const Field3D &x0) {
  SCOREP0();
  if(  ((inner_boundary_flags & INVERT_SET) && mesh->firstX()) ||
       ((outer_boundary_flags & INVERT_SET) && mesh->lastX()) ) {
    Field3D bs = copy(b);
    
    int xbndry = 2;
    if(global_flags & INVERT_BOTH_BNDRY_ONE)
      xbndry = 1;
    
    if((inner_boundary_flags & INVERT_SET) && mesh->firstX()) {
      // Copy x0 inner boundary into bs
      for(int ix=0;ix<xbndry;ix++)
        for(int iy=0;iy<mesh->LocalNy;iy++)
          for(int iz=0;iz<mesh->LocalNz;iz++)
            bs(ix,iy,iz) = x0(ix,iy,iz);
    }
    if((outer_boundary_flags & INVERT_SET) && mesh->lastX()) {
      // Copy x0 outer boundary into bs
      for(int ix=mesh->LocalNx-1;ix>=mesh->LocalNx-xbndry;ix--)
        for(int iy=0;iy<mesh->LocalNy;iy++)
          for(int iz=0;iz<mesh->LocalNz;iz++)
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
  SCOREP0();
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
 * @param[inout] gp   gam from the processor mesh->PE_XIND + 1, and returned to mesh->PE_XIND - 1
 * @param[inout] up   u from processor mesh->PE_XIND + 1, and returned to mesh->PE_XIND - 1
 */
void LaplaceSPT::tridagBack(dcomplex *u, int n,
                             dcomplex *gam, dcomplex &gp, dcomplex &up) {
  SCOREP0();
  int j;

  u[n-1] = u[n-1] - gp*up;

  for(j=n-2;j>=0;j--) {
    u[j] = u[j]-gam[j+1]*u[j+1];
  }
  gp = gam[0];
  up = u[0];
}

/// Simple parallelisation of the Thomas tridiagonal solver algorithm (serial code)
/*!
 * This is a reference code which performs the same operations as the serial code.
 * To invert a single XZ slice (FieldPerp object), data must pass from the innermost
 * processor (mesh->PE_XIND = 0) to the outermost (mesh->PE_XIND = mesh->NXPE-1) and back again.
 *
 * Some parallelism is achieved by running several inversions simultaneously, so while
 * processor #1 is inverting Y=0, processor #0 is starting on Y=1. This works ok as long
 * as the number of slices to be inverted is greater than the number of X processors (MYSUB > mesh->NXPE).
 * If MYSUB < mesh->NXPE then not all processors can be busy at once, and so efficiency will fall sharply.
 *
 * @param[in]    b      RHS values (Ax = b)
 * @param[in]    global_flags  Inversion settings (see boundary.h for values)
 * @param[in]    inner_boundary_flags  Inversion settings for inner boundary (see invert_laplace.hxx for values)
 * @param[in]    outer_boundary_flags  Inversion settings for outer boundary (see invert_laplace.hxx for values)
 * @param[in]    a      This is a 2D matrix which allows solution of A = Delp2 + a
 * @param[out]   data   Structure containing data needed for second half of inversion
 * @param[in]    ccoef  Optional coefficient for first-order derivative
 * @param[in]    d      Optional factor to multiply the Delp2 operator
 */
int LaplaceSPT::start(const FieldPerp &b, SPT_data &data) {
  SCOREP0();
  if(mesh->firstX() && mesh->lastX())
    throw BoutException("Error: SPT method only works for mesh->NXPE > 1\n");

  data.jy = b.getIndex();

  int mm = mesh->LocalNz/2 + 1;
  data.allocate(mm, mesh->LocalNx); // Make sure data is allocated. Already allocated -> does nothing
  
  /// Take FFTs of data

  int ncz = mesh->LocalNz;
  
  for(int ix=0; ix < mesh->LocalNx; ix++) {
    rfft(b[ix], ncz, std::begin(dc1d));
    for(int kz = 0; kz <= maxmode; kz++)
      data.bk(kz, ix) = dc1d[kz];
  }

  BoutReal kwaveFactor = 2.0 * PI / mesh->coordinates()->zlength();

  /// Set matrix elements
  for (int kz = 0; kz <= maxmode; kz++) {
    tridagMatrix(&data.avec(kz, 0), &data.bvec(kz, 0), &data.cvec(kz, 0), &data.bk(kz, 0),
                 kz, kz * kwaveFactor, data.jy, global_flags, inner_boundary_flags,
                 outer_boundary_flags, &Acoef, &Ccoef, &Dcoef);
  }

  data.proc = 0; //< Starts at processor 0
  data.dir = 1;
  
  if(mesh->firstX()) {
    BOUT_OMP(parallel for)
    for(int kz = 0; kz <= maxmode; kz++) {
      dcomplex bet, u0;
      // Start tridiagonal solve
      tridagForward(&data.avec(kz, 0), &data.bvec(kz, 0), &data.cvec(kz, 0),
                    &data.bk(kz, 0), &data.xk(kz, 0), mesh->xend + 1, &data.gam(kz, 0),
                    bet, u0, true);
      // Load intermediate values into buffers
      data.buffer[4*kz]     = bet.real();
      data.buffer[4*kz + 1] = bet.imag();
      data.buffer[4*kz + 2] = u0.real();
      data.buffer[4*kz + 3] = u0.imag();
    }
    
    // Send data
    mesh->sendXOut(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);

  }else if(mesh->PE_XIND == 1) {
    // Post a receive
    data.recv_handle =
        mesh->irecvXIn(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
  }
  
  data.proc++; // Now moved onto the next processor
  if(mesh->NXPE == 2)	
    data.dir = -1; // Special case. Otherwise reversal handled in spt_continue
  
  return 0;
}

/// Shifts the parallelised Thomas algorithm along one processor.
/*!
  Returns non-zero when the calculation is complete.

  @param[inout] data  Structure which keeps track of the calculation
*/
int LaplaceSPT::next(SPT_data &data) {
  SCOREP0();
  if(data.proc < 0) // Already finished
    return 1;
  
  if(mesh->PE_XIND == data.proc) {
    /// This processor's turn to do inversion

    // Wait for data to arrive
    mesh->wait(data.recv_handle);

    if(mesh->lastX()) {
      // Last processor, turn-around
      
      BOUT_OMP(parallel for)
      for(int kz = 0; kz <= maxmode; kz++) {
        dcomplex bet, u0;
        dcomplex gp, up;
	bet = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	u0 = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);
        tridagForward(&data.avec(kz, mesh->xstart), &data.bvec(kz, mesh->xstart),
                      &data.cvec(kz, mesh->xstart), &data.bk(kz, mesh->xstart),
                      &data.xk(kz, mesh->xstart), mesh->xend + 1,
                      &data.gam(kz, mesh->xstart), bet, u0);

        // Back-substitute
	gp = 0.0;
	up = 0.0;
        tridagBack(&data.xk(kz, mesh->xstart), mesh->LocalNx - mesh->xstart,
                   &data.gam(kz, mesh->xstart), gp, up);
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
        tridagForward(&data.avec(kz, mesh->xstart), &data.bvec(kz, mesh->xstart),
                      &data.cvec(kz, mesh->xstart), &data.bk(kz, mesh->xstart),
                      &data.xk(kz, mesh->xstart), mesh->xend - mesh->xstart + 1,
                      &data.gam(kz, mesh->xstart), bet, u0);
        // Load intermediate values into buffers
	data.buffer[4*kz]     = bet.real();
	data.buffer[4*kz + 1] = bet.imag();
	data.buffer[4*kz + 2] = u0.real();
	data.buffer[4*kz + 3] = u0.imag();
      }
      
    }else if(mesh->firstX()) {
      // Back to the start

BOUT_OMP(parallel for)
      for(int kz = 0; kz <= maxmode; kz++) {
	dcomplex gp, up;
	gp = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	up = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);

        tridagBack(&data.xk(kz, 0), mesh->xend + 1, &data.gam(kz, 0), gp, up);
      }

    }else {
      // Middle of X, back-substitution stage

      BOUT_OMP(parallel for)
      for(int kz = 0; kz <= maxmode; kz++) {
	dcomplex gp = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	dcomplex up = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);

        tridagBack(&data.xk(kz, mesh->xstart), mesh->xend - mesh->xstart + 1,
                   &data.gam(kz, mesh->xstart), gp, up);

        data.buffer[4*kz]     = gp.real();
	data.buffer[4*kz + 1] = gp.imag();
	data.buffer[4*kz + 2] = up.real();
	data.buffer[4*kz + 3] = up.imag();
      }
    }

    if(mesh->PE_XIND != 0) { // If not finished yet
      /// Send data
      
      if(data.dir > 0) {
        mesh->sendXOut(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
      }else
        mesh->sendXIn(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
    }

  }else if(mesh->PE_XIND == data.proc + data.dir) {
    // This processor is next, post receive
    
    if(data.dir > 0) {
      data.recv_handle =
          mesh->irecvXIn(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
    }else
      data.recv_handle =
          mesh->irecvXOut(std::begin(data.buffer), 4 * (maxmode + 1), data.comm_tag);
  }
  
  data.proc += data.dir;
  
  if(data.proc == mesh->NXPE-1)
    data.dir = -1; // Reverses direction at the end

  return 0;
}

/// Finishes the parallelised Thomas algorithm
/*!
  @param[inout] data   Structure keeping track of calculation
  @param[in]    global_flags  Inversion flags (same as passed to invert_spt_start)
  @param[in]    inner_boundary_flags  Inversion flags for inner boundary (same as passed to invert_spt_start)
  @param[in]    outer_boundary_flags  Inversion flags for outer boundary (same as passed to invert_spt_start)
  @param[out]   x      The result
*/
void LaplaceSPT::finish(SPT_data &data, FieldPerp &x) {
  SCOREP0();
  int ncx = mesh->LocalNx-1;
  int ncz = mesh->LocalNz;

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

  if(!mesh->firstX()) {
    // Set left boundary to zero (Prevent unassigned values in corners)
    for(int ix=0; ix<mesh->xstart; ix++){
      for(int kz=0;kz<mesh->LocalNz;kz++)
	x(ix,kz) = 0.0;
    }
  }
  if(!mesh->lastX()) {
    // Same for right boundary
    for(int ix=mesh->xend+1; ix<mesh->LocalNx; ix++){
      for(int kz=0;kz<mesh->LocalNz;kz++)
	x(ix,kz) = 0.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// SPT_data helper class

void LaplaceSPT::SPT_data::allocate(int mm, int nx) {
  SCOREP0();
  bk = Matrix<dcomplex>(mm, nx);
  xk = Matrix<dcomplex>(mm, nx);

  gam = Matrix<dcomplex>(mm, nx);

  // Matrix to be solved
  avec = Matrix<dcomplex>(mm, nx);
  bvec = Matrix<dcomplex>(mm, nx);
  cvec = Matrix<dcomplex>(mm, nx);

  buffer = Array<BoutReal>(4 * mm);
}

