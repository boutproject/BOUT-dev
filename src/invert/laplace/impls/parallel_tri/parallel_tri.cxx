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
#include "parallel_tri.hxx"

#include <bout/mesh.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <lapack_routines.hxx>
#include <bout/constants.hxx>
#include <bout/openmpwrap.hxx>
#include <cmath>
#include <bout/sys/timer.hxx>

#include <output.hxx>
#include "boutcomm.hxx"

#include <bout/scorepwrapper.hxx>

LaplaceParallelTri::LaplaceParallelTri(Options *opt, CELL_LOC loc, Mesh *mesh_in)
    : Laplacian(opt, loc, mesh_in), A(0.0), C(1.0), D(1.0), ipt_mean_its(0.), ncalls(0), Borig(50.), Bvals(1000.) {
  A.setLocation(location);
  C.setLocation(location);
  D.setLocation(location);

  OPTION(opt, rtol, 1.e-7);
  OPTION(opt, atol, 1.e-20);
  OPTION(opt, maxits, 100);
  OPTION(opt, B, 1000.0);
  OPTION(opt, om, 1.0);

  static int ipt_solver_count = 1;
  bout::globals::dump.addRepeat(ipt_mean_its,
      "ipt_solver"+std::to_string(ipt_solver_count)+"_mean_its");
  ++ipt_solver_count;

  Borig = B;

  first_call = Matrix<bool>(localmesh->LocalNy,localmesh->LocalNz);
  for(int jy=0; jy<localmesh->LocalNy; jy++){
    for(int jz=0; jz<localmesh->LocalNz; jz++){
      first_call(jy,jz) = true;
    }
  }

  x0saved = Tensor<dcomplex>(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);

  resetSolver();

}

void LaplaceParallelTri::resetSolver(){
  x0saved = 0.0;
  resetMeanIterations();
}

/*!
 * Calcalate stability of the iteration, and amend right-hand side vector minvb to ensure stability.
 */
void LaplaceParallelTri::ensure_stability(const Array<dcomplex> &avec, const Array<dcomplex> &bvec,
                                              const Array<dcomplex> &cvec, const Array<dcomplex> &minvb,
				              const int ncx, Array<dcomplex> &xk1d) {
  SCOREP0();

  Array<dcomplex> xvec, evec;
  xvec = Array<dcomplex>(ncx);
  evec = Array<dcomplex>(ncx);

  Array<dcomplex> sendvec, recvec;
  sendvec = Array<dcomplex>(2);
  recvec = Array<dcomplex>(2);

  // If not on innermost boundary, get information from neighbouring proc and
  // calculate value of solution in halo cell
  if(!localmesh->firstX()) {

    comm_handle recv[1];
    recv[0] = localmesh->irecvXIn(&recvec[0], 2, 0);

    // Need the (xstart-1)-th element
    evec = Array<dcomplex>(ncx);
    for(int i=0; i<ncx; i++){
      evec[i] = 0.0;
    }
    evec[localmesh->xstart-1] = 1.0;

    tridag(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(evec),
             std::begin(xvec), ncx);

    sendvec[0] = xvec[localmesh->xstart];  // element from operator inverse required by neighbour
    sendvec[1] = minvb[localmesh->xstart]; // element from RHS required by neighbour

    localmesh->sendXIn(&sendvec[0],2,1);
    localmesh->wait(recv[0]);

    xk1d[localmesh->xstart-1] = ( recvec[1] + recvec[0]*minvb[localmesh->xstart] )/(1.0 - sendvec[0]*recvec[0]);

  }

  // If not on outermost boundary, get information from neighbouring proc and
  // calculate value of solution in halo cell
  if(!localmesh->lastX()) {

    comm_handle recv[1];
    recv[0] = localmesh->irecvXOut(&recvec[0], 2, 1);

    // Need the (xend+1)-th element
    for(int i=0; i<ncx; i++){
      evec[i] = 0.0;
    }
    evec[localmesh->xend+1] = 1.0; 

    tridag(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(evec),
             std::begin(xvec), ncx);

    sendvec[0] = xvec[localmesh->xend];
    sendvec[1] = minvb[localmesh->xend];

    localmesh->sendXOut(&sendvec[0],2,0);
    localmesh->wait(recv[0]);

    xk1d[localmesh->xend+1] = ( recvec[1] + recvec[0]*minvb[localmesh->xend] )/(1.0 - sendvec[0]*recvec[0]);

  }
}

/// Check whether matrix is diagonally dominant, i.e. whether for every row the absolute
/// value of the diagonal element is greater-or-equal-to the sum of the absolute values
/// of the other elements. Being diagonally dominant is sufficient (but necessary) for
/// the Jacobi iteration to converge.
void LaplaceParallelTri::check_diagonal_dominance(const Array<dcomplex> &avec, const Array<dcomplex> &bvec, const Array<dcomplex> &cvec, const int ncx, const int jy, const int kz) {

    BoutReal on_diag;
    BoutReal off_diag;

    for(int i=0; i<ncx; i++){
      on_diag = abs(bvec[i]);
      off_diag = 0.0;
      if(i > 0){
	off_diag = abs(avec[i-1]);
      }
      if(i < ncx-1){
	off_diag = off_diag + abs(cvec[i]);
      }
      if( off_diag > on_diag){
	output << "Not diagonally dominant on row "<<i<<" jy "<<jy<<" kz "<<kz<<" of proc "<<BoutComm::rank()<<endl;
      }
    }
}

FieldPerp LaplaceParallelTri::solve(const FieldPerp& b) { return solve(b, b); }

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
//FieldPerp LaplaceParallelTri::solve(const FieldPerp& b, const FieldPerp& x0, const FieldPerp& b0) {
FieldPerp LaplaceParallelTri::solve(const FieldPerp& b, const FieldPerp& x0) {

  SCOREP0();
  Timer timer("invert"); ///< Start timer

///  SCOREP_USER_REGION_DEFINE(initvars);
///  SCOREP_USER_REGION_BEGIN(initvars, "init vars",SCOREP_USER_REGION_TYPE_COMMON);

  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  TRACE("LaplaceParallelTri::solve(const, const)");

  FieldPerp x{emptyFrom(b)};

  // Convergence flags
  bool self_in = false;
  bool self_out = false;
  bool neighbour_in = false;
  bool neighbour_out = false;

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
   *        LaplaceParallelTri::solve()
   * bk1d = The 1d array of bk
   * xk   = The fourier transformed of x, where x the output of
   *        LaplaceParallelTri::solve()
   * xk1d = The 1d array of xk
   */
  auto evec = Array<dcomplex>(ncx);
  auto tmp = Array<dcomplex>(ncx);
  auto upperGuardVector = Array<dcomplex>(ncx);
  auto lowerGuardVector = Array<dcomplex>(ncx);
  auto bk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto bk1d = Array<dcomplex>(ncx);
  auto bk1d_eff = Array<dcomplex>(ncx);
  auto xk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto xk1d = Array<dcomplex>(ncx);
  auto xk1dlast = Array<dcomplex>(ncx);
  auto error = Array<dcomplex>(ncx);
  BoutReal error_rel = 1e20, error_abs=1e20, last_error=error_abs;

///  SCOREP_USER_REGION_END(initvars);
///  SCOREP_USER_REGION_DEFINE(initloop);
///  SCOREP_USER_REGION_BEGIN(initloop, "init xk loop",SCOREP_USER_REGION_TYPE_COMMON);

  // Initialise xk to 0 as we only visit 0<= kz <= maxmode in solve
  for (int ix = 0; ix < ncx; ix++) {
    for (int kz = maxmode + 1; kz < ncz / 2 + 1; kz++) {
      xk(ix, kz) = 0.0;
    }
  }
///  SCOREP_USER_REGION_END(initloop);
///  SCOREP_USER_REGION_DEFINE(fftloop);
///  SCOREP_USER_REGION_BEGIN(fftloop, "init fft loop",SCOREP_USER_REGION_TYPE_COMMON);

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
  auto avec_eff = Array<dcomplex>(ncx);
  auto bvec_eff = Array<dcomplex>(ncx);
  auto cvec_eff = Array<dcomplex>(ncx);
  
  auto minvb = Array<dcomplex>(ncx);

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
      //output << "here" << endl;
      rfft(x0[ix], ncz, &bk(ix, 0));

    } else {
      // b is the input
      // bk is the output
      rfft(b[ix], ncz, &bk(ix, 0));
      //rfft(x0[ix], ncz, &xk(ix, 0));
    }
  }
///  SCOREP_USER_REGION_END(fftloop);
///  SCOREP_USER_REGION_DEFINE(mainloop);
///  SCOREP_USER_REGION_BEGIN(mainloop, "main loop",SCOREP_USER_REGION_TYPE_COMMON);

  /* Solve differential equation in x for each fourier mode
   * Note that only the non-degenerate fourier modes are being used (i.e. the
   * offset and all the modes up to the Nyquist frequency)
   */
  for (int kz = 0; kz <= maxmode; kz++) {

///    SCOREP_USER_REGION_DEFINE(kzinit);
///    SCOREP_USER_REGION_BEGIN(kzinit, "kz init",SCOREP_USER_REGION_TYPE_COMMON);
    // set bk1d
    for (int ix = 0; ix < ncx; ix++) {
      // Get bk of the current fourier mode
      bk1d[ix] = bk(ix, kz);

      //xk1d[ix] = xk(ix, kz);
      //xk1dlast[ix] = xk(ix, kz);

      //output << "start1 "<<ix<<" "<<jy<<" "<<kz<<" "<<x0saved(ix,jy,kz)<<endl;
      xk1d[ix] = x0saved(ix, jy, kz);
    }

    int count = 0;

    // Set all convergence flags to false
    self_in = false;
    self_out = false;
    neighbour_in = false;
    neighbour_out = false;

    // Boundary values are "converged" at the start
    // Note: set neighbour's flag (not self's flag) to ensure we do at least
    // one iteration
    if(localmesh->lastX()) { 
      neighbour_out = true;
    }
    if(localmesh->firstX()) { 
      neighbour_in = true;
    }

//      for(int ix = 0; ix<localmesh->LocalNx ; ix++) {
//	output << "before "<<BoutComm::rank()<<" "<<ix<<" "<<jy<<" "<<kz<<" "<<A(ix,jy)<<" "<<bvec[ix]<<endl;//<<" "<<B(ix,kz)<<" "<<C(ix,kz)<<" "<<avec[ix] << " " << bvec[ix] << " " << cvec[ix] << endl;
//      }

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

//      for(int ix = 0; ix<localmesh->LocalNx ; ix++) {
//	output << "after "<<BoutComm::rank()<<" "<<ix<<" "<<jy<<" "<<kz<<" "<<A(ix,jy)<<" "<<bvec[ix]<<endl;//<<" "<<B(ix,kz)<<" "<<C(ix,kz)<<" "<<avec[ix] << " " << bvec[ix] << " " << cvec[ix] << endl;
//      }


    ///////// PERFORM INVERSION /////////
    if (!localmesh->periodicX) {

      // Call tridiagonal solver
      //for(int it = 0; it < maxits; it++){ 
      BoutReal error_last = 1e20;
///      int sub_it = 0;
///      auto lh = Matrix<dcomplex>(3,ncx);
///      auto rh = Matrix<dcomplex>(3,ncx);

///      if( first_call ) {
///	B = Borig;
///      }
///      else {
///	B = Bvals(0,jy,kz);
///      }
      bool allow_B_change = true;

      if( !first_call(jy,0) ) allow_B_change = false;
      //output << first_call << " " << allow_B_change << endl;

      // Patch up internal boundaries
      if(not localmesh->lastX()) { 
	for(int ix = localmesh->xend+1; ix<localmesh->LocalNx ; ix++) {
	  avec_eff[ix] = avec[ix];
	  bvec_eff[ix] = bvec[ix];
	  cvec_eff[ix] = cvec[ix];
	  avec[ix] = 0;
	  bvec[ix] = 1;
	  cvec[ix] = 0;
	  bk1d[ix] = 0;
	}
	avec[localmesh->xend+1] = cvec_eff[localmesh->xend];
      } 
      if(not localmesh->firstX()) { 
	for(int ix = 0; ix<localmesh->xstart ; ix++) {
	  avec_eff[ix] = avec[ix];
	  bvec_eff[ix] = bvec[ix];
	  cvec_eff[ix] = cvec[ix];
	  avec[ix] = 0;
	  bvec[ix] = 1;
	  cvec[ix] = 0;
	  bk1d[ix] = 0;
	}
	cvec[localmesh->xstart-1] = avec_eff[localmesh->xstart];
      }

      // Invert local matrices
      tridag(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(bk1d),
	   std::begin(minvb), ncx);

      // Find edge update vectors
      //
      // Upper interface (nguard vectors, hard-coded to two for now)
      if(not localmesh->lastX()) { 
	// Need the xend-th element
	for(int i=0; i<ncx; i++){
	  evec[i] = 0.0;
	}
	evec[localmesh->LocalNx-2] = 1.0;
	tridag(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(evec),
	     std::begin(tmp), ncx);
	for(int i=0; i<ncx; i++){
	  upperGuardVector[i] = tmp[i];
	}
      }

      // Lower interface (nguard vectors, hard-coded to two for now)
      if(not localmesh->firstX()) { 

	for(int i=0; i<ncx; i++){
	  evec[i] = 0.0;
	}
	evec[1] = 1.0;
	tridag(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(evec),
	     std::begin(tmp), ncx);
	for(int i=0; i<ncx; i++){
	  lowerGuardVector[i] = tmp[i];
	}
      } 

      ensure_stability(avec,bvec,cvec,minvb,ncx,xk1d);

      dcomplex lfac = xk1d[localmesh->xstart-1];
      dcomplex ufac = xk1d[localmesh->xend+1];

	for(int i=0; i<ncx; i++){
	  xk1d[i] = minvb[i];
	}

	if(not localmesh->lastX()) { 
	  //for(int i=localmesh->xstart; i<localmesh->xend+1; i++){
	  for(int i=0; i<ncx; i++){
	    xk1d[i] += upperGuardVector[i]*ufac;
	  }
	  xk1d[localmesh->xend+1] = ufac;
	}

	if(not localmesh->firstX()) { 
	  //for(int i=localmesh->xstart; i<localmesh->xend+1; i++){
	  for(int i=0; i<ncx; i++){
	    xk1d[i] += lowerGuardVector[i]*lfac;
	  }
	  xk1d[localmesh->xstart-1] = lfac;
	} 

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

///    SCOREP_USER_REGION_DEFINE(afterloop);
///    SCOREP_USER_REGION_BEGIN(afterloop, "after faff",SCOREP_USER_REGION_TYPE_COMMON);
    ++ncalls;
    ipt_mean_its = (ipt_mean_its * BoutReal(ncalls-1)
	+ BoutReal(count))/BoutReal(ncalls);
    //output<<"jy="<<jy<<" kz="<<kz<<" count="<<count<<" ncalls="<<ncalls<<" ipt_mean_its="<<ipt_mean_its<<" B="<<B<< endl;
    //Bvals(0,jy,kz) = B;

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
      x0saved(ix, jy, kz) = xk(ix, kz);
    }
///    SCOREP_USER_REGION_END(afterloop);
  }
///  SCOREP_USER_REGION_END(mainloop);

  //std::cout<<"end"<<endl;

  // Done inversion, transform back
  for (int ix = 0; ix < ncx; ix++) {

    if(global_flags & INVERT_ZERO_DC)
      xk(ix, 0) = 0.0;

    irfft(&xk(ix, 0), ncz, x[ix]);

#if CHECK > 2
    for(int kz=0;kz<ncz;kz++)
      if(!finite(x(ix,kz)))
        throw BoutException("Non-finite at %d, %d, %d", ix, jy, kz);
#endif
  }

  //if( first_call ){
  //  bout::globals::dump.add(Bvals, "exponents", false);
  //}
  first_call(jy,0) = false;

  return x; // Result of the inversion
}
