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
  OPTION(opt, omega, 0.0);
  OPTION(opt, new_method, false);

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

/*
 * Get an initial guess for the solution x by solving the system neglecting
 * coupling terms. This may be considered a form of preconditioning.
 * Note that coupling terms are not neglected when they are known from the
 * boundary conditions; consequently this gives the exact solution when using
 * two processors.
 */
void LaplaceParallelTri::get_initial_guess(const int jy, const int kz, Array<dcomplex> &minvb,
					      Tensor<dcomplex> &lowerGuardVector, Tensor<dcomplex> &upperGuardVector,
					      Array<dcomplex> &xk1d) {

SCOREP0();

  int xs = localmesh->xstart;
  int xe = localmesh->xend;
  xk1d[xs] = minvb[xs]/(1.0+upperGuardVector(xs,jy,kz)+lowerGuardVector(xs,jy,kz));
  xk1d[xe] = minvb[xe]/(1.0+upperGuardVector(xe,jy,kz)+lowerGuardVector(xe,jy,kz));

  /*
  int ncx = localmesh->LocalNx;
  int xs = localmesh->xstart;
  int xe = localmesh->xend;
  Array<dcomplex> sendvec, recvec;
  sendvec = Array<dcomplex>(2);
  recvec = Array<dcomplex>(2);

  // If not on innermost boundary, get information from neighbouring proc and
  // calculate value of solution in halo cell
  if(!localmesh->firstX()) {

    comm_handle recv[1];
    recv[0] = localmesh->irecvXIn(&recvec[0], 2, 0);

    sendvec[0] = lowerGuardVector(localmesh->xstart,jy,kz);  // element from operator inverse required by neighbour
    sendvec[1] = minvb[localmesh->xstart]; // element from RHS required by neighbour
    // If last processor, include known boundary terms
    if(localmesh->lastX()) {
      sendvec[1] += lowerGuardVector(localmesh->xstart,jy,kz)*xk1d[localmesh->xend+1];
    }

    localmesh->sendXIn(&sendvec[0],2,1);
    localmesh->wait(recv[0]);

    xk1d[localmesh->xstart-1] = ( recvec[1] + recvec[0]*minvb[localmesh->xstart] )/(1.0 - sendvec[0]*recvec[0]);

  }

  // If not on outermost boundary, get information from neighbouring proc and
  // calculate value of solution in halo cell
  if(!localmesh->lastX()) {

    comm_handle recv[1];
    recv[0] = localmesh->irecvXOut(&recvec[0], 2, 1);

    sendvec[0] = upperGuardVector(localmesh->xend,jy,kz);
    sendvec[1] = minvb[localmesh->xend];
    // If first processor, include known boundary terms
    if(localmesh->firstX()) {
      sendvec[1] += upperGuardVector(localmesh->xend,jy,kz)*xk1d[localmesh->xstart-1];
    }

    localmesh->sendXOut(&sendvec[0],2,0);
    localmesh->wait(recv[0]);

    xk1d[localmesh->xend+1] = ( recvec[1] + recvec[0]*minvb[localmesh->xend] )/(1.0 - sendvec[0]*recvec[0]);

  }

  for(int i=xs; i<xe+1; i++){
    xk1d[i] = minvb[i];
  }
  if(not localmesh->lastX()) {
    for(int i=xs; i<xe+1; i++){
      xk1d[i] += upperGuardVector(i,jy,kz)*xk1d[xe+1];
    }
  }
  if(not localmesh->firstX()) {
    for(int i=xs; i<xe+1; i++){
      xk1d[i] += lowerGuardVector(i,jy,kz)*xk1d[xs-1];
    }
  }
 */
}

/// Check whether matrix is diagonally dominant, i.e. whether for every row the absolute
/// value of the diagonal element is greater-or-equal-to the sum of the absolute values
/// of the other elements. Being diagonally dominant is sufficient (but not necessary) for
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

/// Check whether the reduced matrix is diagonally dominant, i.e. whether for every row the absolute
/// value of the diagonal element is greater-or-equal-to the sum of the absolute values
/// of the other elements. Being diagonally dominant is sufficient (but not necessary) for
/// the Jacobi iteration to converge.
bool LaplaceParallelTri::is_diagonally_dominant(const dcomplex al, const dcomplex au, const dcomplex bl, const dcomplex bu, const int jy, const int kz) {

  bool is_dd = true;
  if(std::fabs(al)+std::fabs(bl)>1.0){
    output<<BoutComm::rank()<<" jy="<<jy<<", kz="<<kz<<", lower row not diagonally dominant"<<endl;
    is_dd = false;
  }
  if(std::fabs(au)+std::fabs(bu)>1.0){
    output<<BoutComm::rank()<<" jy="<<jy<<", kz="<<kz<<", upper row not diagonally dominant"<<endl;
    is_dd = false;
  }
  return is_dd;
}

void get_errors(BoutReal *error_rel,BoutReal *error_abs,const dcomplex x,const dcomplex xlast){
  *error_abs = abs(x - xlast);
  BoutReal xabs = fabs(x);
  if( xabs > 0.0 ){
    *error_rel = *error_abs / xabs;
  }
  else{
    *error_rel = *error_abs;
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

  SCOREP_USER_REGION_DEFINE(initvars);
  SCOREP_USER_REGION_BEGIN(initvars, "init vars",SCOREP_USER_REGION_TYPE_COMMON);

  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  TRACE("LaplaceParallelTri::solve(const, const)");

  FieldPerp x{emptyFrom(b)};

  // Calculation variables
  // proc:       p-1   |          p          |       p+1
  // xloc:     xloc[0] | xloc[1]     xloc[2] | xloc[3]    ...
  // In this method, each processor solves equations on its processor
  // interfaces.  Its lower interface equation (for xloc[1]) is coupled to
  // xloc[0] on the processor below and its upper interface variable xloc[2].
  // Its upper interface equation (for xloc[2]) is coupled to its lower
  // interface variable xloc[1] and xloc[3] processor above.
  // We use these local variables rather than calculate with xk1d directly
  // as the elements of xloc can refer to different elements of xk1d depending
  // on the method used.
  // For example, in the original iteration we have:
  // xloc[0] = xk1d[xstart-1], xloc[1] = xk1d[xstart],
  // xloc[2] = xk1d[xend],     xloc[3] = xk1d[xend+1],
  // but if this is found to be unstable, he must change this to
  // xloc[0] = xk1d[xstart], xloc[1] = xk1d[xstart-1],
  // xloc[2] = xk1d[xend+1], xloc[3] = xk1d[xend].
  // It is easier to change the meaning of local variables and keep the
  // structure of the calculation/communication than it is to change the
  // indexing of xk1d to cover all possible cases.
  //
  auto xloc = Array<dcomplex>(4);
  auto xloclast = Array<dcomplex>(4);
  dcomplex rl, ru, al, au, bl, bu, cl, cu;
  dcomplex alold, auold, blold, buold, rlold, ruold;

  // Convergence flags
  bool self_in = false;
  bool self_out = false;
  bool neighbour_in = false;
  bool neighbour_out = false;

  int jy = b.getIndex();
  int ny = b.getMesh()->LocalNy;

  int ncz = localmesh->LocalNz; // No of z pnts
  int ncx = localmesh->LocalNx; // No of x pnts

  int xs = localmesh->xstart;
  int xe = localmesh->xend;

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
  auto upperGuardVector = Tensor<dcomplex>(ncx, ny, ncz / 2 + 1);
  auto lowerGuardVector = Tensor<dcomplex>(ncx, ny, ncz / 2 + 1);
  auto bk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto bk1d = Array<dcomplex>(ncx);
  auto bk1d_eff = Array<dcomplex>(ncx);
  auto xk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto xk1d = Array<dcomplex>(ncx);
  auto xk1dlast = Array<dcomplex>(ncx);
  auto error = Array<dcomplex>(ncx);
  dcomplex tmp2;
  BoutReal error_rel = 1e20, error_abs=1e20, last_error=error_abs;
  BoutReal error_rel_lower = 1e20, error_abs_lower=1e20;
  BoutReal error_rel_upper = 1e20, error_abs_upper=1e20;
  BoutReal error_rel_lower_last = 1e20, error_abs_lower_last=1e20;
  BoutReal error_rel_upper_last = 1e20, error_abs_upper_last=1e20;
  BoutReal error_rel_lower_two_old = 1e20, error_abs_lower_two_old=1e20;
  BoutReal error_rel_upper_two_old = 1e20, error_abs_upper_two_old=1e20;
  bool lower_stable = false;
  bool upper_stable = false;
  // Down and up coefficients
  dcomplex Bd, Ad, Rd;
  dcomplex Bu, Au, Ru;
  dcomplex Btmp, Atmp, Rtmp;

  SCOREP_USER_REGION_END(initvars);
  SCOREP_USER_REGION_DEFINE(initloop);
  SCOREP_USER_REGION_BEGIN(initloop, "init xk loop",SCOREP_USER_REGION_TYPE_COMMON);

  // Initialise xk to 0 as we only visit 0<= kz <= maxmode in solve
  for (int ix = 0; ix < ncx; ix++) {
    for (int kz = maxmode + 1; kz < ncz / 2 + 1; kz++) {
      xk(ix, kz) = 0.0;
    }
  }
  SCOREP_USER_REGION_END(initloop);
  SCOREP_USER_REGION_DEFINE(fftloop);
  SCOREP_USER_REGION_BEGIN(fftloop, "init fft loop",SCOREP_USER_REGION_TYPE_COMMON);

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
  auto minvb = Array<dcomplex>(ncx);
  bool lowerUnstable = false;
  bool upperUnstable = false;

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
  SCOREP_USER_REGION_END(fftloop);
  SCOREP_USER_REGION_DEFINE(mainloop);
  SCOREP_USER_REGION_BEGIN(mainloop, "main loop",SCOREP_USER_REGION_TYPE_COMMON);

  /* Solve differential equation in x for each fourier mode
   * Note that only the non-degenerate fourier modes are being used (i.e. the
   * offset and all the modes up to the Nyquist frequency)
   */
  for (int kz = 0; kz <= maxmode; kz++) {

    SCOREP_USER_REGION_DEFINE(kzinit);
    SCOREP_USER_REGION_BEGIN(kzinit, "kz init",SCOREP_USER_REGION_TYPE_COMMON);
    // set bk1d
    for (int ix = 0; ix < ncx; ix++) {
      // Get bk of the current fourier mode
      bk1d[ix] = bk(ix, kz);

      //xk1d[ix] = xk(ix, kz);
      //xk1dlast[ix] = xk(ix, kz);

      //output << "start1 "<<ix<<" "<<jy<<" "<<kz<<" "<<x0saved(ix,jy,kz)<<endl;
      xk1d[ix] = x0saved(ix, jy, kz);
      //output << "start2 "<<ix<<" "<<jy<<" "<<kz<<endl;
      xk1dlast[ix] = x0saved(ix, jy, kz);
      //output << "start3 "<<ix<<" "<<jy<<" "<<kz<<endl;
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
	  avec[ix] = 0;
	  bvec[ix] = 1;
	  cvec[ix] = 0;
	  bk1d[ix] = 0;
	}
      } 
      if(not localmesh->firstX()) { 
	for(int ix = 0; ix<localmesh->xstart ; ix++) {
	  avec[ix] = 0;
	  bvec[ix] = 1;
	  cvec[ix] = 0;
	  bk1d[ix] = 0;
	}
      }

	SCOREP_USER_REGION_END(kzinit);
	SCOREP_USER_REGION_DEFINE(invert);
	SCOREP_USER_REGION_BEGIN(invert, "invert local matrices",SCOREP_USER_REGION_TYPE_COMMON);
      // Invert local matrices
      // Calculate Minv*b
      tridag(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(bk1d),
	   std::begin(minvb), ncx);
      // Now minvb is a constant vector throughout the iterations
      //
      // Find edge update vectors
      //
      // Upper interface (nguard vectors, hard-coded to two for now)
      if(not localmesh->lastX()) { 
	// Need the xend-th element
	for(int i=0; i<ncx; i++){
	  evec[i] = 0.0;
	}
	evec[localmesh->xend+1] = 1.0;
	tridag(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(evec),
	     std::begin(tmp), ncx);
	for(int i=0; i<ncx; i++){
	  upperGuardVector(i,jy,kz) = tmp[i];
	}
      } else {
	for(int i=0; i<ncx; i++){
	  upperGuardVector(i,jy,kz) = 0.0;
	}
      }

      // Lower interface (nguard vectors, hard-coded to two for now)
      if(not localmesh->firstX()) { 

	for(int i=0; i<ncx; i++){
	  evec[i] = 0.0;
	}
	evec[localmesh->xstart-1] = 1.0;
	tridag(std::begin(avec), std::begin(bvec), std::begin(cvec), std::begin(evec),
	     std::begin(tmp), ncx);
	for(int i=0; i<ncx; i++){
	  lowerGuardVector(i,jy,kz) = tmp[i];
	}
      } else {
	for(int i=0; i<ncx; i++){
	  lowerGuardVector(i,jy,kz) = 0.0;
	}
      }

	SCOREP_USER_REGION_END(invert);
	SCOREP_USER_REGION_DEFINE(coefs);
	SCOREP_USER_REGION_BEGIN(coefs, "calculate coefs",SCOREP_USER_REGION_TYPE_COMMON);

      //check_diagonal_dominance(avec,bvec,cvec,ncx,jy,kz);
      //get_initial_guess(jy,kz,minvb,lowerGuardVector,upperGuardVector,xk1d);

      // Original method:
      if(not lowerUnstable) {
	xloclast[0] = xk1d[xs-1];
	xloclast[1] = xk1d[xs];
	rl = minvb[xs];
	bl = upperGuardVector(xs,jy,kz);
	al = lowerGuardVector(xs,jy,kz);
      } else {
	xloclast[0] = xk1d[xs];
	xloclast[1] = xk1d[xs-1];
	rl = -minvb[xs]/lowerGuardVector(xs,jy,kz);
	bl = -upperGuardVector(xs,jy,kz)/lowerGuardVector(xs,jy,kz);
	al = 1.0/lowerGuardVector(xs,jy,kz);
      }
      if(not upperUnstable) {
	xloclast[2] = xk1d[xe];
	xloclast[3] = xk1d[xe+1];
	ru = minvb[xe];
	bu = upperGuardVector(xe,jy,kz);
	au = lowerGuardVector(xe,jy,kz);
      } else {
	xloclast[2] = xk1d[xe+1];
	xloclast[3] = xk1d[xe];
	ru = -minvb[xe]/upperGuardVector(xe,jy,kz);
	bu = 1.0/upperGuardVector(xe,jy,kz);
	au = -lowerGuardVector(xe,jy,kz)/upperGuardVector(xe,jy,kz);
      }

      alold = al;
      auold = au;
      blold = bl;
      buold = bu;
      rlold = rl;
      ruold = ru;

      // New method - connect to more distant points
      if(new_method){
      xloclast[0] = xk1d[xs-1];
      xloclast[1] = xk1d[xs];
      xloclast[2] = xk1d[xe];
      xloclast[3] = xk1d[xe+1];

      if(not localmesh->firstX()){
	// Send coefficients down
	Rtmp = rl;
	Atmp = al;
	Btmp = 0.0;
	if( std::fabs(bu) > 1e-14 ){
	  Btmp = bl/bu;
	  Atmp -= Btmp*au;
	  Rtmp -= Btmp*ru;
	}
	// Send these
	Ad = localmesh->communicateXIn(Atmp);
	Bd = localmesh->communicateXIn(Btmp);
	Rd = localmesh->communicateXIn(Rtmp);
      }
      if(not localmesh->lastX()){
	// Send coefficients up
	Rtmp = ru;
	Atmp = 0.0;
	Btmp = bu;
	if( std::fabs(al) > 1e-14 ){
	  Atmp = au/al;
	  Btmp -= Atmp*bl;
	  Rtmp -= Atmp*rl;
	}
	// Send these
	Au = localmesh->communicateXOut(Atmp);
	Bu = localmesh->communicateXOut(Btmp);
	Ru = localmesh->communicateXOut(Rtmp);
      }

      dcomplex Delta;


      if(localmesh->firstX()){
	Ad = 1.0;
	Rd = 0.0;
	Bd = 0.0;
      }
      if(localmesh->lastX()){
	Au = 0.0;
	Ru = 0.0;
	Bu = 1.0;
      }
      Delta = 1.0 - al*Bd - bu*Au + (al*bu - au*bl)*Bd*Au;
      Delta = 1.0 / Delta;
      rl = Delta*( (rlold + alold*Rd + blold*Ru) 
	         + (ruold*blold - rlold*buold)*Au
		 + (auold*blold - alold*buold)*Au*Rd );
      al = Delta*( alold + (auold*blold - alold*buold)*Au )*Ad;	
      bl = Delta * blold * Bu ;	
      ru = Delta*( (ruold + auold*Rd + buold*Ru)
	         + (rlold*auold - ruold*alold)*Bd
		 + (auold*blold - alold*buold)*Bd*Ru );
      bu = Delta*( buold + (auold*blold - alold*buold)*Bd )*Bu;
      au = Delta * auold * Ad ;

///      if(jy==0 and kz==0){
///	output<<alold<<" "<<blold<<" "<<auold<<" "<<buold<<endl;
///	output<<al<<" "<<bl<<" "<<au<<" "<<bu<<endl;
//////	for(int ix=0; ix<ncx;ix++){
/////////	  output<<avec[ix]<<" "<<bvec[ix]<<" "<<cvec[ix]<<endl;
//////	  output<<lowerGuardVector(ix,jy,kz)<<" "<<upperGuardVector(ix,jy,kz)<<" "<<endl;
//////	}
///	output<<" "<<endl;
///      }

      }

	SCOREP_USER_REGION_END(coefs);
      //if(jy==0 and kz==1){
      //output<<"Coefficients: "<<BoutComm::rank()<<" "<<jy<<" "<<kz<<" "<<" "<<rl<<" "<<al<<" "<<bl<<" "<<ru<<" "<<au<<" "<<bu<<endl;
      //output<<"Coefficients: "<<BoutComm::rank()<<" "<<jy<<" "<<kz<<" "<<" "<<al<<" "<<bl<<" "<<au<<" "<<bu<<endl;
      //}
      //output<<"xvec "<<BoutComm::rank()<<" "<<"initial"<<" "<<xloc[0]<<" "<<xloc[1]<<" "<<xloc[2]<<" "<<xloc[3]<<" "<<xloclast[0]<<" "<<xloclast[1]<<" "<<xloclast[2]<<" "<<xloclast[3]<<" "<<error_rel_lower<<" "<<error_rel_lower_last<<" "<<error_rel_lower_two_old<<" "<<error_abs_lower<<" "<<error_abs_lower_last<<" "<<error_abs_lower_two_old<<" "<<error_rel_upper<<" "<<error_rel_upper_last<<" "<<error_rel_upper_two_old<<" "<<error_abs_upper<<" "<<error_abs_upper_last<<" "<<error_abs_upper_two_old<<endl;

///      SCOREP_USER_REGION_END(kzinit);
      SCOREP_USER_REGION_DEFINE(whileloop);
      SCOREP_USER_REGION_BEGIN(whileloop, "while loop",SCOREP_USER_REGION_TYPE_COMMON);
//
      BoutReal om = 0.0;
      if(kz==0) om = omega;

      int upper_offset = 0;
      int lower_offset = 0;
      while(true){ 

	SCOREP_USER_REGION_DEFINE(iteration);
	SCOREP_USER_REGION_BEGIN(iteration, "iteration",SCOREP_USER_REGION_TYPE_COMMON);

	// Only need to update interior points
	xloc[1] = rl;
	xloc[2] = ru;
	if(not localmesh->lastX()) {
	  xloc[1] += bl*xloclast[3];
	  xloc[2] += bu*xloclast[3];
	}
	if(not localmesh->firstX()) {
	  xloc[1] += al*xloclast[0];
	  xloc[2] += au*xloclast[0];
	}

	//xloc[1] = (1.0-om)*xloc[1] + om*xloclast[1];
	//xloc[2] = (1.0-om)*xloc[2] + om*xloclast[2];

	SCOREP_USER_REGION_END(iteration);
	SCOREP_USER_REGION_DEFINE(comms);
	SCOREP_USER_REGION_BEGIN(comms, "communication",SCOREP_USER_REGION_TYPE_COMMON);

	TRACE("set comm flags pack");
	// Set communication flags
	if ( count > 2 and
            (kz==0 or (
	     (error_rel_lower<rtol or error_abs_lower<atol) and
             (error_rel_upper<rtol or error_abs_upper<atol) ))) {
	  // In the next iteration this proc informs its neighbours that its halo cells
	  // will no longer be updated, then breaks.
	  self_in = true;
	  self_out = true;
	}

	// Communication
	// A proc is finished when it is both in- and out-converged.
	// Once this happens, that processor communicates once, then breaks.
	//
	// A proc can be converged in only one of the directions. This happens
	// if it has not met the error tolerance, but one of its neighbours has
	// converged. In this case, that boundary will no longer update (and 
	// communication in that direction should stop), but the other boundary
	// may still be changing.
	//
	// There are four values to consider:
	//   neighbour_in  = whether my in-neighbouring proc has out-converged
	//   self_in       = whether this proc has in-converged
	//   self_out      = whether this proc has out-converged
	//   neighbour_out = whether my out-neighbouring proc has in-converged
	//
	// If neighbour_in = true, I must have been told this by my neighbour. My
	// neighbour has therefore done its one post-converged communication. My in-boundary
	// values are therefore correct, and I am in-converged. My neighbour is not 
	// expecting us to communicate.
	if(!neighbour_in) {
	  // Communicate in
	  neighbour_in = localmesh->communicateXIn(self_in);
	  if(new_method){
	    xloc[0] = localmesh->communicateXIn(xloc[2]);
	  }
	  else{
	    xloc[0] = localmesh->communicateXIn(xloc[1]);
	  }
	  //output<<BoutComm::rank()<<" "<<xloc[0]<<" "<<xloc[1]<<" "<<xloc[2]<<" "<<xloc[3]<<endl;
	}

	// Outward communication
	// See note above for inward communication.
	if(!neighbour_out) {
	  // Communicate out
	  neighbour_out = localmesh->communicateXOut(self_out);
	  if(new_method){
	    xloc[3] = localmesh->communicateXOut(xloc[1]);
	  }
	  else{
	    xloc[3] = localmesh->communicateXOut(xloc[2]);
	  }
	}
	SCOREP_USER_REGION_END(comms);

	// Now I've done my communication, exit if I am both in- and out-converged
	if( self_in and self_out ) {
          //output<<"Breaking, proc "<< BoutComm::rank() << ", count "<<count<<" "<<jy<<" "<<kz<<endl<<std::flush;
	  break;
	}
	SCOREP_USER_REGION_DEFINE(comms_after_break);
	SCOREP_USER_REGION_BEGIN(comms_after_break, "comms after break",SCOREP_USER_REGION_TYPE_COMMON);

	// If my neighbour has converged, I know that I am also converged on that
	// boundary. Set this flag after the break loop above, to ensure we do one
	// iteration using our neighbour's converged value.
	if(neighbour_in) {
	  self_in = true;
	}
	if(neighbour_out) {
	  self_out = true;
	}

	++count;
	SCOREP_USER_REGION_END(comms_after_break);
	if (count>maxits) {
	    output<<alold<<" "<<blold<<" "<<auold<<" "<<buold<<endl;
	    output<<al<<" "<<bl<<" "<<au<<" "<<bu<<endl;
	  //break;
	  // Maximum number of allowed iterations reached.
	  // If the iteration matrix is diagonally-dominant, then convergence is guaranteed, so maxits is set too low.
	  // Otherwise, the method may or may not converge.
          if(is_diagonally_dominant(al,au,bl,bu,jy,kz)){
	    throw BoutException("LaplaceParallelTri error: Not converged within maxits=%i iterations. The iteration matrix is diagonally dominant on processor %i and convergence is guaranteed (if all processors are diagonally dominant). Please increase maxits and retry.",maxits,BoutComm::rank());
	  }
	  else{
	    output<<alold<<" "<<blold<<" "<<auold<<" "<<buold<<endl;
	    output<<al<<" "<<bl<<" "<<au<<" "<<bu<<endl;
	    output<<Ad<<" "<<Bd<<" "<<Au<<" "<<Bu<<endl;
	    throw BoutException("LaplaceParallelTri error: Not converged within maxits=%i iterations. The iteration matrix is not diagonally dominant on processor %i, so there is no guarantee this method will converge. Consider increasing maxits or using a different solver.",maxits,BoutComm::rank());
	  }
	}

	SCOREP_USER_REGION_DEFINE(errors);
	SCOREP_USER_REGION_BEGIN(errors, "calculate errors",SCOREP_USER_REGION_TYPE_COMMON);
//

	// Calculate errors
	error_abs_lower = 0.0;
	error_abs_upper = 0.0;
	error_rel_lower = 0.0;
	error_rel_upper = 0.0;

	// Calcalate errors on left halo and right interior point - this means the
	// errors on neighbouring processors agree exactly without the need for
	// communication.
	get_errors(&error_rel_lower,&error_abs_lower,xloc[1],xloclast[1]);
	get_errors(&error_rel_upper,&error_abs_upper,xloc[2],xloclast[2]);
	
        //output<<"xvec "<<BoutComm::rank()<<" "<<count<<" "<<xloc[0]<<" "<<xloc[1]<<" "<<xloc[2]<<" "<<xloc[3]<<" "<<xloclast[0]<<" "<<xloclast[1]<<" "<<xloclast[2]<<" "<<xloclast[3]<<" "<<error_rel_lower<<" "<<error_rel_lower_last<<" "<<error_rel_lower_two_old<<" "<<error_abs_lower<<" "<<error_abs_lower_last<<" "<<error_abs_lower_two_old<<" "<<error_rel_upper<<" "<<error_rel_upper_last<<" "<<error_rel_upper_two_old<<" "<<error_abs_upper<<" "<<error_abs_upper_last<<" "<<error_abs_upper_two_old<<endl;
	SCOREP_USER_REGION_END(errors);

	SCOREP_USER_REGION_DEFINE(copylast);
	SCOREP_USER_REGION_BEGIN(copylast, "copy to last",SCOREP_USER_REGION_TYPE_COMMON);
        for (int ix = 0; ix < 4; ix++) {
	  xloclast[ix] = xloc[ix];
	}
	SCOREP_USER_REGION_END(copylast);
	
      }
      SCOREP_USER_REGION_END(whileloop);

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

    SCOREP_USER_REGION_DEFINE(afterloop);
    SCOREP_USER_REGION_BEGIN(afterloop, "after faff",SCOREP_USER_REGION_TYPE_COMMON);
    ++ncalls;
    ipt_mean_its = (ipt_mean_its * BoutReal(ncalls-1)
	+ BoutReal(count))/BoutReal(ncalls);

    //if(jy==0 and kz==1){
    //output<<"jy="<<jy<<" kz="<<kz<<" count="<<count<<" ncalls="<<ncalls<<" ipt_mean_its="<<ipt_mean_its<<" B="<<B<< endl;
    //}
    //Bvals(0,jy,kz) = B;

    // Original method:
    if(not lowerUnstable) {
      xk1d[xs-1] = xloc[0];
      xk1d[xs]   = xloc[1];
      xk1dlast[xs-1] = xloclast[0];
      xk1dlast[xs]   = xloclast[1];
    } else {
      xk1d[xs-1] = xloc[1];
      xk1d[xs]   = xloc[0];
      xk1dlast[xs-1] = xloclast[1];
      xk1dlast[xs]   = xloclast[0];
    }
    if(not upperUnstable) {
      xk1d[xe]   = xloc[2];
      xk1d[xe+1] = xloc[3];
      xk1dlast[xe]   = xloclast[2];
      xk1dlast[xe+1] = xloclast[3];
    } else {
      xk1d[xe]   = xloc[3];
      xk1d[xe+1] = xloc[2];
      xk1dlast[xe]   = xloclast[3];
      xk1dlast[xe+1] = xloclast[2];
    }

    if(new_method){
    dcomplex d = 1.0/(buold*alold - blold*auold);
    // If boundary processor, halo cell is already correct, and d is undefined.
    // Lower boundary proc => al = au = 0
    // Upper boundary proc => bl = bu = 0
    if(not localmesh->firstX() and not localmesh->lastX()){
      // General case
      xk1dlast[xs-1] =  d*(buold*(xk1dlast[xs]-rlold) - blold*(xk1dlast[xe]-ruold));
      xk1dlast[xe+1] = -d*(auold*(xk1dlast[xs]-rlold) - alold*(xk1dlast[xe]-ruold));
    } else if(localmesh->firstX() and not localmesh->lastX()) {
      // Lower boundary but not upper boundary
      // xk1dlast[xs-1] = already correct
      xk1dlast[xe+1] = (xk1dlast[xe]-ruold)/buold;
    } else if(localmesh->lastX() and not localmesh->firstX()){
      // Upper boundary but not lower boundary
      // xk1dlast[xe+1] = already correct
      xk1dlast[xs-1] = (xk1dlast[xs]-rlold)/alold;
    } 
    // No "else" case. If both upper and lower boundaries, both xs-1 and xe+1
    // are already correct
    }

    //output<<"Converged "<<xk1dlast[xs-1]<<" "<<xk1dlast[xs]<<" "<<xk1dlast[xe]<<" "<<xk1dlast[xe+1]<<endl;

    // Now that halo cells are converged, use these to calculate whole solution
    for(int i=0; i<ncx; i++){
      xk1d[i] = minvb[i];
    }
    if(not localmesh->lastX()) { 
      for(int i=0; i<ncx; i++){
        xk1d[i] += upperGuardVector(i,jy,kz)*xk1dlast[xe+1];
      }
    }
    if(not localmesh->firstX()) { 
      for(int i=0; i<ncx; i++){
        xk1d[i] += lowerGuardVector(i,jy,kz)*xk1dlast[xs-1];
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
      x0saved(ix, jy, kz) = xk(ix, kz);
    }
    SCOREP_USER_REGION_END(afterloop);
  }
  SCOREP_USER_REGION_END(mainloop);

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
