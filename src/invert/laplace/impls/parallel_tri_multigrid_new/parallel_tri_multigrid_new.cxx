/**************************************************************************
 * Perpendicular Laplacian inversion. Parallel code using FFTs in z
 * and an iterative tridiagonal solver in x.
 *
 **************************************************************************
 * Copyright 2020 Joseph Parker
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
#include "parallel_tri_multigrid_new.hxx"

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

LaplaceParallelTriMGNew::LaplaceParallelTriMGNew(Options *opt, CELL_LOC loc, Mesh *mesh_in)
    : Laplacian(opt, loc, mesh_in), A(0.0), C(1.0), D(1.0), ipt_mean_its(0.), ncalls(0) {
  A.setLocation(location);
  C.setLocation(location);
  D.setLocation(location);

  OPTION(opt, rtol, 1.e-7);
  OPTION(opt, atol, 1.e-20);
  OPTION(opt, maxits, 100);
  OPTION(opt, max_level, 3);
  OPTION(opt, max_cycle, 3);
  OPTION(opt, new_method, false);
  OPTION(opt, use_previous_timestep, false);
  OPTION(opt, algorithm, 0);

  static int ipt_solver_count = 1;
  bout::globals::dump.addRepeat(ipt_mean_its,
      "ipt_solver"+std::to_string(ipt_solver_count)+"_mean_its");
  ++ipt_solver_count;

  first_call = Matrix<bool>(localmesh->LocalNy,localmesh->LocalNz / 2 + 1);

  x0saved = Tensor<dcomplex>(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz / 2 + 1);

  levels = std::vector<Level>(10);

  resetSolver();

}

/*
 * Reset the solver to its initial state
 */
void LaplaceParallelTriMGNew::resetSolver(){
  x0saved = 0.0;
  for(int jy=0; jy<localmesh->LocalNy; jy++){
    for(int kz=0; kz<localmesh->LocalNz / 2 + 1; kz++){
      first_call(jy,kz) = true;
    }
  }
  resetMeanIterations();
}

/*
 * Get an initial guess for the solution x by solving the system neglecting
 * coupling terms. This may be considered a form of preconditioning.
 * Note that coupling terms are not neglected when they are known from the
 * boundary conditions; consequently this gives the exact solution when using
 * two processors.
 */
void LaplaceParallelTriMGNew::get_initial_guess(const int jy, const int kz, Matrix<dcomplex> &minvb,
					      Tensor<dcomplex> &lowerGuardVector, Tensor<dcomplex> &upperGuardVector,
					      Matrix<dcomplex> &xk1d) {

SCOREP0();

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

    sendvec[0] = lowerGuardVector(xs,jy,kz);  // element from operator inverse required by neighbour
    sendvec[1] = minvb(kz,xs); // element from RHS required by neighbour
    // If last processor, include known boundary terms
    if(localmesh->lastX()) {
      sendvec[1] += lowerGuardVector(xs,jy,kz)*xk1d(kz,xe+1);
    }

    localmesh->sendXIn(&sendvec[0],2,1);
    localmesh->wait(recv[0]);

    xk1d(kz,xs-1) = ( recvec[1] + recvec[0]*minvb(kz,xs) )/(1.0 - sendvec[0]*recvec[0]);

  }

  // If not on outermost boundary, get information from neighbouring proc and
  // calculate value of solution in halo cell
  if(!localmesh->lastX()) {

    comm_handle recv[1];
    recv[0] = localmesh->irecvXOut(&recvec[0], 2, 1);

    sendvec[0] = upperGuardVector(xe,jy,kz);
    sendvec[1] = minvb(kz,xe);
    // If first processor, include known boundary terms
    if(localmesh->firstX()) {
      sendvec[1] += upperGuardVector(xe,jy,kz)*xk1d(kz,xs-1);
    }

    localmesh->sendXOut(&sendvec[0],2,0);
    localmesh->wait(recv[0]);

    xk1d(kz,xe+1) = ( recvec[1] + recvec[0]*minvb(kz,xe) )/(1.0 - sendvec[0]*recvec[0]);

  }

  for(int i=0; i<localmesh->LocalNx; i++){
    xk1d(kz,i) = minvb(kz,i);
  }
  if(not localmesh->lastX()) {
    for(int i=0; i<localmesh->LocalNx; i++){
      xk1d(kz,i) += upperGuardVector(i,jy,kz)*xk1d(kz,xe+1);
    }
  }
  if(not localmesh->firstX()) {
    for(int i=0; i<localmesh->LocalNx; i++){
      xk1d(kz,i) += lowerGuardVector(i,jy,kz)*xk1d(kz,xs-1);
    }
  }
}

/*
 * Check whether the reduced matrix is diagonally dominant, i.e. whether for every row the absolute
 * value of the diagonal element is greater-or-equal-to the sum of the absolute values
 * of the other elements. Being diagonally dominant is sufficient (but not necessary) for
 * the Jacobi iteration to converge.
 */
bool LaplaceParallelTriMGNew::is_diagonally_dominant(const dcomplex al, const dcomplex au, const dcomplex bl, const dcomplex bu, const int jy, const int kz) {

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

/*
 * Reconstruct the full solution from the subproblem and the halo cell values
 */
void LaplaceParallelTriMGNew::reconstruct_full_solution(Level &l, const int jy, Matrix<dcomplex> &halos){
  SCOREP0();
  for (int kz = 0; kz < nmode; kz++) {
    for(int i=0; i<l.ncx; i++){
      l.soln(kz,i) = l.minvb(kz,i);
    }
    if(not localmesh->lastX()) {
      for(int i=0; i<l.ncx; i++){
        l.soln(kz,i) += l.upperGuardVector(i,jy,kz)*halos(kz,l.xe+1);
      }
    }
    if(not localmesh->firstX()) {
      for(int i=0; i<l.ncx; i++){
        l.soln(kz,i) += l.lowerGuardVector(i,jy,kz)*halos(kz,l.xs-1);
      }
    }
  }
}

/*
 * Reconstruct the full solution from the subproblem and the halo cell values.
 * We call this at the end of the iteration when a processor has three correct
 * values: xk1d at its own xs (the first interior point), and the xs of its
 * nearest neighbours. These values are in xloc(0), xloc(1) and xloc(3), or
 * xloc(0), xloc(1) and xloc(2) for the final processor.
 * Reconstructing the solution uses the cached values from the local problem,
 * and so needs my upper neighbour's xs point (which we know) and my lower
 * neighbour's xe point (which we don't). We also calculate this point locally
 * using:
 * my_xk1d(xs) = rlold(xs) + alold(xs)*lower_xk1d(xe) + blold(xs)*upper_xk1d(xs)
 */
void LaplaceParallelTriMGNew::reconstruct_full_solution(Level &l, const int jy){
  SCOREP0();

  Array<dcomplex> x_lower, x_upper, x_lower_halo, x_upper_halo;
  x_lower = Array<dcomplex>(nmode);
  x_upper = Array<dcomplex>(nmode);
  x_lower_halo = Array<dcomplex>(nmode);
  x_upper_halo = Array<dcomplex>(nmode);

  for (int kz = 0; kz < nmode; kz++) {

    // Note: Use xloclast to reconstruct interior points, as this mimicks
    // what the iteration step does, ensuring the constructed solution is
    // consistent with xloc at the current step.
    x_lower[kz] = l.xloc(0,kz);
    x_upper[kz] = l.xloc(3,kz);

    if(not localmesh->firstX()){
      x_lower[kz] = (l.xloc(1,kz)-l.rlold[kz]-l.blold(jy,kz)*l.xloc(3,kz))/l.alold(jy,kz);
    }
  }

  for (int kz = 0; kz < nmode; kz++) {
    l.soln(kz,l.xs-1) = x_lower[kz];
    for(int i=l.xs; i<l.xe+1; i++){
      l.soln(kz,i) = l.minvb(kz,i) + l.upperGuardVector(i,jy,kz)*x_upper[kz] + l.lowerGuardVector(i,jy,kz)*x_lower[kz];
    }
    l.soln(kz,l.xe+1) = x_upper[kz];
  }

  /*
  output<<"full solution ";
  for(int i=0;i<l.ncx;i++){
    output<<l.soln(1,i)<<" ";
  }
  output<<endl;
  */

}

// TODO Move to Array
/*
 * Returns true if all values of bool array are true, otherwise returns false.
 */
bool LaplaceParallelTriMGNew::all(const Array<bool> a){
  SCOREP0();
  for(int i=0; i<a.size(); i++){
    if(a[i]==false){
      return false;
    }
  }
  return true;
}

// TODO Move to Array
/*
 * Returns true if any values of bool array are true, otherwise returns false.
 */
bool LaplaceParallelTriMGNew::any(const Array<bool> a){
  SCOREP0();
  for(int i=0; i<a.size(); i++){
    if(a[i]==true){
      return true;
    }
  }
  return false;
}

// TODO Move to Array
/*
 * Returns maximum value of BoutReal array.
 */
BoutReal LaplaceParallelTriMGNew::max(const Array<BoutReal> a){
  SCOREP0();
  BoutReal maxval = a[0];
  for(int i=1; i<a.size(); i++){
    if(a[i]>maxval){
      maxval = a[i];
    }
  }
  return maxval;
}

// TODO Move to Array
/*
 * Returns location of maximum value of BoutReal array.
 */
int LaplaceParallelTriMGNew::maxloc(const Array<BoutReal> a){
  SCOREP0();
  BoutReal maxval = a[0];
  int maxloc = 0;
  for(int i=1; i<a.size(); i++){
    if(a[i]>maxval){
      maxloc = i;
      maxval = a[i];
    }
  }
  return maxloc;
}

FieldPerp LaplaceParallelTriMGNew::solve(const FieldPerp& b) { return solve(b, b); }

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
//FieldPerp LaplaceParallelTriMGNew::solve(const FieldPerp& b, const FieldPerp& x0, const FieldPerp& b0) {
FieldPerp LaplaceParallelTriMGNew::solve(const FieldPerp& b, const FieldPerp& x0) {

  SCOREP0();
  Timer timer("invert"); ///< Start timer

  ///SCOREP_USER_REGION_DEFINE(initvars);
  ///SCOREP_USER_REGION_BEGIN(initvars, "init vars",///SCOREP_USER_REGION_TYPE_COMMON);

  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  TRACE("LaplaceParallelTriMGNew::solve(const, const)");

  FieldPerp x{emptyFrom(b)};

  // Info for halo swaps
  int xproc = localmesh->getXProcIndex();
  int yproc = localmesh->getYProcIndex();
  nproc = localmesh->getNXPE();
  myproc = yproc * nproc + xproc;
  proc_in = myproc - 1;
  proc_out = myproc + 1;

  nmode = maxmode + 1;
  int jy = b.getIndex();

  int ncz = localmesh->LocalNz; // Number of local z points
  int ncx = localmesh->LocalNx; // Number of local x points

  int xs = localmesh->xstart; // First interior point
  int xe = localmesh->xend; // Last interior point

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
  // It is easier to change the meaning of local variables and keep the
  // structure of the calculation/communication than it is to change the
  // indexing of xk1d to cover all possible cases.
  //
  // In the original iteration we have:
  // xloc[0] = xk1d[xstart-1], xloc[1] = xk1d[xstart],
  // xloc[2] = xk1d[xend],     xloc[3] = xk1d[xend+1],
  //
  auto xloc = Matrix<dcomplex>(4,nmode); // Values of x at the processor interface
  auto xloclast = Matrix<dcomplex>(4,nmode); // Values of x at the processor interface from previous iteration 

  BoutReal kwaveFactor = 2.0 * PI / coords->zlength();

  // Should we store coefficients?
  store_coefficients = not (inner_boundary_flags & INVERT_AC_GRAD);
  store_coefficients = store_coefficients && not (outer_boundary_flags & INVERT_AC_GRAD);
  store_coefficients = store_coefficients && not (inner_boundary_flags & INVERT_SET);
  store_coefficients = store_coefficients && not (outer_boundary_flags & INVERT_SET);

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

  /* Allocation for
  * bk   = The fourier transformed of b, where b is one of the inputs in
  *        LaplaceParallelTriMGNew::solve()
  * bk1d = The 1d array of bk
  * xk   = The fourier transformed of x, where x the output of
  *        LaplaceParallelTriMGNew::solve()
  * xk1d = The 1d array of xk
  */
  auto bk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto bk1d = Array<dcomplex>(ncx);
  auto xk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto xk1d = Matrix<dcomplex>(ncz/2+1,ncx);
  auto xk1dlast = Matrix<dcomplex>(ncz/2+1,ncx);

  // Error interpolated onto a finer grid
  auto fine_error = Matrix<dcomplex>(4,ncz/2+1);

  // Define indexing of xloc that depends on method. Doing this now removes
  // branching in tight loops
  // index_in is the index of the x that a processor sends inwards,
  // index_out is the index of the x it sends outwards.
  // In the usual method, these are its values closest to the processor it is
  // sending to, i.e. 1 in and 2 out. In the extended domain method, we send
  // the value of x furthest from the neighbouring proc, i.e. 2 in and 1 out.
  index_in = 1;
  index_out = 2;
  if(new_method){
    index_in = 2;
    index_out = 1;
  }

  ///SCOREP_USER_REGION_END(initvars);
  ///SCOREP_USER_REGION_DEFINE(initloop);
  ///SCOREP_USER_REGION_BEGIN(initloop, "init xk loop",SCOREP_USER_REGION_TYPE_COMMON);

  // Initialise xk to 0 as we only visit 0<= kz <= maxmode in solve
  for (int ix = 0; ix < ncx; ix++) {
    // TODO This looks wrong: should it be kz<ncz?
    // This came from serial tri
    for (int kz = maxmode + 1; kz < ncz / 2 + 1; kz++) {
      xk(ix, kz) = 0.0;
    }
  }
  ///SCOREP_USER_REGION_END(initloop);
  ///SCOREP_USER_REGION_DEFINE(fftloop);
  ///SCOREP_USER_REGION_BEGIN(fftloop, "init fft loop",SCOREP_USER_REGION_TYPE_COMMON);

  /* Coefficents in the tridiagonal solver matrix
  * Following the notation in "Numerical recipes"
  * avec is the lower diagonal of the matrix
  * bvec is the diagonal of the matrix
  * cvec is the upper diagonal of the matrix
  * NOTE: Do not confuse avec, bvec and cvec with the A, C, and D coefficients
  *       above
  */
  auto avec = Matrix<dcomplex>(nmode,ncx);
  auto bvec = Matrix<dcomplex>(nmode,ncx);
  auto cvec = Matrix<dcomplex>(nmode,ncx);
  auto bcmplx = Matrix<dcomplex>(nmode,ncx);

  BOUT_OMP(parallel for)
  for (int ix = 0; ix < ncx; ix++) {
    /* This for loop will set the bk (initialized by the constructor)
    * bk is the z fourier modes of b in z
    * If the INVERT_SET flag is set (meaning that x0 will be used to set the
    * bounadry values),
    */
    if (((ix < inbndry) && (inner_boundary_flags & INVERT_SET) && localmesh->firstX()) ||
    ((ncx - ix - 1 < outbndry) && (outer_boundary_flags & INVERT_SET) &&
    localmesh->lastX())) {
      // Use the values in x0 in the boundary

      // x0 is the input
      // bk is the output
      rfft(x0[ix], ncz, &bk(ix, 0));

    } else {
      // b is the input
      // bk is the output
      rfft(b[ix], ncz, &bk(ix, 0));
      //rfft(x0[ix], ncz, &xk(ix, 0));
    }
  }
  ///SCOREP_USER_REGION_END(fftloop);
  ///SCOREP_USER_REGION_DEFINE(kzinit);
  ///SCOREP_USER_REGION_BEGIN(kzinit, "kz init",///SCOREP_USER_REGION_TYPE_COMMON);

  /* Solve differential equation in x for each fourier mode
  * Note that only the non-degenerate fourier modes are being used (i.e. the
  * offset and all the modes up to the Nyquist frequency)
  */
  for (int kz = 0; kz <= maxmode; kz++) {
    for (int ix = 0; ix < ncx; ix++) {
      // Get bk of the current fourier mode
      bcmplx(kz,ix) = bk(ix, kz);
      // Initialize xk1d with cached values (note: in Fourier space)
      xk1d(kz,ix) = x0saved(ix, jy, kz);
      xk1dlast(kz,ix) = x0saved(ix, jy, kz);
    }
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
  //if(first_call(jy,0) || not store_coefficients ){
    for (int kz = 0; kz <= maxmode; kz++) {
      tridagMatrix(&avec(kz,0), &bvec(kz,0), &cvec(kz,0), &bcmplx(kz,0),
      jy,
      // wave number index
      kz,
      // wave number (different from kz only if we are taking a part
      // of the z-domain [and not from 0 to 2*pi])
      kz * kwaveFactor, global_flags, inner_boundary_flags,
      outer_boundary_flags, &A, &C, &D);

      // Patch up internal boundaries
      if(not localmesh->lastX()) { 
	for(int ix = localmesh->xend+1; ix<localmesh->LocalNx ; ix++) {
	  avec(kz,ix) = 0;
	  bvec(kz,ix) = 1;
	  cvec(kz,ix) = 0;
	  bcmplx(kz,ix) = 0;
	}
      } 
      if(not localmesh->firstX()) { 
	for(int ix = 0; ix<localmesh->xstart ; ix++) {
	  avec(kz,ix) = 0;
	  bvec(kz,ix) = 1;
	  cvec(kz,ix) = 0;
	  bcmplx(kz,ix) = 0;
	}
      }
    }
  //}
  ///SCOREP_USER_REGION_END(kzinit);
  ///SCOREP_USER_REGION_DEFINE(initlevels);
  ///SCOREP_USER_REGION_BEGIN(initlevels, "init levels",///SCOREP_USER_REGION_TYPE_COMMON);


  // Initialize levels. Note that the finest grid (level 0) has a different
  // routine to coarse grids (which generally depend on the grid one step
  // finer than itself).
  //
  // If the operator to invert doesn't change from one timestep to another,
  // much of the information for each level may be stored. In particular,
  // if using conventional multigrid (algorithm=0), the coefficients a/b/cvec
  // are stored. Data that cannot be cached (e.g. the changing right-hand
  // sides) is calculated in init_rhs below.
  if(first_call(jy,0) || not store_coefficients ){

    init(levels[0], ncx, jy, avec, bvec, cvec, xs, xe);

    int ncx_coarse = ncx; //(xe-xs+1)/2 + xs + ncx - xe - 1;
    if(max_level>0){
      for(int l = 1; l<=max_level; l++){
	ncx_coarse = (ncx_coarse-4)/2+4;
	init(levels[l], levels[l-1], ncx_coarse, xs, ncx_coarse-3,l,jy); //FIXME assumes mgy=2
      }
    }
  }

  // Compute coefficients that depend on the right-hand side and which
  // therefore change every time.
  init_rhs(levels[0], jy, bcmplx);

  ///SCOREP_USER_REGION_END(initlevels);

  ///SCOREP_USER_REGION_DEFINE(setsoln);
  ///SCOREP_USER_REGION_BEGIN(setsoln, "set level 0 solution",///SCOREP_USER_REGION_TYPE_COMMON);
  for(int kz=0; kz<nmode; kz++){
    for(int ix=0; ix<ncx; ix++){
      levels[0].rvec(kz,ix) = bcmplx(kz,ix);
      levels[0].soln(kz,ix) = xk1d(kz,ix);
      levels[0].solnlast(kz,ix) = xk1d(kz,ix);
    }
    if(algorithm!=0){
      levels[0].xloc(0,kz) = xk1d(kz,xs-1);
      levels[0].xloc(1,kz) = xk1d(kz,xs);
      levels[0].xloc(2,kz) = xk1d(kz,xe);
      levels[0].xloc(3,kz) = xk1d(kz,xe+1);
      levels[0].xloclast(0,kz) = xk1d(kz,xs-1);
      levels[0].xloclast(1,kz) = xk1d(kz,xs);
      levels[0].xloclast(2,kz) = xk1d(kz,xe);
      levels[0].xloclast(3,kz) = xk1d(kz,xe+1);
    }
  }
  ///SCOREP_USER_REGION_END(setsoln);
  ///SCOREP_USER_REGION_DEFINE(initwhileloop);
  ///SCOREP_USER_REGION_BEGIN(initwhileloop, "init while loop",///SCOREP_USER_REGION_TYPE_COMMON);

  int count = 0;
  int subcount = 0;
  int current_level = 0;
  bool down = true;
  auto converged = Array<bool>(nmode);

  auto total = Array<BoutReal>(nmode);
  auto totalold = Array<BoutReal>(nmode);
  auto error_rel = Array<BoutReal>(nmode);
  for(int kz=0; kz<nmode; kz++){
    converged[kz] = false;
    total[kz] = 1e20;
    totalold[kz] = 1e20;
    error_rel[kz] = 0.0;
  }
  int ml;

  ///SCOREP_USER_REGION_END(initwhileloop);
  ///SCOREP_USER_REGION_DEFINE(whileloop);
  ///SCOREP_USER_REGION_BEGIN(whileloop, "while loop",///SCOREP_USER_REGION_TYPE_COMMON);
  /*
  output<<"BEFORE LOOP"<<endl;
  output<<"xloc ";
  for(int i=0;i<4;i++){
    output<<levels[0].xloc(i,1)<<" ";
  }
  output<<endl;
  */

  while(true){

    output<<"START OF LOOP "<<count<<" level "<<current_level<<endl;
    calculate_residual(levels[current_level],converged,jy);

    gauss_seidel_red_black(levels[current_level],converged,jy);
    // For debugging:
    //if(current_level==0){
    //}

    ///SCOREP_USER_REGION_DEFINE(l0rescalc);
    ///SCOREP_USER_REGION_BEGIN(l0rescalc, "level 0 residual calculation",///SCOREP_USER_REGION_TYPE_COMMON);
    if(current_level==0 and subcount==max_cycle-1){
      // Not necessay, but for diagnostics
      for(int kz=0; kz<nmode; kz++){
        if(!converged[kz]){
          totalold[kz] = total[kz];
        }
      }
      calculate_residual(levels[current_level],converged,jy);
      calculate_total_residual(total,error_rel,converged,levels[current_level]);
    }

    ///SCOREP_USER_REGION_END(l0rescalc);
    ///SCOREP_USER_REGION_DEFINE(increment);
    ///SCOREP_USER_REGION_BEGIN(increment, "increment counters",///SCOREP_USER_REGION_TYPE_COMMON);
    ++count;
    ++subcount;
    ///SCOREP_USER_REGION_END(increment);

    // Force at least max_cycle iterations at each level
    // Do not skip with tolerence to minimize comms
    if(subcount < max_cycle){
    }
    else if( all(converged) and current_level==0 ){
    /*
      ml = maxloc(total);
      output<<"Exit "<<count<<" "<<ml<<" "<<total[ml]<<" "<<total[ml]/totalold[ml]<<endl;
    output<<"xloc final"<<endl;
    for(int ix=0; ix<4;ix++){
      output<<" "<<levels[current_level].xloc(ix,kzp).real() << " ";
    }
    output<<endl;
    */
      break;
    }
    else if( not down ){
      output<<"up step"<<endl;
      //calculate_residual(levels[current_level],converged,jy);
      refine(levels[current_level],fine_error,converged);
      current_level--;
      update_solution(levels[current_level],fine_error,converged);

      /*
      if(algorithm!=0 and current_level==0){
	output<<"Extended domain loop"<<endl;
	// If extended domain method, need to set xloc(0) and xloc(3) to be the far
	// interior point of neighbouring domain, which is not local. Therefore get
	// this now via comms.
	// Note this is a repeat of code in Jacobi function, and could be encapsulated.
	if(new_method){
	  struct Message { dcomplex value; bool done; };
	  Array<Message> message_send, message_recv;
	  message_send = Array<Message>(nmode);
	  message_recv = Array<Message>(nmode);
	  MPI_Comm comm = BoutComm::get();
	  int err;
          // Communication receive from left
          if(!localmesh->firstX()){
            for (int kz = 0; kz <= maxmode; kz++) {
              if(!converged[kz]){
	        message_send[kz].value = levels[0].xloc(index_in,kz);
                ///	  message_send[kz].done  = self_in[kz];
              }
            }
            err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 1, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 0, comm, MPI_STATUS_IGNORE);
            for (int kz = 0; kz <= maxmode; kz++) {
              if(!converged[kz]){
	        levels[0].xloc(0,kz) = message_recv[kz].value;
///	        neighbour_in[kz] = message_recv[kz].done;
              }
            }
          }

          // Communicate receive from right
          if(!localmesh->lastX()){
            for (int kz = 0; kz <= maxmode; kz++) {
              if(!converged[kz]){
	        message_send[kz].value = levels[0].xloc(index_out,kz);
///             message_send[kz].done  = self_out[kz];
              }
            }
            err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 0, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 1, comm, MPI_STATUS_IGNORE);
            for (int kz = 0; kz < nmode; kz++) {
              if(!converged[kz]){
	        levels[0].xloc(3,kz) = message_recv[kz].value;
///	        neighbour_out[kz] = message_recv[kz].done;
              }
            }
          }
	}
      }
      */
      subcount=0;

      if(current_level==0){
        down = true;
      }
    }
    else if( down && max_level > 0 ){

      output<<"down step"<<endl;
      calculate_residual(levels[current_level],converged,jy);
      current_level++;
      coarsen(levels[current_level],levels[current_level-1].residual,converged);
      output<<"AFTER COARSEN"<<endl;
      subcount=0;

      // If we are on the coarsest grid, stop trying to coarsen further
      if(current_level==max_level){
        down = false;
      }
    }
    else{
      // When only using one level, need to ensure subcount < max_count
      subcount=0;
    }

    // Throw error if we are performing too many iterations
    // TODO Throw error message appropriate for the method
    if (count>maxits) {
      throw BoutException("LaplaceParallelTriMGNew error: Not converged within maxits=%i iterations. The iteration matrix is not diagonally dominant on processor %i, so there is no guarantee this method will converge. Consider increasing maxits or using a different solver.",maxits,BoutComm::rank());
      //break;
      /*
      // Maximum number of allowed iterations reached.
      // If the iteration matrix is diagonally-dominant, then convergence is
      // guaranteed, so maxits is set too low.
      // Otherwise, the method may or may not converge.
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!is_diagonally_dominant(al(jy,kz),au(jy,kz),bl(jy,kz),bu(jy,kz),jy,kz)){
	  throw BoutException("LaplaceParallelTriMGNew error: Not converged within maxits=%i iterations. The iteration matrix is not diagonally dominant on processor %i, so there is no guarantee this method will converge. Consider increasing maxits or using a different solver.",maxits,BoutComm::rank());
	}
      }
      throw BoutException("LaplaceParallelTriMGNew error: Not converged within maxits=%i iterations. The iteration matrix is diagonally dominant on processor %i and convergence is guaranteed (if all processors are diagonally dominant). Please increase maxits and retry.",maxits,BoutComm::rank());
      //output<<alold(jy,kz)<<" "<<blold(jy,kz)<<" "<<auold(jy,kz)<<" "<<buold(jy,kz)<<endl;
      //output<<al(jy,kz)<<" "<<bl(jy,kz)<<" "<<au(jy,kz)<<" "<<bu(jy,kz)<<endl;
      //output<<Ad<<" "<<Bd<<" "<<Au<<" "<<Bu<<endl;
      */
    }
  }
  ///SCOREP_USER_REGION_END(whileloop);
  ///SCOREP_USER_REGION_DEFINE(afterloop);
  ///SCOREP_USER_REGION_BEGIN(afterloop, "after faff",///SCOREP_USER_REGION_TYPE_COMMON);
  //
  /*
  output<<"AFTER LOOP"<<endl;
  output<<"xloc ";
  for(int i=0;i<4;i++){
    output<<levels[0].xloc(i,1)<<" ";
  }
  output<<endl;
  */

  reconstruct_full_solution(levels[0],jy);

  ++ncalls;
  ipt_mean_its = (ipt_mean_its * BoutReal(ncalls-1)
  + BoutReal(count))/BoutReal(ncalls);

  for(int kz=0; kz<nmode; kz++){
    for(int i=0; i<ncx; i++){
      xk1d(kz,i) = levels[0].soln(kz,i);
      x0saved(i,jy,kz) = levels[0].soln(kz,i);
    }
  }

  for (int kz = 0; kz <= maxmode; kz++) {
    // If the global flag is set to INVERT_KX_ZERO
    if ((global_flags & INVERT_KX_ZERO) && (kz == 0)) {
      dcomplex offset(0.0);
      for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
	offset += xk1d(kz,ix);
      }
      offset /= static_cast<BoutReal>(localmesh->xend - localmesh->xstart + 1);
      for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
	xk1d(kz,ix) -= offset;
      }
    }

    // Store the solution xk for the current fourier mode in a 2D array
    for (int ix = 0; ix < ncx; ix++) {
      xk(ix, kz) = xk1d(kz,ix);
    }
    first_call(jy,kz) = false;
  }
  ///SCOREP_USER_REGION_END(afterloop);

  ///SCOREP_USER_REGION_DEFINE(fftback);
  ///SCOREP_USER_REGION_BEGIN(fftback, "fft back",///SCOREP_USER_REGION_TYPE_COMMON);
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
  ///SCOREP_USER_REGION_END(fftback);
  return x; // Result of the inversion
}

/*
 * Perform an iteration on the reduced grid
 */
void LaplaceParallelTriMGNew::jacobi(Level &l, const int jy, const Array<bool> &converged){

  SCOREP0();
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      
      // Set xloclast to xloc 
      l.xloclast(0,kz) = l.xloc(0,kz);
      l.xloclast(1,kz) = l.xloc(1,kz);
      l.xloclast(2,kz) = l.xloc(2,kz);
      l.xloclast(3,kz) = l.xloc(3,kz);

      // Update my interior points
      l.xloc(1,kz) = l.rl[kz] + l.al(jy,kz)*l.xloclast(0,kz) + l.bl(jy,kz)*l.xloclast(3,kz);
      l.xloc(2,kz) = l.ru[kz] + l.au(jy,kz)*l.xloclast(0,kz) + l.bu(jy,kz)*l.xloclast(3,kz);

      // If I am a boundary processor, apply boundary conditions.
      if( localmesh->firstX() ){
	//l.xloc(0,kz) = ( l.minvb(kz,l.xs-1) + l.upperGuardVector(l.xs-1,jy,kz)*l.xloclast(3,kz) ) / (1.0 - l.lowerGuardVector(l.xs-1,jy,kz));
        l.xloc(0,kz) = -l.cvec(jy,kz,l.xs-1)*l.xloc(1,kz) / l.bvec(jy,kz,l.xs-1);
      }
      if( localmesh->lastX() ){
	//l.xloc(3,kz) = ( l.minvb(kz,l.xe+1) + l.lowerGuardVector(l.xe+1,jy,kz)*l.xloclast(0,kz) ) / (1.0 - l.upperGuardVector(l.xe+1,jy,kz));
        l.xloc(3,kz) = -l.avec(jy,kz,l.xe+1)*l.xloc(2,kz) / l.bvec(jy,kz,l.xe+1);
      }
    }
  }

  ///SCOREP_USER_REGION_DEFINE(comms);
  ///SCOREP_USER_REGION_BEGIN(comms, "communication",SCOREP_USER_REGION_TYPE_COMMON);

  // Communication
  auto sendvec = Array<dcomplex>(nmode);
  auto recvvec = Array<dcomplex>(nmode);
  MPI_Comm comm = BoutComm::get();
  int err;

  // Communication receive from left
  if(!localmesh->firstX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	sendvec[kz] = l.xloc(index_in,kz);
      }
    }
    err = MPI_Sendrecv(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, &recvvec[0], nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, comm, MPI_STATUS_IGNORE);
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	l.xloc(0,kz) = recvvec[kz];
      }
    }
  }

  // Communicate receive from right
  if(!localmesh->lastX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	sendvec[kz] = l.xloc(index_out,kz);
      }
    }
    err = MPI_Sendrecv(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, &recvvec[0], nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, comm, MPI_STATUS_IGNORE);
    for (int kz = 0; kz < nmode; kz++) {
      if(!converged[kz]){
	l.xloc(3,kz) = recvvec[kz];
      }
    }
  }
  ///SCOREP_USER_REGION_END(comms);
}  

/*
 * Perform Gauss-Seidel with red-black colouring on the reduced system
 */
void LaplaceParallelTriMGNew::gauss_seidel_red_black(Level &l, const Array<bool> &converged, const int jy){

  //output<<"GAUSS SEIDEL BLACK RED"<<endl;
  SCOREP0();

  if(l.included){
    //output<<"running GS"<<endl;

    Array<dcomplex> sendvec, recvecin, recvecout;
    sendvec = Array<dcomplex>(nmode);
    recvecin = Array<dcomplex>(nmode);
    recvecout = Array<dcomplex>(nmode);
    MPI_Request rreqin, rreqout, sreqin, sreqout;

    output<<"before"<<endl;
      for(int i=0;i<4;i++){
	output<<l.xloc(i,1)<<" ";
      }
      output<<endl;

    // BLACK SWEEP
    // 
    // Red processors communication only
    if(l.red){
      //output<<"red comm in black sweep"<<endl;
      // Post receives
      if(not localmesh->firstX()){
	//output<<"posting recv from "<<l.proc_in<<endl;
	MPI_Irecv(&recvecin[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 0, BoutComm::get(), &rreqin);
      }
      if(not localmesh->lastX()){ // this is always be true is we force an even core count
	//output<<"posting recv from "<<l.proc_out<<endl;
	MPI_Irecv(&recvecout[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 1, BoutComm::get(), &rreqout);
      }

      // Receive and put data in arrays
      if(!localmesh->firstX()){
	MPI_Wait(&rreqin,MPI_STATUS_IGNORE);
	for(int kz=0; kz < nmode; kz++){
	  if(!converged[kz]){
	    l.xloc(0,kz) = recvecin[kz];
	  }
	}
      }
      if(!localmesh->lastX()){
	MPI_Wait(&rreqout,MPI_STATUS_IGNORE);
	for(int kz=0; kz < nmode; kz++){
	  if(!converged[kz]){
	    l.xloc(3,kz) = recvecout[kz];
	  }
	}
      }
    }

    // Black processors: work and communication
    if(l.black){
      //output<<"black work and comm"<<endl;
      // Black processors do work
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!converged[kz]){ // TODO worthwhile?
	  if(not localmesh->lastX()){
	    //l.xloc(1,kz) = l.rl[kz] + l.al(jy,kz)*l.xloc(0,kz) + l.bl(jy,kz)*l.xloc(3,kz);
	    l.xloc(1,kz) = ( l.rr(1,kz) - l.ar(jy,1,kz)*l.xloc(0,kz) - l.cr(jy,1,kz)*l.xloc(3,kz) ) / l.br(jy,1,kz);
	  }
	  else{
	    // Due to extra point on final proc, indexing of last term is 2, not 3
	    // TODO Index is constant, so could remove this branching by defining index_end
	    //l.xloc(1,kz) = l.rl[kz] + l.al(jy,kz)*l.xloc(0,kz) + l.bl(jy,kz)*l.xloc(2,kz);
	    l.xloc(1,kz) = ( l.rr(1,kz) - l.ar(jy,1,kz)*l.xloc(0,kz) - l.cr(jy,1,kz)*l.xloc(2,kz) ) / l.br(jy,1,kz);
	  }
	}
      }
      // Send same data up and down
      if(not localmesh->firstX()){
	//output<<"sending to  "<<l.proc_in<<endl;
	//MPI_Isend(&l.xloc(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, BoutComm::get(), &sreqin);
	MPI_Send(&l.xloc(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, BoutComm::get());
      }
      if(not localmesh->lastX()){
	//output<<"sending to  "<<l.proc_out<<endl;
	//MPI_Isend(&l.xloc(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, BoutComm::get(), &sreqout);
	MPI_Send(&l.xloc(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, BoutComm::get());
      }
    }

    /*
    output<<"after black"<<endl;
      for(int i=0;i<4;i++){
	output<<l.xloc(i,1)<<" ";
      }
    output<<endl;
    */

    // RED SWEEP
    //
    // Black processors only comms
    if(l.black or localmesh->lastX()){
      // Post receives
      if(not localmesh->firstX()){
	//output<<"posting recv from  "<<l.proc_in<<endl;
	MPI_Irecv(&recvecin[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 0, BoutComm::get(), &rreqin);
      }
      if(not localmesh->lastX()){ // this is always be true is we force an even core count
	//output<<"posting recv from  "<<l.proc_out<<endl;
	MPI_Irecv(&recvecout[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 1, BoutComm::get(), &rreqout);
      }

      // Receive and put data in arrays
      if(!localmesh->firstX()){
	MPI_Wait(&rreqin,MPI_STATUS_IGNORE);
	for(int kz=0; kz < nmode; kz++){
	  if(!converged[kz]){
	    l.xloc(0,kz) = recvecin[kz];
	  }
	}
      }
      if(!localmesh->lastX()){
	MPI_Wait(&rreqout,MPI_STATUS_IGNORE);
	for(int kz=0; kz < nmode; kz++){
	  if(!converged[kz]){
	    l.xloc(3,kz) = recvecout[kz];
	  }
	}
      }
    }

    // Red processors do work and comms
    if(l.red){
      //output<<"calling red"<<endl;
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!converged[kz]){ // TODO is guarding work still worth it without x loops?
	  //l.xloc(1,kz) = l.rl[kz] + l.al(jy,kz)*l.xloc(0,kz) + l.bl(jy,kz)*l.xloc(3,kz);
	  l.xloc(1,kz) = ( l.rr(1,kz) - l.ar(jy,1,kz)*l.xloc(0,kz) - l.cr(jy,1,kz)*l.xloc(3,kz) ) / l.br(jy,1,kz);
	}
      }
    }
    if(localmesh->lastX()){
      //output<<"calling endpoint red"<<endl;
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!converged[kz]){ // TODO is guarding work still worth it without x loops?
	  //l.xloc(2,kz) = l.ru[kz] + l.au(jy,kz)*l.xloc(1,kz) + l.bu(jy,kz)*l.xloc(3,kz);
	  if(l.current_level==0){
	    l.xloc(2,kz) = ( l.rr(2,kz) - l.ar(jy,2,kz)*l.xloc(1,kz) - l.cr(jy,2,kz)*l.xloc(3,kz) ) / l.br(jy,2,kz);
	  }
	  else{
	    l.xloc(2,kz) = ( l.rr(2,kz) - l.ar(jy,2,kz)*l.xloc(0,kz) - l.cr(jy,2,kz)*l.xloc(3,kz) ) / l.br(jy,2,kz);
	  }
	}
      }
    }

    //if(l.red or (localmesh->lastX() and l.current_level>0)){
    if(not l.black){ // red, or last proc when not on level zero
      // Send same data up and down
      if(not localmesh->firstX() and l.red){ // excludes last proc
	//output<<"sending to  "<<l.proc_in<<endl;
	//MPI_Isend(&l.xloc(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, BoutComm::get(), &sreqin);
	MPI_Send(&l.xloc(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, BoutComm::get());
      }
      else if(not localmesh->firstX()){
	//output<<"sending to  "<<l.proc_in<<endl;
	//MPI_Isend(&l.xloc(2,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, BoutComm::get(), &sreqin);
	MPI_Send(&l.xloc(2,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, BoutComm::get());
      }
      if(not localmesh->lastX()){
	//output<<"sending to  "<<l.proc_out<<endl;
	//MPI_Isend(&l.xloc(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, BoutComm::get(), &sreqout);
	MPI_Send(&l.xloc(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, BoutComm::get());
      }
    }

    /*
    output<<"after red"<<endl;
      for(int i=0;i<4;i++){
	output<<l.xloc(i,1)<<" ";
      }
    output<<endl;
    */

    if(l.current_level==0){
      // Update boundaries to match interior points
      // Do this after communication
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!converged[kz]){
	  if(localmesh->firstX()){
	    l.xloc(0,kz) = - l.cvec(jy,kz,l.xs-1)*l.xloc(1,kz) / l.bvec(jy,kz,l.xs-1);
	  }
	  if(localmesh->lastX()){
	    //TODO this is a bug in the mg branch?
	    //l.xloc(3,kz) = - l.avec(jy,kz,l.xe+1)*l.xloc(1,kz) / l.bvec(jy,kz,l.xe+1);
	    l.xloc(3,kz) = - l.avec(jy,kz,l.xe+1)*l.xloc(2,kz) / l.bvec(jy,kz,l.xe+1);
	  }
	}
      }
    }

    /*
    if(!localmesh->lastX()){
      MPI_Wait(&sreqout,MPI_STATUS_IGNORE);
    }
    if(!localmesh->firstX()){
      MPI_Wait(&sreqin,MPI_STATUS_IGNORE);
    }
    */

    
    output<<"after GS"<<endl;
    output<<"xloc ";
      for(int i=0;i<4;i++){
	output<<l.xloc(i,1)<<" ";
      }
      output<<endl;
  }
  else{
    output<<"skipping GS"<<endl;
  }
   
}  

/*
 * Perform a Gauss--Seidel iteration with red black colouring explicitly on the full system 
 * Note that this assumes that each processor has an even number of points
 */
void LaplaceParallelTriMGNew::gauss_seidel_red_black_full_system_comp_comm_overlap(Level &l, const Array<bool> &converged, const int jy){

  SCOREP0();
  Array<dcomplex> sendvecred, recvecred, sendvecblack, recvecblack;
  sendvecred = Array<dcomplex>(nmode);
  recvecred = Array<dcomplex>(nmode);
  sendvecblack = Array<dcomplex>(nmode);
  recvecblack = Array<dcomplex>(nmode);
  MPI_Request rredreq, sredreq, rblackreq, sblackreq;

  // Communicate: final grid point needs data from proc above
  // Overlap comm/comp by not posting recv until needed in final
  // loop iteration
  if(!localmesh->lastX()){
    MPI_Irecv(&recvecred[0], nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, BoutComm::get(), &rredreq);
  }
  // Communicate: to synchronize, first grid point needs data from proc below
  if(!localmesh->firstX()){
    MPI_Irecv(&recvecblack[0], nmode, MPI_DOUBLE_COMPLEX, proc_in, 2, BoutComm::get(), &rblackreq);
  }

  // Red sweep:

  // Do ix = xs first, so we can start the send
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      int ix = l.xs;
      l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1) ) / l.bvec(jy,kz,ix);
      if(!localmesh->firstX()){
        sendvecred[kz] = l.soln(kz,l.xs);
      }
    }
  }

  // Can send lower guard data now
  if(!localmesh->firstX()){
    MPI_Isend(&sendvecred[0], nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, BoutComm::get(), &sredreq);
  }

  // Red sweep: remaining points
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      for (int ix = l.xs+2; ix < l.xe+1; ix+=2) {
	l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1) ) / l.bvec(jy,kz,ix);
      }
    }
  }

  // Black sweep: even points
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      for (int ix = l.xs+1; ix < l.xe; ix+=2) {
	l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1) ) / l.bvec(jy,kz,ix);
      }
    }
  }

  // Black sweep: wait and calculate final interior point
  bool wait_red = false;
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      if(!localmesh->lastX()){
        if(not wait_red){
          wait_red = true;
	  MPI_Wait(&rredreq,MPI_STATUS_IGNORE);
        }
        l.soln(kz,l.xe+1) = recvecred[kz];
      }
      l.soln(kz,l.xe) = ( l.rvec(kz,l.xe) - l.avec(jy,kz,l.xe)*l.soln(kz,l.xe-1) - l.cvec(jy,kz,l.xe)*l.soln(kz,l.xe+1) ) / l.bvec(jy,kz,l.xe);
    }
  }

  if(!localmesh->lastX()){
    for(int kz=0; kz < nmode; kz++){
      if(!converged[kz]){
	sendvecblack[kz] = l.soln(kz,l.xe);
      }
    }
    MPI_Isend(&sendvecblack[0], nmode, MPI_DOUBLE_COMPLEX, proc_out, 2, BoutComm::get(), &sblackreq);
  }

  bool wait_black = false; // have we posted the wait?
  if(!localmesh->firstX()){
    for(int kz=0; kz < nmode; kz++){
      if(!converged[kz]){
        if(not wait_black){
	  wait_black = true;
	  MPI_Wait(&rblackreq,MPI_STATUS_IGNORE);
	}
	l.soln(kz,l.xs-1) = recvecblack[kz];
      }
    }
  }

  if(l.current_level==0){
    // Update boundaries to match interior points
    // Do this after communication, otherwise this breaks on 1 interior pt per proc
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	if(localmesh->firstX()){
	  for (int ix = l.xs-1; ix > 0; ix--) {
	    l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1) ) / l.bvec(jy,kz,ix);
	  }
	  l.soln(kz,0) = ( l.rvec(kz,0) - l.cvec(jy,kz,0)*l.soln(kz,1) ) / l.bvec(jy,kz,0);
	}
	if(localmesh->lastX()){
	  for (int ix = l.xe; ix < l.ncx-1; ix++) {
	    l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1) ) / l.bvec(jy,kz,ix);
	  }
	  l.soln(kz,l.ncx-1) = ( l.rvec(kz,l.ncx-1) - l.avec(jy,kz,l.ncx-1)*l.soln(kz,l.ncx-2) ) / l.bvec(jy,kz,l.ncx-1);
	}
      }
    }
  }

  if(!localmesh->lastX()){
    MPI_Wait(&sblackreq,MPI_STATUS_IGNORE);
  }
  if(!localmesh->firstX()){
    MPI_Wait(&sredreq,MPI_STATUS_IGNORE);
  }
}  

// Write info about levels to screen
void LaplaceParallelTriMGNew::levels_info(const Level l, const int jy){

  SCOREP0();
  output<<endl;
  output<<"Level "<<l.current_level<<endl;
  output<<"xs = "<<l.xs<<" , xe = "<<l.xe<<", ncx = "<<l.ncx<<endl;

  int kz = 0;
  for(int ix = 0; ix<l.ncx; ix++){
    output<<l.avec(jy,kz,ix)<<" "<<l.bvec(jy,kz,ix)<<" "<<l.cvec(jy,kz,ix)<<endl;
  }

  if(!localmesh->firstX() and l.current_level>0){
    output<<"comm "<<l.acomm[0]<<" "<<l.bcomm[0]<<" "<<l.ccomm[0]<<endl;
  }

  if(algorithm!=0){
    output<<"l coefs "<<l.al(jy,kz)<<" "<<l.bl(jy,kz)<<" "<<l.rl[kz]<<endl;
    output<<"u coefs "<<l.au(jy,kz)<<" "<<l.bu(jy,kz)<<" "<<l.ru[kz]<<endl;

  }
}

// Initialization routine for coarser grids. Initialization depends on the grid
// one step finer, lup.
void LaplaceParallelTriMGNew::init(Level &l, const Level lup, int ncx, const int xs, const int xe, const int current_level, const int jy){

  output<<"INIT"<<endl;
  SCOREP0();
  l.xs = xs;
  l.xe = xe;
  l.ncx = ncx;
  int ny = localmesh->LocalNy;
  l.current_level = current_level;

  Array<dcomplex> sendvec, recvecin, recvecout;
  sendvec = Array<dcomplex>(3*nmode);
  recvecin = Array<dcomplex>(3*nmode);
  recvecout = Array<dcomplex>(3*nmode);

  l.myproc = lup.myproc;
  // Whether this proc is involved in the multigrid calculation
  l.included = ( l.myproc%int((pow(2,current_level))) == 0 ) or localmesh->lastX();
  // Whether this proc is involved in the calculation on the grid one level more refined
  l.included_up = lup.included;

  // TODO Probably no longer a problem:
  if(l.xe-l.xs<1){
    throw BoutException("LaplaceParallelTriMGNew error: Coarse grids must contain at least two points on every processor. Please set max_level smaller than %i, use fewer processors, or increase x resolution.",l.current_level);
    // Note: Grids are initialized from finest to coarsest, so l.current_level
    // is the finest grid that fails.
  }

  if(l.included){

    // Colouring of processor for Gauss-Seidel
    l.red = ((l.myproc/int((pow(2,current_level))))%2 == 0);
    l.black = ((l.myproc/int((pow(2,current_level))))%2 == 1);
    // The last processor is a special case. It is always included because of
    // the final grid point, which is treated explicitly. Otherwise it should
    // not be included in either the red or black work.
    if(localmesh->lastX()){
      l.red = false;
      l.black = false;
    }
    output<<"red "<<l.red<<" black "<<l.black<<endl;

    output<<"current level "<< current_level<<endl;
    output<<"skip step "<<int(pow(2,current_level))<<endl;
    // My neighbouring procs
    l.proc_in = l.myproc - int(pow(2,current_level));
    if(localmesh->lastX()){
      l.proc_in += 1;
    }
    int p = l.myproc + int(pow(2,current_level));
    l.proc_out = (p < nproc-1) ? p : nproc - 1;
    output<<"proc_in "<<l.proc_in<<" proc_out "<<l.proc_out<<endl;

    l.avec = Tensor<dcomplex>(ny,nmode,ncx);
    l.bvec = Tensor<dcomplex>(ny,nmode,ncx);
    l.cvec = Tensor<dcomplex>(ny,nmode,ncx);
    l.rvec = Matrix<dcomplex>(nmode,ncx);
    l.residual = Matrix<dcomplex>(nmode,4);
    l.soln = Matrix<dcomplex>(nmode,ncx);
    l.solnlast = Matrix<dcomplex>(nmode,ncx);
    l.xloc = Matrix<dcomplex>(4,nmode);
    l.xloclast = Matrix<dcomplex>(4,nmode);
    l.rl = Array<dcomplex>(nmode);
    l.ru = Array<dcomplex>(nmode);

    for(int kz=0; kz<nmode; kz++){
      for(int ix=0; ix<ncx; ix++){
       l.avec(jy,kz,ix) = lup.avec(jy,kz,ix); 
       l.bvec(jy,kz,ix) = lup.bvec(jy,kz,ix); 
       l.cvec(jy,kz,ix) = lup.cvec(jy,kz,ix); 
      }
    }

    // Coefficients for the reduced iterations
    l.ar = Tensor<dcomplex>(ny,4,nmode);
    l.br = Tensor<dcomplex>(ny,4,nmode);
    l.cr = Tensor<dcomplex>(ny,4,nmode);
    l.rr = Matrix<dcomplex>(4,nmode);

    output<<"BEFORE LOOP"<<endl;

    for(int kz = 0; kz < nmode; kz++){
      if(localmesh->firstX()){
	l.ar(jy,1,kz) = 0.5*lup.ar(jy,1,kz);
	l.br(jy,1,kz) = 0.5*lup.br(jy,1,kz) + 0.25*lup.cr(jy,1,kz) + 0.25*lup.ar(jy,3,kz) + 0.125*lup.br(jy,3,kz);
	l.cr(jy,1,kz) = 0.25*lup.cr(jy,1,kz) + 0.125*lup.br(jy,3,kz) + 0.25*lup.cr(jy,3,kz);
      }
      else{
	l.ar(jy,1,kz) = 0.25*lup.ar(jy,0,kz) + 0.125*lup.br(jy,0,kz) + 0.25*lup.ar(jy,1,kz);
	l.br(jy,1,kz) = 0.125*lup.br(jy,0,kz) + 0.25*lup.cr(jy,0,kz) + 0.25*lup.ar(jy,1,kz) + 0.5*lup.br(jy,1,kz) + 0.25*lup.cr(jy,1,kz) + 0.25*lup.ar(jy,3,kz) + 0.125*lup.br(jy,3,kz);
	l.cr(jy,1,kz) = 0.25*lup.cr(jy,1,kz) + 0.125*lup.br(jy,3,kz) + 0.25*lup.cr(jy,3,kz);
      }

      // Last proc does calculation on index 2 as well as index 1.
      // If current_level=1, the point to my left on the level above it my
      // index 1. Otherwise, it is my index 0. 
      if(localmesh->lastX()){
	if(l.current_level==1){
	  l.ar(jy,2,kz) = 0.25*lup.ar(jy,1,kz) + 0.125*lup.br(jy,1,kz) + 0.25*lup.ar(jy,2,kz);
	  l.br(jy,2,kz) = 0.125*lup.br(jy,1,kz) + 0.25*lup.cr(jy,1,kz) + 0.25*lup.ar(jy,2,kz) + 0.5*lup.br(jy,2,kz);
	  l.cr(jy,2,kz) = 0.5*lup.cr(jy,2,kz);
	}
	else{
	  l.ar(jy,2,kz) = 0.25*lup.ar(jy,0,kz) + 0.125*lup.br(jy,0,kz) + 0.25*lup.ar(jy,2,kz);
	  l.br(jy,2,kz) = 0.125*lup.br(jy,0,kz) + 0.25*lup.cr(jy,0,kz) + 0.25*lup.ar(jy,2,kz) + 0.5*lup.br(jy,2,kz);
	  l.cr(jy,2,kz) = 0.5*lup.cr(jy,2,kz);
	}
      }

      // Need to communicate my index 1 to this level's neighbours
      // Index 2 if last proc.
      if(not localmesh->lastX()){
	sendvec[kz] = l.ar(jy,1,kz);
	sendvec[kz+nmode] = l.br(jy,1,kz);
	sendvec[kz+2*nmode] = l.cr(jy,1,kz);
      }
      else {
	sendvec[kz] = l.ar(jy,2,kz);
	sendvec[kz+nmode] = l.br(jy,2,kz);
	sendvec[kz+2*nmode] = l.cr(jy,2,kz);
      }
    }
    int kzp=1;

    MPI_Comm comm = BoutComm::get();
    int err;

    // Communicate in
    if(!localmesh->firstX()){
      err = MPI_Sendrecv(&sendvec[0], 3*nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, &recvecin[0], 3*nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 0, comm, MPI_STATUS_IGNORE);
    }

    // Communicate out
    if(!localmesh->lastX()){
      err = MPI_Sendrecv(&sendvec[0], 3*nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, &recvecout[0], 3*nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 1, comm, MPI_STATUS_IGNORE);
    }

    for(int kz=0; kz<nmode; kz++){
      if(not localmesh->firstX()){
	l.ar(jy,0,kz) = recvecin[kz];
	l.br(jy,0,kz) = recvecin[kz+nmode];
	l.cr(jy,0,kz) = recvecin[kz+2*nmode];
      }
      if(not localmesh->lastX()){
	l.ar(jy,3,kz) = recvecout[kz];
	l.br(jy,3,kz) = recvecout[kz+nmode];
	l.cr(jy,3,kz) = recvecout[kz+2*nmode];
      }
    }

    output<<l.ar(jy,0,kzp)<<" "<<l.ar(jy,1,kzp)<<" "<<l.ar(jy,2,kzp)<<" "<<l.ar(jy,3,kzp)<<endl;
    output<<l.br(jy,0,kzp)<<" "<<l.br(jy,1,kzp)<<" "<<l.br(jy,2,kzp)<<" "<<l.br(jy,3,kzp)<<endl;
    output<<l.cr(jy,0,kzp)<<" "<<l.cr(jy,1,kzp)<<" "<<l.cr(jy,2,kzp)<<" "<<l.cr(jy,3,kzp)<<endl;
    //levels_info(l,jy);
    output<<"END INIT"<<endl;
  }
}

// Init routine for finest level
void LaplaceParallelTriMGNew::init(Level &l, const int ncx, const int jy, const Matrix<dcomplex> avec, const Matrix<dcomplex> bvec, const Matrix<dcomplex> cvec, const int xs, const int xe){

  // Basic definitions for conventional multigrid
  SCOREP0();
  l.xs = xs;
  l.xe = xe;
  l.ncx = ncx;
  int ny = localmesh->LocalNy;
  l.current_level = 0;

  // Processor information
  l.myproc = myproc; // unique id
  l.proc_in = myproc-1; // in-neighbour
  l.proc_out = myproc+1; // out-neighbour
  l.included = true; // whether processor is included in this level's calculation
  // Colouring of processor for Gauss-Seidel
  l.red = (l.myproc%2 == 0);
  l.black = (l.myproc%2 == 1);
  output<<"red "<<l.red<<" black "<<l.black<<endl;

  l.avec = Tensor<dcomplex>(ny,nmode,ncx);
  l.bvec = Tensor<dcomplex>(ny,nmode,ncx);
  l.cvec = Tensor<dcomplex>(ny,nmode,ncx);
  l.rvec = Matrix<dcomplex>(nmode,ncx);

  // Coefficients for the reduced iterations
  l.ar = Tensor<dcomplex>(ny,4,nmode);
  l.br = Tensor<dcomplex>(ny,4,nmode);
  l.cr = Tensor<dcomplex>(ny,4,nmode);
  l.rr = Matrix<dcomplex>(4,nmode);

  l.residual = Matrix<dcomplex>(nmode,4);
  l.soln = Matrix<dcomplex>(nmode,ncx);
  l.solnlast = Matrix<dcomplex>(nmode,ncx);

  for(int kz=0; kz<nmode; kz++){
    for(int ix=0; ix<ncx; ix++){
     l.avec(jy,kz,ix) = avec(kz,ix); 
     l.bvec(jy,kz,ix) = bvec(kz,ix); 
     l.cvec(jy,kz,ix) = cvec(kz,ix); 
    }
    for(int ix=0; ix<4; ix++){
      l.residual(kz,ix) = 0.0;
    }
  }
  // end basic definitions

  auto evec = Array<dcomplex>(ncx);
  auto tmp = Array<dcomplex>(ncx);

  Array<dcomplex> sendvec, recvecin, recvecout;
  sendvec = Array<dcomplex>(3*nmode);
  recvecin = Array<dcomplex>(3*nmode);
  recvecout = Array<dcomplex>(3*nmode);

  // Define sizes of local coefficients
  l.xloc = Matrix<dcomplex>(4,nmode);
  l.xloclast = Matrix<dcomplex>(4,nmode);
  l.minvb = Matrix<dcomplex>(nmode,ncx);
  l.upperGuardVector = Tensor<dcomplex>(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz / 2 + 1);
  l.lowerGuardVector = Tensor<dcomplex>(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz / 2 + 1);
  l.al = Matrix<dcomplex>(localmesh->LocalNy, localmesh->LocalNz / 2 + 1);
  l.bl = Matrix<dcomplex>(localmesh->LocalNy, localmesh->LocalNz / 2 + 1);
  l.au = Matrix<dcomplex>(localmesh->LocalNy, localmesh->LocalNz / 2 + 1);
  l.bu = Matrix<dcomplex>(localmesh->LocalNy, localmesh->LocalNz / 2 + 1);

  l.alold = Matrix<dcomplex>(localmesh->LocalNy, localmesh->LocalNz / 2 + 1);
  l.blold = Matrix<dcomplex>(localmesh->LocalNy, localmesh->LocalNz / 2 + 1);
  l.auold = Matrix<dcomplex>(localmesh->LocalNy, localmesh->LocalNz / 2 + 1);
  l.buold = Matrix<dcomplex>(localmesh->LocalNy, localmesh->LocalNz / 2 + 1);

  l.rl = Array<dcomplex>(localmesh->LocalNz / 2 + 1);
  l.ru = Array<dcomplex>(localmesh->LocalNz / 2 + 1);
  l.rlold = Array<dcomplex>(localmesh->LocalNz / 2 + 1);
  l.ruold = Array<dcomplex>(localmesh->LocalNz / 2 + 1);

  // Coefs used to compute rl from domain below
  l.r1 = Matrix<dcomplex>(localmesh->LocalNy, nmode);
  l.r2 = Matrix<dcomplex>(localmesh->LocalNy, nmode);

  for (int kz = 0; kz <= maxmode; kz++) {

    ///SCOREP_USER_REGION_DEFINE(invert);
    ///SCOREP_USER_REGION_BEGIN(invert, "invert local matrices",///SCOREP_USER_REGION_TYPE_COMMON);

    // Invert local matrices to find upper/lower guard vectos.
    // Note Minv*b is calculated in init_rhs.
    //
    // Upper interface
    if(not localmesh->lastX()) {
      // Need the xend-th element
      for(int i=0; i<ncx; i++){
	evec[i] = 0.0;
      }
      evec[l.xe+1] = 1.0;
      tridag(&avec(kz,0), &bvec(kz,0), &cvec(kz,0), std::begin(evec),
      std::begin(tmp), ncx);
      for(int i=0; i<ncx; i++){
	l.upperGuardVector(i,jy,kz) = tmp[i];
      }
    } else {
      for(int i=0; i<ncx; i++){
	l.upperGuardVector(i,jy,kz) = 0.0;
      }
    }

    // Lower interface
    if(not localmesh->firstX()) { 
      for(int i=0; i<ncx; i++){
	evec[i] = 0.0;
      }
      evec[l.xs-1] = 1.0;
      tridag(&avec(kz,0), &bvec(kz,0), &cvec(kz,0), std::begin(evec),
      std::begin(tmp), ncx);
      for(int i=0; i<ncx; i++){
	l.lowerGuardVector(i,jy,kz) = tmp[i];
      }
    } else {
      for(int i=0; i<ncx; i++){
	l.lowerGuardVector(i,jy,kz) = 0.0;
      }
    }

    ///SCOREP_USER_REGION_END(invert);
    ///SCOREP_USER_REGION_DEFINE(coefs);
    ///SCOREP_USER_REGION_BEGIN(coefs, "calculate coefs",///SCOREP_USER_REGION_TYPE_COMMON);

    l.bl(jy,kz) = l.upperGuardVector(l.xs,jy,kz);
    l.al(jy,kz) = l.lowerGuardVector(l.xs,jy,kz);

    l.bu(jy,kz) = l.upperGuardVector(l.xe,jy,kz);
    l.au(jy,kz) = l.lowerGuardVector(l.xe,jy,kz);

    l.alold(jy,kz) = l.al(jy,kz);
    l.auold(jy,kz) = l.au(jy,kz);
    l.blold(jy,kz) = l.bl(jy,kz);
    l.buold(jy,kz) = l.bu(jy,kz);

    // First compute coefficients that depend on the matrix to be inverted
    // and which therefore might be constant throughout a run.

    // Boundary processor values to be overwritten when relevant
    dcomplex Ad, Bd, Au, Bu, Atmp, Btmp;
    Ad = 1.0;
    Bd = 0.0;
    Au = 0.0;
    Bu = 1.0;
    // TODO these should be one way sends, not swaps
    if(not localmesh->firstX()){
      // Send coefficients down
      Atmp = l.al(jy,kz);
      Btmp = 0.0;
      if( std::fabs(l.bu(jy,kz)) > 1e-14 ){
	Btmp = l.bl(jy,kz)/l.bu(jy,kz);
	Atmp -= Btmp*l.au(jy,kz);
      }
      // Send these
      Ad = localmesh->communicateXIn(Atmp);
      Bd = localmesh->communicateXIn(Btmp);
    }
    if(not localmesh->lastX()){
      // Send coefficients up
      Atmp = 0.0;
      Btmp = l.bu(jy,kz);
      if( std::fabs(l.al(jy,kz)) > 1e-14 ){
	Atmp = l.au(jy,kz)/l.al(jy,kz);
	Btmp -= Atmp*l.bl(jy,kz);
      }
      // Send these
      Au = localmesh->communicateXOut(Atmp);
      Bu = localmesh->communicateXOut(Btmp);
    }

    dcomplex Delta;
    Delta = 1.0 / (1.0 - l.al(jy,kz)*Bd);
    l.al(jy,kz) = Delta * l.alold(jy,kz) * Ad;
    l.bl(jy,kz) = Delta * l.blold(jy,kz);

    l.r1(jy,kz) = Delta * l.alold(jy,kz);
    l.r2(jy,kz) = Delta;

    // lastX is a special case having two points on the level 0 grid
    if(localmesh->lastX()){
      // Note that if the first proc is also the last proc, then both alold and
      // auold are zero, and l.au = l.auold is already correct.
      if(not localmesh->lastX()){
	l.au(jy,kz) = l.auold(jy,kz)/l.alold(jy,kz);
	l.bu(jy,kz) = l.buold(jy,kz) - l.au(jy,kz) * l.blold(jy,kz); // NB depends on previous line
      }

      // Use BCs to replace x(xe+1) = -avec(xe+1) x(xe) / bvec(xe+1)
      //  => only bl changes
      l.bl(jy,kz) = - l.avec(jy,kz,l.xe+1)*l.bl(jy,kz) / l.bvec(jy,kz,l.xe+1);
    }

    // Now set coefficients for reduced iterations (shared by all levels)
    l.ar(jy,1,kz) = -l.al(jy,kz);
    l.br(jy,1,kz) = 1.0;
    l.cr(jy,1,kz) = -l.bl(jy,kz);

    l.ar(jy,2,kz) = -l.au(jy,kz);
    l.br(jy,2,kz) = 1.0;
    l.cr(jy,2,kz) = -l.bu(jy,kz);

    sendvec[kz] = l.ar(jy,1,kz);
    sendvec[kz+nmode] = l.br(jy,1,kz);
    sendvec[kz+2*nmode] = l.cr(jy,1,kz);


  ///SCOREP_USER_REGION_END(coefs);
  } // end of kz loop

  // Synchonize reduced coefficients with neighbours
  int kzp=1;
  //output<<l.cr(jy,0,kzp)<<" "<<l.cr(jy,1,kzp)<<" "<<l.cr(jy,2,kzp)<<" "<<l.cr(jy,3,kzp)<<endl;

  MPI_Comm comm = BoutComm::get();
  int err;

  // Communicate in
  if(!localmesh->firstX()){
    err = MPI_Sendrecv(&sendvec[0], 3*nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, &recvecin[0], 3*nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 0, comm, MPI_STATUS_IGNORE);
  }

  // Communicate out
  if(!localmesh->lastX()){
    err = MPI_Sendrecv(&sendvec[0], 3*nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, &recvecout[0], 3*nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 1, comm, MPI_STATUS_IGNORE);
  }

  for(int kz=0; kz<nmode; kz++){
    if(not localmesh->firstX()){
      l.ar(jy,0,kz) = recvecin[kz];
      l.br(jy,0,kz) = recvecin[kz+nmode];
      l.cr(jy,0,kz) = recvecin[kz+2*nmode];
    }
    if(not localmesh->lastX()){
      l.ar(jy,3,kz) = recvecout[kz];
      l.br(jy,3,kz) = recvecout[kz+nmode];
      l.cr(jy,3,kz) = recvecout[kz+2*nmode];
    }
  }

  output<<l.ar(jy,0,kzp)<<" "<<l.ar(jy,1,kzp)<<" "<<l.ar(jy,2,kzp)<<" "<<l.ar(jy,3,kzp)<<endl;
  output<<l.br(jy,0,kzp)<<" "<<l.br(jy,1,kzp)<<" "<<l.br(jy,2,kzp)<<" "<<l.br(jy,3,kzp)<<endl;
  output<<l.cr(jy,0,kzp)<<" "<<l.cr(jy,1,kzp)<<" "<<l.cr(jy,2,kzp)<<" "<<l.cr(jy,3,kzp)<<endl;



  //levels_info(l,jy);
}

// Init routine for finest level information that cannot be cached
void LaplaceParallelTriMGNew::init_rhs(Level &l, const int jy, const Matrix<dcomplex> bcmplx){

  SCOREP0();

  int ncz = localmesh->LocalNz;
  auto rlold = Array<dcomplex>(nmode);
  auto ruold = Array<dcomplex>(nmode);
  auto Rd = Array<dcomplex>(nmode);
  auto Ru = Array<dcomplex>(nmode);
  auto Rsendup = Array<dcomplex>(ncz+2); // 2*nmode?
  auto Rsenddown = Array<dcomplex>(ncz+2);
  auto Rrecvup = Array<dcomplex>(ncz+2);
  auto Rrecvdown = Array<dcomplex>(ncz+2);
  auto evec = Array<dcomplex>(l.ncx);
  auto tmp = Array<dcomplex>(l.ncx);
  int err;

  //output<<"INIT RHS"<<endl;
  for (int kz = 0; kz <= maxmode; kz++) {

    ///SCOREP_USER_REGION_DEFINE(invertforrhs);
    ///SCOREP_USER_REGION_BEGIN(invertforrhs, "invert local matrices for rhs",///SCOREP_USER_REGION_TYPE_COMMON);

    // Invert local matrices
    // Calculate Minv*b
    tridag(&l.avec(jy,kz,0), &l.bvec(jy,kz,0), &l.cvec(jy,kz,0), &bcmplx(kz,0),
    &l.minvb(kz,0), l.ncx);
    // Now minvb is a constant vector throughout the iterations

    ///SCOREP_USER_REGION_END(invertforrhs);
    ///SCOREP_USER_REGION_DEFINE(coefsforrhs);
    ///SCOREP_USER_REGION_BEGIN(coefsforrhs, "calculate coefs for rhs",///SCOREP_USER_REGION_TYPE_COMMON);

    l.rl[kz] = l.minvb(kz,l.xs);
    l.ru[kz] = l.minvb(kz,l.xe);
    l.rlold[kz] = l.rl[kz];
    l.ruold[kz] = l.ru[kz];

    // Boundary processor values to be overwritten when relevant
    Rd[kz] = 0.0;
    Ru[kz] = 0.0;
    if(not localmesh->firstX()){
      // Send coefficients down
      Rsenddown[kz] = l.rl[kz];
      if( std::fabs(l.buold(jy,kz)) > 1e-14 ){
	Rsenddown[kz] -= l.ru[kz]*l.blold(jy,kz)/l.buold(jy,kz);
      }
      Rd[kz] = localmesh->communicateXIn(Rsenddown[kz]);
    }
    if(not localmesh->lastX()){
      // Send coefficients up
      Rsendup[kz] = l.ru[kz];
      if( std::fabs(l.alold(jy,kz)) > 1e-14 ){
	Rsendup[kz] -= l.rl[kz]*l.auold(jy,kz)/l.alold(jy,kz);
      }
      Ru[kz] = localmesh->communicateXOut(Rsendup[kz]);
    }
    ///SCOREP_USER_REGION_END(coefsforrhs);
  } // end of kz loop

  // Communicate vector in kz
  if(not localmesh->firstX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      Rsenddown[kz+maxmode] = l.xloclast(2,kz);
    }
  }
  if(not localmesh->lastX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      Rsendup[kz+maxmode] = l.xloclast(1,kz);
    }
  }

  if(not localmesh->firstX()){
    err = MPI_Sendrecv(&Rsenddown[0], 2*nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, &Rrecvdown[0], 2*nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, BoutComm::get(), MPI_STATUS_IGNORE);
  }
  if(not localmesh->lastX()){
    err = MPI_Sendrecv(&Rsendup[0], 2*nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, &Rrecvup[0], 2*nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, BoutComm::get(), MPI_STATUS_IGNORE);
  }

  if(not localmesh->firstX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      Rd[kz] = Rrecvdown[kz];
      l.xloclast(0,kz) = Rrecvdown[kz+maxmode];
    }
  }
  if(not localmesh->lastX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      Ru[kz] = Rrecvup[kz];
      l.xloclast(3,kz) = Rrecvup[kz+maxmode];
    }
  }

  for (int kz = 0; kz <= maxmode; kz++) {
    l.rl[kz] = l.r1(jy,kz)*Rd[kz] + l.r2(jy,kz)*l.rlold[kz];
    l.rr(1,kz) = l.rl[kz];

    // Special case for multiple points on last proc
    // Note that if the first proc is also the last proc, then both alold and
    // auold are zero, and l.ru = l.ruold is already correct.
    if(localmesh->lastX() and not localmesh->firstX()){
      l.ru[kz] = l.ruold[kz] - l.auold(jy,kz)*rlold[kz]/l.alold(jy,kz);
    }
    l.rr(2,kz) = l.ru[kz];
  }
  //levels_info(l,jy);
}

/*
 * Sum and communicate total residual for the reduced system
 */
void LaplaceParallelTriMGNew::calculate_total_residual(Array<BoutReal> &error_abs, Array<BoutReal> &error_rel, Array<bool> &converged, Level &l){

  SCOREP0();
  // Communication arrays:
  // residual in (0 .. nmode-1)
  // solution in (nmode .. 2*nmode-1)
  auto total = Array<BoutReal>(2*nmode); // global summed residual
  auto subtotal = Array<BoutReal>(2*nmode); // local contribution to residual

  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      total[kz] = 0.0;
      total[kz+nmode] = 0.0;

      // Only xs and xe have nonzero residuals
      subtotal[kz] = pow(l.residual(kz,1).real(),2) + pow(l.residual(kz,1).imag(),2) + pow(l.residual(kz,2).real(),2) + pow(l.residual(kz,2).imag(),2);

      // TODO This approximation will increase iteration count. The alternatives are:
      // + reconstructing the solution and calculating properly
      // + multiply approximation by (interior points/2) - this can be done 
      //   at runtime by changing rtol
      // Strictly this should be all contributions to the solution, but this under-approximation saves work.
      subtotal[kz+nmode] = pow(l.xloc(1,kz).real(),2) + pow(l.xloc(1,kz).imag(),2) + pow(l.xloc(2,kz).real(),2) + pow(l.xloc(2,kz).imag(),2);
    }
  }

  // Communication needed to ensure processors break on same iteration
  MPI_Allreduce(&subtotal[0], &total[0], 2*nmode, MPI_DOUBLE, MPI_SUM, BoutComm::get());

  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      error_abs[kz] = sqrt(total[kz]);
      error_rel[kz] = error_abs[kz]/sqrt(total[kz+nmode]);
      if( error_abs[kz] < atol or error_rel[kz] < rtol ){
	//output<<"HERE kz "<<kz<<endl;
	converged[kz] = true;
	l.xloclast(0,kz) = l.xloc(0,kz);
	l.xloclast(1,kz) = l.xloc(1,kz);
	l.xloclast(2,kz) = l.xloc(2,kz);
	l.xloclast(3,kz) = l.xloc(3,kz);
      }
    }
  }
}

/*
 * Calculate residual on a reduced x grid. By construction, the residual is 
 * zero, except at the points on the reduced grid.
 */
void LaplaceParallelTriMGNew::calculate_residual(Level &l, const Array<bool> &converged,const int jy){

  SCOREP0();
  if(l.included){
    output<<"CALCULATE_RESIDUAL"<<endl;
    for(int kz=0; kz<nmode; kz++){
      if(!converged[kz]){
	if(not localmesh->lastX()){
	  //l.residual(kz,1) = l.rl[kz] + l.al(jy,kz)*l.xloc(0,kz) - l.xloc(1,kz) + l.bl(jy,kz)*l.xloc(3,kz);
	  l.residual(kz,1) = l.rr(1,kz) - l.ar(jy,1,kz)*l.xloc(0,kz) - l.br(jy,1,kz)*l.xloc(1,kz) - l.cr(jy,1,kz)*l.xloc(3,kz)  ;
	  l.residual(kz,2) = 0.0; //hack
	}
	else{
	  //l.residual(kz,1) = l.rl[kz] + l.al(jy,kz)*l.xloc(0,kz) - l.xloc(1,kz) + l.bl(jy,kz)*l.xloc(2,kz);
	  //l.residual(kz,2) = l.ru[kz] + l.au(jy,kz)*l.xloc(1,kz) - l.xloc(2,kz) + l.bu(jy,kz)*l.xloc(3,kz);
	  if(l.current_level==0){
	    l.residual(kz,1) = l.rr(1,kz) - l.ar(jy,1,kz)*l.xloc(0,kz) - l.br(jy,1,kz)*l.xloc(1,kz) - l.cr(jy,1,kz)*l.xloc(2,kz)  ;
	    l.residual(kz,2) = l.rr(2,kz) - l.ar(jy,2,kz)*l.xloc(1,kz) - l.br(jy,2,kz)*l.xloc(2,kz) - l.cr(jy,2,kz)*l.xloc(3,kz)  ;
	  }
	  else{
	    l.residual(kz,2) = l.rr(2,kz) - l.ar(jy,2,kz)*l.xloc(0,kz) - l.br(jy,2,kz)*l.xloc(2,kz) - l.cr(jy,2,kz)*l.xloc(3,kz)  ;
	  }
	}
      }
    }

    // Communication
    auto sendvec = Array<dcomplex>(nmode);
    auto recvvec = Array<dcomplex>(nmode);
    MPI_Comm comm = BoutComm::get();
    int err;

    // Communicate in
    if(!localmesh->firstX()){
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!converged[kz]){
	  sendvec[kz] = l.residual(kz,1);
	}
      }
      err = MPI_Sendrecv(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, &recvvec[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 0, comm, MPI_STATUS_IGNORE);
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!converged[kz]){
	  l.residual(kz,0) = recvvec[kz];
	}
      }
    }

    // Communicate out
    if(!localmesh->lastX()){
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!converged[kz]){
	  sendvec[kz] = l.residual(kz,1);
	}
      }
      err = MPI_Sendrecv(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, &recvvec[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 1, comm, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < nmode; kz++) {
	if(!converged[kz]){
	  l.residual(kz,3) = recvvec[kz];
	}
      }
    }

    output<<"residual ";
    output<<l.residual(1,0)<<" "<<l.residual(1,1)<<" "<<l.residual(1,2)<<" "<<l.residual(1,3)<<endl;
  }
  else{
    output<<"skipping residual calculation"<<endl;
  }
}

/*
 * Calculate residual on a full x grid
 */
void LaplaceParallelTriMGNew::calculate_residual_full_system(Level &l, const Array<bool> &converged,const int jy){

  SCOREP0();
  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){

      // Only calculate residual in the physical boundaries, not the halos
      if(localmesh->firstX()){
	l.residual(kz,0) = 0.0; //l.rvec(kz,0) - l.bvec(jy,kz,0)*l.soln(kz,0) - l.cvec(jy,kz,0)*l.soln(kz,1);
	for(int ix=1; ix<l.xs; ix++){
	  l.residual(kz,ix) = l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.bvec(jy,kz,ix)*l.soln(kz,ix) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1);
	}
      }

      // Calculate residual in interior points
      for(int ix=l.xs; ix<l.xe+1; ix++){
	l.residual(kz,ix) = l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.bvec(jy,kz,ix)*l.soln(kz,ix) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1);
      }

      // Only calculate residual in the physical boundaries, not the halos
      if(localmesh->lastX()){
	for(int ix=l.xe+1; ix<l.ncx-1; ix++){
	  l.residual(kz,ix) = l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.bvec(jy,kz,ix)*l.soln(kz,ix) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1);
	}
	l.residual(kz,l.ncx-1) = 0.0; //l.rvec(kz,l.ncx-1) - l.avec(jy,kz,l.ncx-1)*l.soln(kz,l.ncx-2) - l.bvec(jy,kz,l.ncx-1)*l.soln(kz,l.ncx-1);
      }
    }
  }

  // Communication
  auto sendvec = Array<dcomplex>(nmode);
  auto recvvec = Array<dcomplex>(nmode);
  MPI_Comm comm = BoutComm::get();
  int err;

  // Communicate in
  if(!localmesh->firstX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	sendvec[kz] = l.residual(kz,l.xs);
      }
    }
    err = MPI_Sendrecv(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, &recvvec[0], nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, comm, MPI_STATUS_IGNORE);
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	l.residual(kz,l.xs-1) = recvvec[kz];
      }
    }
  }

  // Communicate out
  if(!localmesh->lastX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	sendvec[kz] = l.residual(kz,l.xe);
      }
    }
    err = MPI_Sendrecv(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, &recvvec[0], nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, comm, MPI_STATUS_IGNORE);
    for (int kz = 0; kz < nmode; kz++) {
      if(!converged[kz]){
	l.residual(kz,l.xe+1) = recvvec[kz];
      }
    }
  }
}

/*
 * Coarsen the fine residual
 */
void LaplaceParallelTriMGNew::coarsen(Level &l, const Matrix<dcomplex> &fine_residual, const Array<bool> &converged){

  SCOREP0();
  if(l.included){ // whether this processor is included in multigrid?
    output<<"fine residual";
    output<<fine_residual(1,0)<<" "<<fine_residual(1,1)<<" "<<fine_residual(1,2)<<" "<<fine_residual(1,3)<<endl;
    for(int kz=0; kz<nmode; kz++){
      if(!converged[kz]){
	if(not localmesh->lastX()){
	  l.residual(kz,1) = 0.25*fine_residual(kz,0) + 0.5*fine_residual(kz,1) + 0.25*fine_residual(kz,3);
	}
	else{
	  // NB point(kz,1) on last proc only used on level=0
	  l.residual(kz,2) = 0.25*fine_residual(kz,1) + 0.5*fine_residual(kz,2) + 0.25*fine_residual(kz,3);
	}
	for(int ix=0; ix<4; ix++){
	  l.xloc(kz,ix) = 0.0;
	  l.xloclast(kz,ix) = 0.0;
	}
	//l.rvec(kz,ix) = l.residual(kz,ix);
	l.rr(1,kz) = l.residual(kz,1);
	if(localmesh->lastX()){
	  l.rr(2,kz) = l.residual(kz,2);
	}
      }
    }
    output<<"residual after coarsening";
    output<<l.residual(1,0)<<" "<<l.residual(1,1)<<" "<<l.residual(1,2)<<" "<<l.residual(1,3)<<endl;
  }
}

/*
 * Update the solution on the refined grid by adding the error calculated on the coarser grid.
 */
void LaplaceParallelTriMGNew::update_solution(Level &l, const Matrix<dcomplex> &fine_error, const Array<bool> &converged){

  Array<dcomplex> sendvec, recvecin, recvecout;
  sendvec = Array<dcomplex>(nmode);
  recvecin = Array<dcomplex>(nmode);
  recvecout = Array<dcomplex>(nmode);
  MPI_Comm comm = BoutComm::get();
  int err;

  SCOREP0();
  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      for(int ix=0; ix<4; ix++){
	l.xloc(ix,kz) += fine_error(ix,kz);
      }
      sendvec[kz] = l.xloc(1,kz);
    }
  }

  // Fine error is only correct locally, so need to synchronize xloc now

  // Communicate in
  if(!localmesh->firstX()){
    output<<"sending "<<sendvec[1]<<" to "<<l.proc_in<<endl;
    err = MPI_Sendrecv(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, &recvecin[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 0, comm, MPI_STATUS_IGNORE);
    output<<"recving "<<recvecin[1]<<" from "<<l.proc_in<<endl;
  }

  // Communicate out
  if(!localmesh->lastX()){
    output<<"sending "<<sendvec[1]<<" to "<<l.proc_out<<endl;
    err = MPI_Sendrecv(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, &recvecout[0], nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 1, comm, MPI_STATUS_IGNORE);
    output<<"recving "<<recvecout[1]<<" from "<<l.proc_out<<endl;
  }

  for(int kz=0; kz<nmode; kz++){
    if(not localmesh->firstX()){
      l.xloc(0,kz) = recvecin[kz];
    }
    if(not localmesh->lastX()){
      l.xloc(3,kz) = recvecout[kz];
    }
  }

}

/*
 * Refine the reduced system.
 * There are three types of proc to cover:
 *  + procs included at this level. Calculate error and send contributions to neighbours
 *  + procs not included at this level but included at the refined level. Receive contributions
 *  + procs included neither on this level or the level above. Do nothing
 */
void LaplaceParallelTriMGNew::refine(Level &l, Matrix<dcomplex> &fine_error, const Array<bool> &converged){

  SCOREP0();
  Array<dcomplex> sendvec, recvecin, recvecout;
  sendvec = Array<dcomplex>(nmode);
  recvecin = Array<dcomplex>(nmode);
  recvecout = Array<dcomplex>(nmode);
  comm_handle recv[2];

  // Need to update this for my actual neighbours!

  if(l.included){
    for(int kz=0; kz<nmode; kz++){
      if(!converged[kz]){
	fine_error(1,kz) = l.xloc(1,kz);
	sendvec[kz] = l.xloc(1,kz);
	if(localmesh->lastX()){
	  fine_error(2,kz) = l.xloc(2,kz);
	}
      }
    }

    if(!localmesh->lastX()){
      localmesh->sendXOut(&sendvec[0],nmode,0);
    }
    if(!localmesh->firstX()){
      localmesh->sendXIn(&sendvec[0],nmode,1);
    }

  }
  else if(l.included_up){
    if(!localmesh->firstX()){
      recv[0] = localmesh->irecvXIn(&recvecin[0], nmode, 0);
    }
    if(!localmesh->lastX()){
      recv[1] = localmesh->irecvXOut(&recvecout[0], nmode, 1);
    }

    if(!localmesh->firstX()){
      localmesh->wait(recv[0]);
      for(int kz=0; kz < nmode; kz++){
	if(!converged[kz]){
	  fine_error(1,kz) = 0.5*recvecin[kz];
	}
      }
    }
    if(!localmesh->lastX()){
      localmesh->wait(recv[1]);
      for(int kz=0; kz < nmode; kz++){
	if(!converged[kz]){
	  fine_error(1,kz) += 0.5*recvecout[kz];
	}
      }
    }
  }
  // Special case where we need to fill (1,kz) on final proc
  if(localmesh->lastX() and l.current_level==1){
    for(int kz=0; kz < nmode; kz++){
      if(!converged[kz]){
	fine_error(1,kz) = 0.5*(l.xloc(0,kz) + l.xloc(2,kz));
      }
    }
  }

  output<<"fine error";
  output<<fine_error(0,1)<<" "<<fine_error(1,1)<<" "<<fine_error(2,1)<<" "<<fine_error(3,1)<<endl;
}

/*
 * Synchronize the values of a reduced field(4,nmode) between processors that
 * are neighbours on level l. This assumes each processor's value of
 * field(1,:) is correct, and puts the in-neighbour's value into field(0,:)
 * and out-neighbour's value into field(3,:).
 * NB We do not skipped converged modes in kz.
 */
void LaplaceParallelTriMGNew::synchronize_reduced_field(const Level &l, Matrix<dcomplex> &field){

  if(l.included){

    //Array<dcomplex> sendvec, recvecin, recvecout;
    //sendvec = Array<dcomplex>(nmode);
    //recvecin = Array<dcomplex>(nmode);
    //recvecout = Array<dcomplex>(nmode);
    //comm_handle recv[2];
    MPI_Comm comm = BoutComm::get();
    int err;

    // Communicate in
    if(!localmesh->firstX()){
      err = MPI_Sendrecv(&field(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 1, &field(3,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 0, comm, MPI_STATUS_IGNORE);
    }

    // Communicate out
    if(!localmesh->lastX()){
      err = MPI_Sendrecv(&field(1,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, &field(0,0), nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 1, comm, MPI_STATUS_IGNORE);
    }

    output<<"field ";
    int kzp=1;
    output<<field(kzp,0)<<" "<<field(kzp,1)<<" "<<field(kzp,2)<<" "<<field(kzp,3)<<endl;
  }

}
