/**************************************************************************
 * Perpendicular Laplacian inversion. Serial code using FFT
 * and tridiagonal solver.
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
#include "parallel_tri_multigrid.hxx"

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

LaplaceParallelTriMG::LaplaceParallelTriMG(Options *opt, CELL_LOC loc, Mesh *mesh_in)
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
void LaplaceParallelTriMG::resetSolver(){
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
void LaplaceParallelTriMG::get_initial_guess(const int jy, const int kz, Matrix<dcomplex> &minvb,
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
bool LaplaceParallelTriMG::is_diagonally_dominant(const dcomplex al, const dcomplex au, const dcomplex bl, const dcomplex bu, const int jy, const int kz) {

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
void LaplaceParallelTriMG::reconstruct_full_solution(Level &l, const int jy, Matrix<dcomplex> &halos){
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
 * Reconstruct the full solution from the subproblem and the halo cell values
 */
void LaplaceParallelTriMG::reconstruct_full_solution(Level &l, const int jy){
  SCOREP0();
  for (int kz = 0; kz < nmode; kz++) {
    l.soln(kz,l.xs-1) = l.xloc(0,kz);
    for(int i=l.xs; i<l.xe+1; i++){
      l.soln(kz,i) = l.minvb(kz,i) + l.upperGuardVector(i,jy,kz)*l.xloclast(3,kz) + l.lowerGuardVector(i,jy,kz)*l.xloclast(0,kz);
    }
    l.soln(kz,l.xe+1) = l.xloc(3,kz);
  }
}

bool LaplaceParallelTriMG::all(const Array<bool> a){
  SCOREP0();
  for(int i=0; i<a.size(); i++){
    if(a[i]==false){
      return false;
    }
  }
  return true;
}

bool LaplaceParallelTriMG::any(const Array<bool> a){
  SCOREP0();
  for(int i=0; i<a.size(); i++){
    if(a[i]==true){
      return true;
    }
  }
  return false;
}

BoutReal LaplaceParallelTriMG::max(const Array<BoutReal> a){
  SCOREP0();
  BoutReal maxval = 0.0;
  for(int i=0; i<a.size(); i++){
    if(a[i]>maxval){
      maxval = a[i];
    }
  }
  return maxval;
}

int LaplaceParallelTriMG::maxloc(const Array<BoutReal> a){
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

FieldPerp LaplaceParallelTriMG::solve(const FieldPerp& b) { return solve(b, b); }

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
//FieldPerp LaplaceParallelTriMG::solve(const FieldPerp& b, const FieldPerp& x0, const FieldPerp& b0) {
FieldPerp LaplaceParallelTriMG::solve(const FieldPerp& b, const FieldPerp& x0) {

  SCOREP0();
  Timer timer("invert"); ///< Start timer

  SCOREP_USER_REGION_DEFINE(initvars);
  SCOREP_USER_REGION_BEGIN(initvars, "init vars",SCOREP_USER_REGION_TYPE_COMMON);

  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  TRACE("LaplaceParallelTriMG::solve(const, const)");

  FieldPerp x{emptyFrom(b)};

  // Info for halo swaps
  int xproc = localmesh->getXProcIndex();
  int yproc = localmesh->getYProcIndex();
  int myproc = yproc * localmesh->getNXPE() + xproc;
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
  auto xloc = Matrix<dcomplex>(4,nmode);
  auto xloclast = Matrix<dcomplex>(4,nmode);

  BoutReal kwaveFactor = 2.0 * PI / coords->zlength();

  // Should we store coefficients?
  store_coefficients = not (inner_boundary_flags & INVERT_AC_GRAD);
  store_coefficients = store_coefficients && not (outer_boundary_flags & INVERT_AC_GRAD);
  store_coefficients = store_coefficients && not (inner_boundary_flags & INVERT_SET);
  store_coefficients = store_coefficients && not (outer_boundary_flags & INVERT_SET);
  //store_coefficients = false;

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
  *        LaplaceParallelTriMG::solve()
  * bk1d = The 1d array of bk
  * xk   = The fourier transformed of x, where x the output of
  *        LaplaceParallelTriMG::solve()
  * xk1d = The 1d array of xk
  */
  auto bk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto bk1d = Array<dcomplex>(ncx);
  auto xk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto xk1d = Matrix<dcomplex>(ncz/2+1,ncx);
  auto xk1dlast = Matrix<dcomplex>(ncz/2+1,ncx);

  // Error interpolated onto a finer grid
  auto fine_error = Matrix<dcomplex>(ncz/2+1,ncx);

  // Define indexing of xloc that depends on method. Doing this now removes
  // branching in tight loops
  index_in = 1;
  index_out = 2;
  if(new_method){
    index_in = 2;
    index_out = 1;
  }

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
  SCOREP_USER_REGION_END(fftloop);
  SCOREP_USER_REGION_DEFINE(kzinit);
  SCOREP_USER_REGION_BEGIN(kzinit, "kz init",SCOREP_USER_REGION_TYPE_COMMON);

  /* Solve differential equation in x for each fourier mode
  * Note that only the non-degenerate fourier modes are being used (i.e. the
  * offset and all the modes up to the Nyquist frequency)
  */
  for (int kz = 0; kz <= maxmode; kz++) {
    // set bk1d
    for (int ix = 0; ix < ncx; ix++) {
      // Get bk of the current fourier mode
      bcmplx(kz,ix) = bk(ix, kz);
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
  SCOREP_USER_REGION_END(kzinit);
  SCOREP_USER_REGION_DEFINE(initlevels);
  SCOREP_USER_REGION_BEGIN(initlevels, "init levels",SCOREP_USER_REGION_TYPE_COMMON);


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
  init_rhs(levels[0], jy, bcmplx);
  SCOREP_USER_REGION_END(initlevels);

  SCOREP_USER_REGION_DEFINE(setsoln);
  SCOREP_USER_REGION_BEGIN(setsoln, "set level 0 solution",SCOREP_USER_REGION_TYPE_COMMON);
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
  SCOREP_USER_REGION_END(setsoln);
  SCOREP_USER_REGION_DEFINE(initwhileloop);
  SCOREP_USER_REGION_BEGIN(initwhileloop, "init while loop",SCOREP_USER_REGION_TYPE_COMMON);

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

  SCOREP_USER_REGION_END(initwhileloop);
  SCOREP_USER_REGION_DEFINE(whileloop);
  SCOREP_USER_REGION_BEGIN(whileloop, "while loop",SCOREP_USER_REGION_TYPE_COMMON);

  while(true){

    output<<"xloc count "<<count<<" "<<current_level<<endl;
    for(int ix=0; ix<4;ix++){
      output<<" "<<levels[current_level].xloc(ix,1).real() << " ";
    }
    output<<endl;
    output<<"xloclast "<<endl;
    for(int ix=0; ix<4;ix++){
      output<<" "<<levels[current_level].xloclast(ix,1).real() << " ";
    }
    output<<endl;

    if(algorithm==0 or current_level!=0){ 
      //jacobi_full_system(levels[current_level],jy);
      //gauss_seidel_red_black_full_system(levels[current_level],converged,jy);
      // This version is ~20% faster
      gauss_seidel_red_black_full_system_comp_comm_overlap(levels[current_level],converged,jy);
    }
    else{
      jacobi(levels[current_level], jy, converged);

      /*
      output<<"xloc count "<<count<<endl;
      for(int ix=0; ix<4;ix++){
	output<<" "<<levels[current_level].xloc(ix,1).real() << " ";
      }
      output<<endl;
      output<<"xloclast "<<endl;
      for(int ix=0; ix<4;ix++){
	output<<" "<<levels[current_level].xloclast(ix,1).real() << " ";
      }
      output<<endl;
      */

      //reconstruct_full_solution(levels[0],jy);
      /*
      output<< "soln "<<count<<" "<<jy<<" "<<levels[current_level].xs<<" "<<levels[current_level].xe<<" "<<levels[current_level].ncx<<" "<<current_level<<endl;
      for(int ix=0; ix<levels[current_level].ncx;ix++){
	output<<" "<<levels[current_level].soln(1,ix).real() << " ";
      }
      output<<endl;
      output<< "solnlast "<<count<<" "<<jy<<" "<<levels[current_level].xs<<" "<<levels[current_level].xe<<" "<<levels[current_level].ncx<<" "<<current_level<<endl;
      for(int ix=0; ix<levels[current_level].ncx;ix++){
	output<<" "<<levels[current_level].solnlast(1,ix).real() << " ";
      }
      output<<endl;
      */
    }

    SCOREP_USER_REGION_DEFINE(l0rescalc);
    SCOREP_USER_REGION_BEGIN(l0rescalc, "level 0 residual calculation",SCOREP_USER_REGION_TYPE_COMMON);
    if(current_level==0 and subcount==max_cycle-1){
      // Not necessay, but for diagnostics
      for(int kz=0; kz<nmode; kz++){
	if(!converged[kz]){
	  totalold[kz] = total[kz];
	}
      }
      if(algorithm==0 or current_level!=0){ 
	calculate_residual_full_system(levels[current_level],converged,jy);
	calculate_total_residual_full_system(total,error_rel,converged,levels[current_level]);

	/*
	  int kz = 1;
	  output<<"Residual soln "<<jy<<" "<<count<<" "<<current_level;

	  output<<" "<<converged[kz];
	  output<<" "<<total[kz]<<endl;
	  for(int ix=0; ix<levels[0].ncx;ix++){
	    if(ix<levels[current_level].ncx){
	      output<<" "<<levels[current_level].residual(kz,ix).real();
	    }
	    else{
	      output<<" "<<0;
	    }
	  }
	  output<<endl;
	  */
      }
      else{
        calculate_residual(levels[current_level],converged,jy);
        calculate_total_residual(total,error_rel,converged,levels[current_level]);

	/*
	output<<"Totals "<<endl;
	for(int kz=0; kz<nmode; kz++){
	  output<< total[kz]<<" ";
	}
	output<<endl;
	*/
	//if(jy==3){ 
	{
	  /*
	  int kz=1;
	  output<<"Residual xloc "<<jy<<" "<<count<<" "<<current_level;

	  //output<<" "<<converged[kz];
	  output<<" "<<total[kz]<<" "<<error_rel[kz]<<" "
	    << max(total)<< " "<< atol <<" "<< max(error_rel)<<" "<< rtol <<endl;
	      //output<<" "<<levels[current_level].residual(kz,ix).real()<<" "<<levels[current_level].residual(kz,ix).imag();
	  output<<" "<<levels[current_level].residual(kz,levels[current_level].xs-1).real();
	  output<<" "<<levels[current_level].residual(kz,levels[current_level].xs).real();
	  output<<" "<<levels[current_level].residual(kz,levels[current_level].xe).real();
	  output<<" "<<levels[current_level].residual(kz,levels[current_level].xe+1).real();
	  output<<endl;
	  */
	  /*
	  ml = maxloc(total);
          output<<"maxloc abs "<<" "<<ml<<" "<<total[ml]<<" "<<total[ml]/totalold[ml]<<" "<<converged[ml]<<" "<<error_rel[ml]<<endl;
	  ml = maxloc(error_rel);
          output<<"maxloc rel "<<" "<<ml<<" "<<error_rel[ml]<<" "<<converged[ml]<<" "<<total[ml]<<endl;
	  */
	}
	/*
	output<<"xloc after correction "<<count<<endl;
	for(int ix=0; ix<4;ix++){
	  output<<" "<<levels[current_level].xloc(ix,1).real() << " ";
	}
	output<<endl;
	output<<"xloclast "<<endl;
	for(int ix=0; ix<4;ix++){
	  output<<" "<<levels[current_level].xloclast(ix,1).real() << " ";
	}
	output<<endl;
	*/
	//reconstruct_full_solution(levels[0],jy);
	/*
	output<< "soln "<<count<<" "<<jy<<" "<<levels[current_level].xs<<" "<<levels[current_level].xe<<" "<<levels[current_level].ncx<<" "<<current_level<<endl;
	for(int ix=0; ix<levels[current_level].ncx;ix++){
	  output<<" "<<levels[current_level].soln(1,ix).real() << " ";
	}
	output<<endl;
	output<< "solnlast "<<count<<" "<<jy<<" "<<levels[current_level].xs<<" "<<levels[current_level].xe<<" "<<levels[current_level].ncx<<" "<<current_level<<endl;
	for(int ix=0; ix<levels[current_level].ncx;ix++){
	  output<<" "<<levels[current_level].solnlast(1,ix).real() << " ";
	}
	output<<endl;
	*/
      }
    }

    SCOREP_USER_REGION_END(l0rescalc);
    SCOREP_USER_REGION_DEFINE(increment);
    SCOREP_USER_REGION_BEGIN(increment, "increment counters",SCOREP_USER_REGION_TYPE_COMMON);
    ++count;
    ++subcount;
    SCOREP_USER_REGION_END(increment);

    SCOREP_USER_REGION_DEFINE(solneqsolnlast);
    SCOREP_USER_REGION_BEGIN(solneqsolnlast, "soln = soln last",SCOREP_USER_REGION_TYPE_COMMON);
    if(algorithm==0 or current_level!=0){ 
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!converged[kz]){
	  for (int ix = 0; ix < levels[current_level].ncx; ix++) {
	    levels[current_level].solnlast(kz,ix) = levels[current_level].soln(kz,ix);
	  }
	  /*
	  levels[current_level].xloclast(0,kz) = levels[current_level].xloc(0,kz);
	  levels[current_level].xloclast(1,kz) = levels[current_level].xloc(1,kz);
	  levels[current_level].xloclast(2,kz) = levels[current_level].xloc(2,kz);
	  levels[current_level].xloclast(3,kz) = levels[current_level].xloc(3,kz);
	  */
	}
      }
    }
    SCOREP_USER_REGION_END(solneqsolnlast);
    //ml = maxloc(total);
    //output<<""<<jy<<" "<<count<<" "<<ml<<" "<<total[ml]<<" "<<total[ml]/totalold[ml]<<endl;
    //output<<current_level<<" "<<jy<<" "<<count<<" "<<0<<" "<<total[0]<<" "<<total[0]/totalold[0]<<" "<<error_rel[0]<<endl;

    // Force at least max_cycle iterations at each level
    // Do not skip with tolerence to minimize comms
    //if(subcount < max_cycle and total>rtol){
    if(subcount < max_cycle){
    }
    else if( all(converged) and current_level==0 ){
    //else if( max(error_rel) < rtol and current_level==0 ){
      //ml = maxloc(total);
      //output<<"Exit "<<jy<<" "<<count<<" "<<ml<<" "<<total[ml]<<" "<<total[ml]/totalold[ml]<<endl;
      /*
      output<<jy<<" "<<count<<" ";
      for(int kz=0;kz<nmode;kz++){
	output<<total[kz];
      }
      output<<endl;
      */
      //
      /*
      {
      int kz=0;
      output<<count<<" "<<current_level;

      output<<" "<<total;
      for(int ix=0; ix<levels[0].ncx;ix++){
	if(ix<levels[current_level].ncx){
	  output<<" "<<levels[current_level].residual(kz,ix).real();
	}
	else{
	  output<<" "<<0;
	}
      }
      output<<endl;
      }
      */
      break;
    }
    //else if( total < rtol or current_level==max_level){
    else if( not down ){
      calculate_residual_full_system(levels[current_level],converged,jy);
      //refine(xloc,xloclast);
      refine_full_system(levels[current_level],fine_error,converged);
      current_level--;
      update_solution(levels[current_level],fine_error,converged);
      if(algorithm!=0 and current_level==0){
	for(int kz=0; kz<nmode; kz++){
	  if(not converged[kz]){
	    levels[0].xloc(0,kz) = levels[0].soln(kz,levels[0].xs-1);
	    levels[0].xloc(1,kz) = levels[0].soln(kz,levels[0].xs);
	    levels[0].xloc(2,kz) = levels[0].soln(kz,levels[0].xe);
	    levels[0].xloc(3,kz) = levels[0].soln(kz,levels[0].xe+1);
	  }
	}
      }
      subcount=0;

      if(current_level==0){
	down = true;
      }
    }
    //else if( current_level < max_level ){
    else if( down && max_level > 0 ){

      // Coarsening requires data from the grid BEFORE it is made coarser
      //coarsen(levels[current_level],xloc,xloclast,jy);

      // Calculate residual on finer grid
      if(algorithm!=0 and current_level==0){ 
	reconstruct_full_solution(levels[0],jy);
      }
      calculate_residual_full_system(levels[current_level],converged,jy);
      current_level++;
      coarsen_full_system(levels[current_level],levels[current_level-1].residual,converged);
      subcount=0;

      if(current_level==max_level){
	down = false;
      }
    //{
    //int kz=0;
    /*
    for(int ix=0; ix<levels[0].ncx;ix++){
      if(ix<levels[current_level].ncx){
	output<<" "<<levels[current_level].residual(kz,ix);
      }
      else{
	output<<" "<<0;
      }
    }
    output<<endl;
    */
    //}
    }
    else{
      // When only using one level, need to ensure subcount < max_count
      subcount=0;
    }

    // Implicitly, else = carry on doing iterations at this level
    if (count>maxits) {
      break;
      /*
      // Maximum number of allowed iterations reached.
      // If the iteration matrix is diagonally-dominant, then convergence is
      // guaranteed, so maxits is set too low.
      // Otherwise, the method may or may not converge.
      for (int kz = 0; kz <= maxmode; kz++) {
	if(!is_diagonally_dominant(al(jy,kz),au(jy,kz),bl(jy,kz),bu(jy,kz),jy,kz)){
	  throw BoutException("LaplaceParallelTriMG error: Not converged within maxits=%i iterations. The iteration matrix is not diagonally dominant on processor %i, so there is no guarantee this method will converge. Consider increasing maxits or using a different solver.",maxits,BoutComm::rank());
	}
      }
      throw BoutException("LaplaceParallelTriMG error: Not converged within maxits=%i iterations. The iteration matrix is diagonally dominant on processor %i and convergence is guaranteed (if all processors are diagonally dominant). Please increase maxits and retry.",maxits,BoutComm::rank());
      //output<<alold(jy,kz)<<" "<<blold(jy,kz)<<" "<<auold(jy,kz)<<" "<<buold(jy,kz)<<endl;
      //output<<al(jy,kz)<<" "<<bl(jy,kz)<<" "<<au(jy,kz)<<" "<<bu(jy,kz)<<endl;
      //output<<Ad<<" "<<Bd<<" "<<Au<<" "<<Bu<<endl;
      */
    }
  }
  SCOREP_USER_REGION_END(whileloop);
  SCOREP_USER_REGION_DEFINE(afterloop);
  SCOREP_USER_REGION_BEGIN(afterloop, "after faff",SCOREP_USER_REGION_TYPE_COMMON);

  if(algorithm!=0){
    reconstruct_full_solution(levels[0],jy);
    /*
    output<< "soln "<<count<<" "<<jy<<" "<<levels[current_level].xs<<" "<<levels[current_level].xe<<" "<<levels[current_level].ncx<<" "<<current_level<<endl;
    for(int ix=0; ix<levels[current_level].ncx;ix++){
      output<<" "<<levels[current_level].soln(2,ix).real() << " "<< levels[current_level].soln(2,ix).imag() <<endl;
    }
    */
  /*
  else{
    output<< "xloc "<<count<<" "<<jy<<" "<<levels[current_level].xs<<" "<<levels[current_level].xe<<" "<<levels[current_level].ncx<<current_level<<endl;
    for(int ix=0; ix<4;ix++){
      output<<" "<<levels[current_level].xloclast(ix,2).real() << " "<< levels[current_level].xloclast(ix,2).imag() <<endl;
    }
  }
  output<< endl;
  */
    for (int kz = 0; kz <= maxmode; kz++) {
      xloclast(0,kz) = levels[0].soln(kz,levels[0].xs-1);
      xloclast(1,kz) = levels[0].soln(kz,levels[0].xs);
      xloclast(2,kz) = levels[0].soln(kz,levels[0].xe);
      xloclast(3,kz) = levels[0].soln(kz,levels[0].xe+1);
      xloc(0,kz) = levels[0].soln(kz,levels[0].xs-1);
      xloc(1,kz) = levels[0].soln(kz,levels[0].xs);
      xloc(2,kz) = levels[0].soln(kz,levels[0].xe);
      xloc(3,kz) = levels[0].soln(kz,levels[0].xe+1);
    }
  }

  ++ncalls;
  ipt_mean_its = (ipt_mean_its * BoutReal(ncalls-1)
  + BoutReal(count))/BoutReal(ncalls);

  if(algorithm!=0){
    for (int kz = 0; kz <= maxmode; kz++) {
      // Original method:
      xk1d(kz,xs-1) = xloc(0,kz);
      xk1d(kz,xs)   = xloc(1,kz);
      xk1dlast(kz,xs-1) = xloclast(0,kz);
      xk1dlast(kz,xs)   = xloclast(1,kz);
      xk1d(kz,xe)   = xloc(2,kz);
      xk1d(kz,xe+1) = xloc(3,kz);
      xk1dlast(kz,xe)   = xloclast(2,kz);
      xk1dlast(kz,xe+1) = xloclast(3,kz);

      if(new_method){
	SCOREP_USER_REGION_DEFINE(getnewmethodsoln);
	SCOREP_USER_REGION_BEGIN(getnewmethodsoln, "get new method soln",SCOREP_USER_REGION_TYPE_COMMON);
	dcomplex d = 1.0/(levels[0].buold(jy,kz)*levels[0].alold(jy,kz) - levels[0].blold(jy,kz)*levels[0].auold(jy,kz));
	// If boundary processor, halo cell is already correct, and d is undefined.
	// Lower boundary proc => al = au = 0
	// Upper boundary proc => bl = bu = 0
	if(not localmesh->firstX() and not localmesh->lastX()){
	  // General case
	  xk1dlast(kz,xs-1) =  d*(levels[0].buold(jy,kz)*(xk1dlast(kz,xs)-levels[0].rlold[kz]) - levels[0].blold(jy,kz)*(xk1dlast(kz,xe)-levels[0].ruold[kz]));
	  xk1dlast(kz,xe+1) = -d*(levels[0].auold(jy,kz)*(xk1dlast(kz,xs)-levels[0].rlold[kz]) - levels[0].alold(jy,kz)*(xk1dlast(kz,xe)-levels[0].ruold[kz]));
	} else if(localmesh->firstX() and not localmesh->lastX()) {
	  // Lower boundary but not upper boundary
	  // xk1dlast[xs-1] = already correct
	  xk1dlast(kz,xe+1) = (xk1dlast(kz,xe)-levels[0].ruold[kz])/levels[0].buold(jy,kz);
	} else if(localmesh->lastX() and not localmesh->firstX()){
	  // Upper boundary but not lower boundary
	  // xk1dlast[xe+1] = already correct
	  xk1dlast(kz,xs-1) = (xk1dlast(kz,xs)-levels[0].rlold[kz])/levels[0].alold(jy,kz);
	} 
	// No "else" case. If both upper and lower boundaries, both xs-1 and xe+1
	// are already correct
	SCOREP_USER_REGION_END(getnewmethodsoln);
      }
    }

    // Now that halo cells are converged, use these to calculate whole solution
    reconstruct_full_solution(levels[0],jy,xk1dlast);
  }

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
  SCOREP_USER_REGION_END(afterloop);

  SCOREP_USER_REGION_DEFINE(fftback);
  SCOREP_USER_REGION_BEGIN(fftback, "fft back",SCOREP_USER_REGION_TYPE_COMMON);
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
  SCOREP_USER_REGION_END(fftback);
  return x; // Result of the inversion
}

void LaplaceParallelTriMG::jacobi(Level &l, const int jy, const Array<bool> &converged){

  SCOREP0();
  struct Message { dcomplex value; bool done; };
  Array<Message> message_send, message_recv;
  message_send = Array<Message>(nmode);
  message_recv = Array<Message>(nmode);
  MPI_Comm comm = BoutComm::get();
  int err;

  // Only need to update interior points
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      l.xloclast(0,kz) = l.xloc(0,kz);
      l.xloclast(1,kz) = l.xloc(1,kz);
      l.xloclast(2,kz) = l.xloc(2,kz);
      l.xloclast(3,kz) = l.xloc(3,kz);
      if( localmesh->firstX() ){
	l.xloc(0,kz) = ( l.minvb(kz,l.xs-1) + l.upperGuardVector(l.xs-1,jy,kz)*l.xloclast(3,kz) ) / (1.0 - l.lowerGuardVector(l.xs-1,jy,kz));
      }
      l.xloc(1,kz) = l.rl[kz] + l.al(jy,kz)*l.xloclast(0,kz) + l.bl(jy,kz)*l.xloclast(3,kz);
      l.xloc(2,kz) = l.ru[kz] + l.au(jy,kz)*l.xloclast(0,kz) + l.bu(jy,kz)*l.xloclast(3,kz);
      if( localmesh->lastX() ){
	l.xloc(3,kz) = ( l.minvb(kz,l.xe+1) + l.lowerGuardVector(l.xe+1,jy,kz)*l.xloclast(0,kz) ) / (1.0 - l.upperGuardVector(l.xe+1,jy,kz));
      }
    }
  }

    //output<<"after work jy, count "<<jy<<" "<<count<<endl;
    SCOREP_USER_REGION_DEFINE(comms);
    SCOREP_USER_REGION_BEGIN(comms, "communication",SCOREP_USER_REGION_TYPE_COMMON);

    // Communication

      // TODO These for loops do buffer (un)packing for data we don't care about
      // Guard? Or move to work loop?
      if(!localmesh->firstX()){
      for (int kz = 0; kz <= maxmode; kz++) {
///	if(!neighbour_in[kz]){
	  message_send[kz].value = l.xloc(index_in,kz);
///	  message_send[kz].done  = self_in[kz];
///	}
      }
      err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 1, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 0, comm, MPI_STATUS_IGNORE);
      for (int kz = 0; kz <= maxmode; kz++) {
///	if(!self_in[kz]){
	  l.xloc(0,kz) = message_recv[kz].value;
///	  neighbour_in[kz] = message_recv[kz].done;
///	}
      }
///    }
}

    // Communicate out
    // See note above for inward communication.
//    TODO Guard comms
///    if(!all(neighbour_out)) {
      //output<<"neighbour_out proc "<<BoutComm::rank()<<endl;
      if(!localmesh->lastX()){
      for (int kz = 0; kz <= maxmode; kz++) {
	message_send[kz].value = l.xloc(index_out,kz);
///	message_send[kz].done  = self_out[kz];
      }
      err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 0, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 1, comm, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < nmode; kz++) {
	l.xloc(3,kz) = message_recv[kz].value;
///	neighbour_out[kz] = message_recv[kz].done;
      }
      }
///    }
    SCOREP_USER_REGION_END(comms);
}  

void LaplaceParallelTriMG::gauss_seidel_red_black_full_system(Level &l, const Array<bool> &converged, const int jy){

  SCOREP0();
  Array<dcomplex> sendvec, recvec;
  sendvec = Array<dcomplex>(nmode);
  recvec = Array<dcomplex>(nmode);
  comm_handle recv[1];

  // Red sweep: odd points
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      for (int ix = l.xs; ix < l.xe+1; ix+=2) {
	l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1) ) / l.bvec(jy,kz,ix);
      }
    }
  }

  // Communicate: final grid point needs data from proc above
  if(!localmesh->lastX()){
    recv[0] = localmesh->irecvXOut(&recvec[0], nmode, 1);
  }
  if(!localmesh->firstX()){
    for(int kz=0; kz < nmode; kz++){
      if(!converged[kz]){
	sendvec[kz] = l.soln(kz,l.xs);
      }
    }
    localmesh->sendXIn(&sendvec[0],nmode,1);
  }
  if(!localmesh->lastX()){
    localmesh->wait(recv[0]);
    for(int kz=0; kz < nmode; kz++){
      if(!converged[kz]){
	l.soln(kz,l.xe+1) = recvec[kz];
      }
    }
  }

  // Black sweep: even points
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      for (int ix = l.xs+1; ix <= l.xe; ix+=2) {
	l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1) ) / l.bvec(jy,kz,ix);
      }
    }
  }

  // Communicate: to synchronize, first grid point needs data from proc below
  if(!localmesh->firstX()){
    recv[0] = localmesh->irecvXIn(&recvec[0], nmode, 1);
  }
  if(!localmesh->lastX()){
    for(int kz=0; kz < nmode; kz++){
      if(!converged[kz]){
	sendvec[kz] = l.soln(kz,l.xe);
      }
    }
    localmesh->sendXOut(&sendvec[0],nmode,1);
  }
  if(!localmesh->firstX()){
    localmesh->wait(recv[0]);
    for(int kz=0; kz < nmode; kz++){
      if(!converged[kz]){
	l.soln(kz,l.xs-1) = recvec[kz];
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

  /*
  for (int ix = 0; ix < l.ncx; ix++) {
    output << l.soln(0,ix).real()<<" ";
  }
  output<<endl;
  */

}  


/*
 * Perform a Gauss--Seidel iteration with red black colouring explicitly on the full system 
 * Note that this assumes that each processor has an ever number of points
 */
void LaplaceParallelTriMG::gauss_seidel_red_black_full_system_comp_comm_overlap(Level &l, const Array<bool> &converged, const int jy){

  SCOREP0();
  Array<dcomplex> sendvecred, recvecred, sendvecblack, recvecblack;
  sendvecred = Array<dcomplex>(nmode);
  recvecred = Array<dcomplex>(nmode);
  sendvecblack = Array<dcomplex>(nmode);
  recvecblack = Array<dcomplex>(nmode);
  //comm_handle recv[1], send[1];
  MPI_Request rredreq, sredreq, rblackreq, sblackreq;

  // Communicate: final grid point needs data from proc above
  // Overlap comm/comp by not posting recv until needed in final
  // loop iteration
  if(!localmesh->lastX()){
    //recv[0] = localmesh->irecvXOut(&recvec[0], nmode, 1);
    MPI_Irecv(&recvecred[0], nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, BoutComm::get(), &rredreq);
  }
  // Communicate: to synchronize, first grid point needs data from proc below
  if(!localmesh->firstX()){
    //recv[0] = localmesh->irecvXIn(&recvec[0], nmode, 1);
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

  bool wait_red = false; // have we posted the wait?
  // Black sweep: even points
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      for (int ix = l.xs+1; ix < l.xe; ix+=2) {
	l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1) ) / l.bvec(jy,kz,ix);
      }
    }
  }
  for (int kz = 0; kz <= maxmode; kz++) {
    if(!converged[kz]){
      if(!localmesh->lastX()){
        if(not wait_red){
          wait_red = true;
          //localmesh->wait(recv[0]);
	  MPI_Wait(&rredreq,MPI_STATUS_IGNORE);
        }
        l.soln(kz,l.xe+1) = recvecred[kz];
      }
      {
        int ix = l.xe;
	l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1) ) / l.bvec(jy,kz,ix);
      }
    }
  }

  if(!localmesh->lastX()){
    for(int kz=0; kz < nmode; kz++){
      if(!converged[kz]){
	sendvecblack[kz] = l.soln(kz,l.xe);
      }
    }
    //localmesh->sendXOut(&sendvec[0],nmode,2);
    MPI_Isend(&sendvecblack[0], nmode, MPI_DOUBLE_COMPLEX, proc_out, 2, BoutComm::get(), &sblackreq);
  }

  bool wait_black = false; // have we posted the wait?
  if(!localmesh->firstX()){
    for(int kz=0; kz < nmode; kz++){
      if(!converged[kz]){
        if(not wait_black){
	  wait_black = true;
	  //localmesh->wait(recv[0]);
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
  /*
  for (int ix = 0; ix < l.ncx; ix++) {
    output << l.soln(0,ix).real()<<" ";
  }
  output<<endl;
  */

}  

// Perform a Jacobi iteration explicitly on the full system 
void LaplaceParallelTriMG::jacobi_full_system(Level &l, const int jy){

  SCOREP0();
  struct Message { dcomplex value; bool done; };
  Array<Message> message_send, message_recv;
  message_send = Array<Message>(nmode);
  message_recv = Array<Message>(nmode);
  MPI_Comm comm = BoutComm::get();
  int err;

  for (int kz = 0; kz <= maxmode; kz++) {
    // TODO guard work for converged kz
    l.soln(kz,l.xs-1) = l.solnlast(kz,l.xs-1);
    for (int ix = l.xs; ix < l.xe+1; ix++) {
      l.soln(kz,ix) = ( l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.solnlast(kz,ix-1) - l.cvec(jy,kz,ix)*l.solnlast(kz,ix+1) ) / l.bvec(jy,kz,ix);
    }

    if(!localmesh->firstX()){
      for (int kz = 0; kz <= maxmode; kz++) {
///	if(!neighbour_in[kz]){
	  message_send[kz].value = l.soln(kz,l.xs);
///	  message_send[kz].done  = self_in[kz];
///	}
      }
      err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 1, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 0, comm, MPI_STATUS_IGNORE);
      for (int kz = 0; kz <= maxmode; kz++) {
///	if(!self_in[kz]){
	  l.soln(kz,l.xs-1) = message_recv[kz].value;
///	  neighbour_in[kz] = message_recv[kz].done;
///	}
      }
///    }
}

    // Communicate out
    // See note above for inward communication.
//    TODO Guard comms
///    if(!all(neighbour_out)) {
      //output<<"neighbour_out proc "<<BoutComm::rank()<<endl;
    if(!localmesh->lastX()){
      for (int kz = 0; kz <= maxmode; kz++) {
	message_send[kz].value = l.soln(kz,l.xe);
///	message_send[kz].done  = self_out[kz];
      }
      err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 0, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 1, comm, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < nmode; kz++) {
	l.soln(kz,l.xe+1) = message_recv[kz].value;
///	neighbour_out[kz] = message_recv[kz].done;
      }
      }


    // Update boundaries to match interior points
    // Do this after communication, otherwise this breaks on 1 interior pt per proc
    if(l.current_level == 0){
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

// Write info about levels to screen
void LaplaceParallelTriMG::levels_info(const Level l, const int jy){

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
void LaplaceParallelTriMG::init(Level &l, const Level lup, int ncx, const int xs, const int xe, const int current_level, const int jy){

  SCOREP0();
  l.xs = xs;
  l.xe = xe;
  l.ncx = ncx;
  l.current_level = current_level;
  int ny = localmesh->LocalNy;

  if(l.xe-l.xs<1){
    throw BoutException("LaplaceParallelTriMG error: Coarse grids must contain at least two points on every processor. Please set max_level smaller than %i, use fewer processors, or increase x resolution.",l.current_level);
    // Note: Grids are initialized from finest to coarsest, so l.current_level
    // is the finest grid that fails.
  }

  l.avec = Tensor<dcomplex>(ny,nmode,ncx);
  l.bvec = Tensor<dcomplex>(ny,nmode,ncx);
  l.cvec = Tensor<dcomplex>(ny,nmode,ncx);
  l.rvec = Matrix<dcomplex>(nmode,ncx);
  l.residual = Matrix<dcomplex>(nmode,ncx);
  l.soln = Matrix<dcomplex>(nmode,ncx);
  l.solnlast = Matrix<dcomplex>(nmode,ncx);
  if(algorithm!=0){
    l.xloc = Matrix<dcomplex>(4,nmode);
    l.xloclast = Matrix<dcomplex>(4,nmode);
  }

  // For coarsening, need to communicate coefficients for halo cells
  l.acomm = Array<dcomplex>(nmode);
  l.bcomm = Array<dcomplex>(nmode);
  l.ccomm = Array<dcomplex>(nmode);

  Array<dcomplex> sendvec, recvec;
  sendvec = Array<dcomplex>(3*nmode);
  recvec = Array<dcomplex>(3*nmode);
  comm_handle recv[1];

  if(!localmesh->firstX()){
    recv[0] = localmesh->irecvXIn(&recvec[0], 3*nmode, 1);
  }
  if(!localmesh->lastX()){
    for(int kz=0; kz < nmode; kz++){
      sendvec[kz] = lup.avec(jy,kz,lup.xe);
    }
    for(int kz=0; kz < nmode; kz++){
      sendvec[nmode+kz] = lup.bvec(jy,kz,lup.xe);
    }
    for(int kz=0; kz < nmode; kz++){
      sendvec[2*nmode+kz] = lup.cvec(jy,kz,lup.xe);
    }
    localmesh->sendXOut(&sendvec[0],3*nmode,1);
  }
  if(!localmesh->firstX()){
    localmesh->wait(recv[0]);
    for(int kz=0; kz < nmode; kz++){
      l.acomm[kz] = recvec[kz];
      l.bcomm[kz] = recvec[kz+nmode];
      l.ccomm[kz] = recvec[kz+2*nmode];
    }
  }

  for(int kz = 0; kz < nmode; kz++){
    for(int ix = 0; ix<l.xs; ix++){
      if(localmesh->firstX()){
	l.avec(jy,kz,ix) = 0.5*lup.avec(jy,kz,ix);
	l.bvec(jy,kz,ix) = 0.5*lup.bvec(jy,kz,ix);
	l.cvec(jy,kz,ix) = 0.5*lup.cvec(jy,kz,ix);
      }
      else{
	l.avec(jy,kz,ix) = lup.avec(jy,kz,ix);
	l.bvec(jy,kz,ix) = lup.bvec(jy,kz,ix);
	l.cvec(jy,kz,ix) = lup.cvec(jy,kz,ix);
      }
    }
    // interior points
    for(int ixc = l.xs; ixc<l.xe+1; ixc++){
      int ixf = 2*(ixc-l.xs)+l.xs;
      if(localmesh->firstX() and ixc == l.xs){
	// No lumping in avec for first interior point:
	// The gap between this point and the first guard cell is NOT doubled when the mesh is refined. 
	l.avec(jy,kz,ixc) = 0.5*lup.avec(jy,kz,ixf);
	l.bvec(jy,kz,ixc) = 0.5*lup.bvec(jy,kz,ixf) + 0.25*lup.cvec(jy,kz,ixf) + 0.25*lup.avec(jy,kz,ixf+1) + 0.125*lup.bvec(jy,kz,ixf+1);
	l.cvec(jy,kz,ixc) = 0.25*lup.cvec(jy,kz,ixf) + 0.125*lup.bvec(jy,kz,ixf+1) + 0.25*lup.cvec(jy,kz,ixf+1);
      }
      else if(ixc == l.xs){
	// need a/b/cvec from the proc below
	l.avec(jy,kz,ixc) = 0.25*l.acomm[kz] + 0.125*l.bcomm[kz] + 0.25*lup.avec(jy,kz,ixf) ;
	l.bvec(jy,kz,ixc) = 0.125*l.bcomm[kz] + 0.25*l.ccomm[kz] + 0.25*lup.avec(jy,kz,ixf) + 0.5*lup.bvec(jy,kz,ixf) + 0.25*lup.cvec(jy,kz,ixf) + 0.25*lup.avec(jy,kz,ixf+1) + 0.125*lup.bvec(jy,kz,ixf+1);
	l.cvec(jy,kz,ixc) = 0.25*lup.cvec(jy,kz,ixf) + 0.125*lup.bvec(jy,kz,ixf+1) +  0.25*lup.cvec(jy,kz,ixf+1); 
      }
      else{
	l.avec(jy,kz,ixc) = 0.25*lup.avec(jy,kz,ixf-1) + 0.125*lup.bvec(jy,kz,ixf-1) + 0.25*lup.avec(jy,kz,ixf) ;
	l.bvec(jy,kz,ixc) = 0.125*lup.bvec(jy,kz,ixf-1) + 0.25*lup.cvec(jy,kz,ixf-1) + 0.25*lup.avec(jy,kz,ixf) + 0.5*lup.bvec(jy,kz,ixf) + 0.25*lup.cvec(jy,kz,ixf) + 0.25*lup.avec(jy,kz,ixf+1) + 0.125*lup.bvec(jy,kz,ixf+1);
	l.cvec(jy,kz,ixc) = 0.25*lup.cvec(jy,kz,ixf) + 0.125*lup.bvec(jy,kz,ixf+1) +  0.25*lup.cvec(jy,kz,ixf+1); 
      }
    }
    for(int ixc = l.xe+1; ixc<l.ncx; ixc++){
      // Index on fine grid
      int ixf = ixc + lup.ncx - l.ncx;
      if( localmesh->lastX() ){
	if( ixc == l.xe+1){
	  // Lump avec on first physical boundary point:
	  // The grid spacing has been doubled here
	  l.avec(jy,kz,ixc) =  0.25*lup.avec(jy,kz,ixf-1) + 0.125*lup.bvec(jy,kz,ixf-1) + 0.25*lup.avec(jy,kz,ixf);
	  l.bvec(jy,kz,ixc) = 0.125*lup.bvec(jy,kz,ixf-1) +  0.25*lup.cvec(jy,kz,ixf-1) + 0.25*lup.avec(jy,kz,ixf) + 0.5*lup.bvec(jy,kz,ixf);
	  l.cvec(jy,kz,ixc) = 0.5*lup.cvec(jy,kz,ixf);
	}
	else{
	  l.avec(jy,kz,ixc) = 0.5*lup.avec(jy,kz,ixf);
	  l.bvec(jy,kz,ixc) = 0.5*lup.bvec(jy,kz,ixf);
	  l.cvec(jy,kz,ixc) = 0.5*lup.cvec(jy,kz,ixf);
	}
      }
      else{
	l.avec(jy,kz,ixc) = lup.avec(jy,kz,ixf);
	l.bvec(jy,kz,ixc) = lup.bvec(jy,kz,ixf);
	l.cvec(jy,kz,ixc) = lup.cvec(jy,kz,ixf);
      }
    }
  }

  //levels_info(l,jy);

}

// Init routine for finest level
void LaplaceParallelTriMG::init(Level &l, const int ncx, const int jy, const Matrix<dcomplex> avec, const Matrix<dcomplex> bvec, const Matrix<dcomplex> cvec, const int xs, const int xe){

  // Basic definitions for conventional multigrid
  SCOREP0();
  l.xs = xs;
  l.xe = xe;
  l.ncx = ncx;
  l.current_level = 0;
  int ny = localmesh->LocalNy;

  l.avec = Tensor<dcomplex>(ny,nmode,ncx);
  l.bvec = Tensor<dcomplex>(ny,nmode,ncx);
  l.cvec = Tensor<dcomplex>(ny,nmode,ncx);
  l.rvec = Matrix<dcomplex>(nmode,ncx);
  l.residual = Matrix<dcomplex>(nmode,ncx);
  l.soln = Matrix<dcomplex>(nmode,ncx);
  l.solnlast = Matrix<dcomplex>(nmode,ncx);

  for(int kz=0; kz<nmode; kz++){
    for(int ix=0; ix<ncx; ix++){
     l.avec(jy,kz,ix) = avec(kz,ix); 
     l.bvec(jy,kz,ix) = bvec(kz,ix); 
     l.cvec(jy,kz,ix) = cvec(kz,ix); 
     l.residual(kz,ix) = 0.0;
    }
  }
  // end basic definitions

  if(algorithm!=0){
    auto rlold = Array<dcomplex>(nmode);
    auto ruold = Array<dcomplex>(nmode);
    //auto Rd = Array<dcomplex>(ncz/2+1);
    //auto Ru = Array<dcomplex>(ncz/2+1);
    //auto Rsendup = Array<dcomplex>(ncz+2);
    //auto Rsenddown = Array<dcomplex>(ncz+2);
    //auto Rrecvup = Array<dcomplex>(ncz+2);
    //auto Rrecvdown = Array<dcomplex>(ncz+2);
    auto evec = Array<dcomplex>(ncx);
    auto tmp = Array<dcomplex>(ncx);

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

    // Coefs used to compute rl, ru in extended domain method
    l.r1 = Matrix<dcomplex>(localmesh->LocalNy, nmode);
    l.r2 = Matrix<dcomplex>(localmesh->LocalNy, nmode);
    l.r3 = Matrix<dcomplex>(localmesh->LocalNy, nmode);
    l.r4 = Matrix<dcomplex>(localmesh->LocalNy, nmode);
    l.r5 = Matrix<dcomplex>(localmesh->LocalNy, nmode);
    l.r6 = Matrix<dcomplex>(localmesh->LocalNy, nmode);
    l.r7 = Matrix<dcomplex>(localmesh->LocalNy, nmode);
    l.r8 = Matrix<dcomplex>(localmesh->LocalNy, nmode);

    for (int kz = 0; kz <= maxmode; kz++) {

      SCOREP_USER_REGION_DEFINE(invert);
      SCOREP_USER_REGION_BEGIN(invert, "invert local matrices",SCOREP_USER_REGION_TYPE_COMMON);

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

      SCOREP_USER_REGION_END(invert);
      SCOREP_USER_REGION_DEFINE(coefs);
      SCOREP_USER_REGION_BEGIN(coefs, "calculate coefs",SCOREP_USER_REGION_TYPE_COMMON);

      l.bl(jy,kz) = l.upperGuardVector(l.xs,jy,kz);
      l.al(jy,kz) = l.lowerGuardVector(l.xs,jy,kz);

      l.bu(jy,kz) = l.upperGuardVector(l.xe,jy,kz);
      l.au(jy,kz) = l.lowerGuardVector(l.xe,jy,kz);

      l.alold(jy,kz) = l.al(jy,kz);
      l.auold(jy,kz) = l.au(jy,kz);
      l.blold(jy,kz) = l.bl(jy,kz);
      l.buold(jy,kz) = l.bu(jy,kz);

      // New method - connect to more distant points
      if(new_method){

	// First compute coefficients that depend on the matrix to be inverted
	// and which therefore might be constant throughout a run.

	// Boundary processor values to be overwritten when relevant
	dcomplex Ad, Bd, Au, Bu, Atmp, Btmp;
	Ad = 1.0;
	Bd = 0.0;
	Au = 0.0;
	Bu = 1.0;
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
	Delta = 1.0 - l.al(jy,kz)*Bd - l.bu(jy,kz)*Au + (l.al(jy,kz)*l.bu(jy,kz) - l.au(jy,kz)*l.bl(jy,kz))*Bd*Au;
	Delta = 1.0 / Delta;
	l.al(jy,kz) = Delta*( l.alold(jy,kz) + (l.auold(jy,kz)*l.blold(jy,kz) - l.alold(jy,kz)*l.buold(jy,kz))*Au )*Ad;
	l.bl(jy,kz) = Delta * l.blold(jy,kz) * Bu ;
	l.au(jy,kz) = Delta * l.auold(jy,kz) * Ad ;
	l.bu(jy,kz) = Delta*( l.buold(jy,kz) + (l.auold(jy,kz)*l.blold(jy,kz) - l.alold(jy,kz)*l.buold(jy,kz))*Bd )*Bu;

	dcomplex d = l.auold(jy,kz)*l.blold(jy,kz) - l.alold(jy,kz)*l.buold(jy,kz);
	l.r1(jy,kz) = Delta*(l.alold(jy,kz) + d*Au);
	l.r2(jy,kz) = Delta*( 1.0 - l.buold(jy,kz)*Au );
	l.r3(jy,kz) = Delta*l.blold(jy,kz)*Au;
	l.r4(jy,kz) = Delta*l.blold(jy,kz);
	l.r5(jy,kz) = Delta*l.auold(jy,kz);
	l.r6(jy,kz) = Delta*l.auold(jy,kz)*Bd;
	l.r7(jy,kz) = Delta*( 1.0 - l.alold(jy,kz)*Bd );
	l.r8(jy,kz) = Delta*(l.buold(jy,kz) + d*Bd);

      } // new method
      SCOREP_USER_REGION_END(coefs);
    } // end of kz loop
  }

  //levels_info(l,jy);

}

// Init routine for finest level information that cannot be cached
void LaplaceParallelTriMG::init_rhs(Level &l, const int jy, const Matrix<dcomplex> bcmplx){

  SCOREP0();

  if(algorithm!=0){
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

    for (int kz = 0; kz <= maxmode; kz++) {

      SCOREP_USER_REGION_DEFINE(invertforrhs);
      SCOREP_USER_REGION_BEGIN(invertforrhs, "invert local matrices for rhs",SCOREP_USER_REGION_TYPE_COMMON);

      // Invert local matrices
      // Calculate Minv*b
      tridag(&l.avec(jy,kz,0), &l.bvec(jy,kz,0), &l.cvec(jy,kz,0), &bcmplx(kz,0),
      &l.minvb(kz,0), l.ncx);
      // Now minvb is a constant vector throughout the iterations

      SCOREP_USER_REGION_END(invertforrhs);
      SCOREP_USER_REGION_DEFINE(coefsforrhs);
      SCOREP_USER_REGION_BEGIN(coefsforrhs, "calculate coefs for rhs",SCOREP_USER_REGION_TYPE_COMMON);

      l.rl[kz] = l.minvb(kz,l.xs);
      l.ru[kz] = l.minvb(kz,l.xe);
      l.rlold[kz] = l.rl[kz];
      l.ruold[kz] = l.ru[kz];

      // New method - connect to more distant points
      if(new_method){

	// Compute coefficients that depend on the right-hand side and
	// which therefore change every time.

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
      } // new method
      SCOREP_USER_REGION_END(coefsforrhs);
    } // end of kz loop

    // Communicate vector in kz
    if(new_method){
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
	l.rl[kz] = l.r1(jy,kz)*Rd[kz] + l.r2(jy,kz)*l.rlold[kz] + l.r3(jy,kz)*l.ruold[kz] + l.r4(jy,kz)*Ru[kz] ;
	l.ru[kz] = l.r5(jy,kz)*Rd[kz] + l.r6(jy,kz)*l.rlold[kz] + l.r7(jy,kz)*l.ruold[kz] + l.r8(jy,kz)*Ru[kz] ;
      }
    }
  }

  //levels_info(l,jy);

}

/*
 * Sum and communicate total residual
 */
void LaplaceParallelTriMG::calculate_total_residual(Array<BoutReal> &error_abs, Array<BoutReal> &error_rel, Array<bool> &converged, Level &l){

  SCOREP0();
  // Communication arrays:
  // residual in (0 .. nmode-1)
  // solution in (nmode .. 2*nmode-1)
  auto subtotal = Array<BoutReal>(2*nmode);
  auto total = Array<BoutReal>(2*nmode);

  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      error_abs[kz] = 0.0;
      total[kz] = 0.0;
      total[kz+nmode] = 0.0;
      // Only xs and xe have nonzero residuals
      subtotal[kz] = pow(l.residual(kz,l.xs).real(),2) + pow(l.residual(kz,l.xs).imag(),2) + pow(l.residual(kz,l.xe).real(),2) + pow(l.residual(kz,l.xe).imag(),2);
      // Strictly this should be all contributions to the solution, but this under-approximation saves work. Could multiply by (interior points/2)
      subtotal[kz+nmode] = pow(l.xloc(1,kz).real(),2) + pow(l.xloc(1,kz).imag(),2) + pow(l.xloc(2,kz).real(),2) + pow(l.xloc(2,kz).imag(),2);
    }
  }

  // Communication needed to ensure processorsbreak iterations at same point
  // TODO do this rarely
  MPI_Allreduce(&subtotal[0], &total[0], 2*nmode, MPI_DOUBLE, MPI_SUM, BoutComm::get());

  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      error_abs[kz] = sqrt(total[kz]);
      error_rel[kz] = error_abs[kz]/sqrt(total[kz+nmode]);
      //if( error_rel[kz] < rtol){
      if( error_abs[kz] < atol or error_rel[kz] < rtol ){
	//output<<"CONVERGED! "<<kz<<" "<<error_abs[kz]<<" "<<error_rel[kz]<<endl;
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
 * Sum and communicate total residual
 */
void LaplaceParallelTriMG::calculate_total_residual_full_system(Array<BoutReal> &error_abs, Array<BoutReal> &error_rel, Array<bool> &converged, const Level l){

  SCOREP0();
  // Communication arrays:
  // residual in (0 .. nmode-1)
  // solution in (nmode .. 2*nmode-1)
  auto subtotal = Array<BoutReal>(2*nmode);
  auto total = Array<BoutReal>(2*nmode);

  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      error_abs[kz] = 0.0;
      total[kz] = 0.0;
      total[kz+nmode] = 0.0;
      subtotal[kz] = 0.0;
      subtotal[kz+nmode] = 0.0;
      for(int ix=0; ix<l.ncx;ix++){
	// Contributions to total residual from interior points and physical
	// boundary points only
	if( (ix < l.xs and localmesh->firstX()) or
	    (ix >= l.xs and ix <= l.xe) or
	    (ix > l.xe and localmesh->lastX()) ){
	  subtotal[kz] += pow(l.residual(kz,ix).real(),2) + pow(l.residual(kz,ix).imag(),2);
	  subtotal[kz+nmode] += pow(l.soln(kz,ix).real(),2) + pow(l.soln(kz,ix).imag(),2);
	}  
      }
    }
  }

  // Communication needed to ensure processorsbreak iterations at same point
  // TODO do this rarely
  MPI_Allreduce(&subtotal[0], &total[0], 2*nmode, MPI_DOUBLE, MPI_SUM, BoutComm::get());

  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      error_abs[kz] = sqrt(total[kz]);
      error_rel[kz] = error_abs[kz]/sqrt(total[kz+nmode]);
      //if( error_rel[kz] < rtol){
      if( error_abs[kz] < atol or error_rel[kz] < rtol ){
	converged[kz] = true;
      }
    }
  }
}

// Calculate residual on a reduced x grid. By construction, the residual is 
// zero, except at the points on the reduced grid.
void LaplaceParallelTriMG::calculate_residual(Level &l, const Array<bool> &converged,const int jy){

  SCOREP0();
  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      //if(localmesh->firstX()){
	l.residual(kz,0) = 0.0;
	l.residual(kz,1) = 0.0;
	l.residual(kz,l.xs) = l.avec(jy,kz,l.xs)*(l.xloc(0,kz)-l.xloclast(0,kz));
	for(int ix=l.xs+1; ix<l.xe; ix++){
	  l.residual(kz,ix) = 0.0;
	}
	l.residual(kz,l.xe) = l.cvec(jy,kz,l.xe)*(l.xloc(3,kz)-l.xloclast(3,kz));
	l.residual(kz,l.xe+1) = 0.0;
	l.residual(kz,l.xe+2) = 0.0;
      //}
    }
  }

  struct Message { dcomplex value; bool done; };
  Array<Message> message_send, message_recv;
  message_send = Array<Message>(nmode);
  message_recv = Array<Message>(nmode);
  MPI_Comm comm = BoutComm::get();
  int err;
  if(!localmesh->firstX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	message_send[kz].value = l.residual(kz,l.xs);
      }
    }
    err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 1, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 0, comm, MPI_STATUS_IGNORE);
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	l.residual(kz,l.xs-1) = message_recv[kz].value;
      }
    }
  }

  // Communicate out
  if(!localmesh->lastX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	message_send[kz].value = l.residual(kz,l.xe);
      }
    }
    err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 0, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 1, comm, MPI_STATUS_IGNORE);
    for (int kz = 0; kz < nmode; kz++) {
      if(!converged[kz]){
	l.residual(kz,l.xe+1) = message_recv[kz].value;
      }
    }
  }
}

// Calculate residual on a full x grid
void LaplaceParallelTriMG::calculate_residual_full_system(Level &l, const Array<bool> &converged,const int jy){

  SCOREP0();
  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      if(localmesh->firstX()){
	l.residual(kz,0) = l.rvec(kz,0) - l.bvec(jy,kz,0)*l.soln(kz,0) - l.cvec(jy,kz,0)*l.soln(kz,1);
	for(int ix=1; ix<l.xs; ix++){
	  l.residual(kz,ix) = l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.bvec(jy,kz,ix)*l.soln(kz,ix) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1);
	}
      }
      for(int ix=l.xs; ix<l.xe+1; ix++){
	l.residual(kz,ix) = l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.bvec(jy,kz,ix)*l.soln(kz,ix) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1);
      }
      if(localmesh->lastX()){
	for(int ix=l.xe+1; ix<l.ncx-1; ix++){
	  l.residual(kz,ix) = l.rvec(kz,ix) - l.avec(jy,kz,ix)*l.soln(kz,ix-1) - l.bvec(jy,kz,ix)*l.soln(kz,ix) - l.cvec(jy,kz,ix)*l.soln(kz,ix+1);
	}
	l.residual(kz,l.ncx-1) = l.rvec(kz,l.ncx-1) - l.avec(jy,kz,l.ncx-1)*l.soln(kz,l.ncx-2) - l.bvec(jy,kz,l.ncx-1)*l.soln(kz,l.ncx-1);
      }
    }
  }

  struct Message { dcomplex value; bool done; };
  Array<Message> message_send, message_recv;
  message_send = Array<Message>(nmode);
  message_recv = Array<Message>(nmode);
  MPI_Comm comm = BoutComm::get();
  int err;
  if(!localmesh->firstX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	message_send[kz].value = l.residual(kz,l.xs);
      }
    }
    err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 1, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_in, 0, comm, MPI_STATUS_IGNORE);
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	l.residual(kz,l.xs-1) = message_recv[kz].value;
      }
    }
  }

  // Communicate out
  if(!localmesh->lastX()){
    for (int kz = 0; kz <= maxmode; kz++) {
      if(!converged[kz]){
	message_send[kz].value = l.residual(kz,l.xe);
      }
    }
    err = MPI_Sendrecv(&message_send[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 0, &message_recv[0], nmode*sizeof(Message), MPI_BYTE, proc_out, 1, comm, MPI_STATUS_IGNORE);
    for (int kz = 0; kz < nmode; kz++) {
      if(!converged[kz]){
	l.residual(kz,l.xe+1) = message_recv[kz].value;
      }
    }
  }
}

void LaplaceParallelTriMG::coarsen_full_system(Level &l, const Matrix<dcomplex> &fine_residual, const Array<bool> &converged){

  SCOREP0();
  int ixc, ixf;
  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      if(localmesh->firstX()){
	for(int ix=0; ix<l.xs; ix++){
	  l.residual(kz,ix) = 0.5*fine_residual(kz,ix);
	  //output<<ix<<endl;
	}
      }
      ixc = l.xs;
      ixf = l.xs;
      if(localmesh->firstX()){
	l.residual(kz,ixc)   = 0.5*fine_residual(kz,ixf) + 0.25*fine_residual(kz,ixf+1);
      }
      else{
	l.residual(kz,ixc)   =  0.25*fine_residual(kz,ixf-1) + 0.5*fine_residual(kz,ixf)   + 0.25*fine_residual(kz,ixf+1);
      }
      for(int ixc=l.xs+1; ixc<l.xe+1; ixc++){
	ixf = 2*(ixc-l.xs)+l.xs;
	l.residual(kz,ixc)   =  0.25*fine_residual(kz,ixf-1) + 0.5*fine_residual(kz,ixf)   + 0.25*fine_residual(kz,ixf+1);
	//output<<ixc<<" "<<ixf<<endl;
      }
      if(localmesh->lastX()){
	// first boundary point
	ixc = l.xe+1;
	ixf = l.xs+2*(l.xe+1-l.xs);
	l.residual(kz,ixc) = 0.25*fine_residual(kz,ixf-1) + 0.5*fine_residual(kz,ixf);
	//output<<"fbp "<<ixc<<" "<<ixf<<endl;
	// FIXME this assumes mgx=2
	for(int ixc=l.xe+2; ixc<l.ncx; ixc++){
	  ixf = l.xs+2*(l.xe-l.xs)+(ixc-l.xe)+1;
	  l.residual(kz,ixc) = 0.5*fine_residual(kz,ixf);
	  //output<<ixc<<" "<<ixf<<endl;
	}
      }
      for(int ix=0; ix<l.ncx; ix++){
	l.rvec(kz,ix) = l.residual(kz,ix);
	l.soln(kz,ix) = 0.0;
	l.solnlast(kz,ix) = 0.0;
      }
    }
  }

  /*
  output<<"coarse residual ";
  for(int ix=0; ix<l.ncx; ix++){
    output<<l.rvec(2,ix)<<" ";
  }
  output<<endl;
  */

}

void LaplaceParallelTriMG::coarsen(const Level l, Matrix<dcomplex> &xloc, Matrix<dcomplex> &xloclast, int jy){

  SCOREP0();
  MPI_Comm comm = BoutComm::get();
  Array<dcomplex> tmpsend, tmprecv;
  MPI_Request request[1];
  tmpsend = Array<dcomplex>(2*nmode);
  tmprecv = Array<dcomplex>(2*nmode);

  if(!localmesh->firstX()){
    MPI_Irecv(&tmprecv[0], 2*nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, comm, &request[0]);
  }

  for(int kz=0; kz<nmode; kz++){

    // Reconstruct required x point
    // xloc[1] and xloc[3] are the same point, no manipulation needed
    // An xloc[0] on a physical boundary does not move, no manipulation
    // Otherwise xloc[0] must be received from the processor below
    // xloc[2] always moves 1 (fine) grid point to the left, so must be recalculated and sent upwards
    xloc(2,kz) = l.minvb(kz,l.xe-1);
    xloclast(2,kz) = l.minvb(kz,l.xe-1);
    if(!localmesh->lastX()){
      xloc(2,kz) += l.upperGuardVector(l.xe-1,jy,kz)*xloc(3,kz);
      xloclast(2,kz) += l.upperGuardVector(l.xe-1,jy,kz)*xloclast(3,kz);
    }
    if(!localmesh->firstX()){
      xloc(2,kz) += l.lowerGuardVector(l.xe-1,jy,kz)*xloc(0,kz);
      xloclast(2,kz) += l.lowerGuardVector(l.xe-1,jy,kz)*xloclast(0,kz);
    }
  }
  if(!localmesh->lastX()){
    // Send upwards
    for(int kz=0; kz<nmode; kz++){
      tmpsend[kz] = xloc(2,kz);
    }
    for(int kz=0; kz<nmode; kz++){
      tmpsend[nmode+kz] = xloclast(2,kz);
    }
    MPI_Isend(&tmpsend[0], 2*nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, comm, &request[0]);
  }
  if(!localmesh->firstX()){
    MPI_Wait(&request[0],MPI_STATUS_IGNORE);
    for(int kz=0; kz<nmode; kz++){
      xloc(0,kz) = tmprecv[kz];
    }
    for(int kz=0; kz<nmode; kz++){
      xloclast(0,kz) = tmprecv[nmode+kz];
    }
  }
}

// Update the solution on the refined grid by adding the error calculated on the coarser grid.
void LaplaceParallelTriMG::update_solution(Level &l, const Matrix<dcomplex> &fine_error, const Array<bool> &converged){

  SCOREP0();
  /*
  output<<"sbu ";
  for(int ix=0; ix<l.ncx; ix++){
    output<<l.soln(0,ix)<<" ";
  }
  output<<endl;
  output<<"fiu ";
  for(int ix=0; ix<l.ncx; ix++){
    output<<fine_error(0,ix)<<" ";
  }
  output<<endl;
  */
  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      for(int ix=0; ix<l.ncx; ix++){
	l.soln(kz,ix) += fine_error(kz,ix);
      }
    }
  }
  /*
  output<<"sau ";
  for(int ix=0; ix<l.ncx; ix++){
    output<<l.soln(0,ix)<<" ";
  }
  output<<endl;
  */
}

void LaplaceParallelTriMG::refine_full_system(Level &l, Matrix<dcomplex> &fine_error, const Array<bool> &converged){

  SCOREP0();
  /*
  output<<"soln ";
  for(int ix=0; ix<l.ncx; ix++){
    output<<l.soln(0,ix)<<" ";
  }
  output<<endl;
  */
  
  for(int kz=0; kz<nmode; kz++){
    if(!converged[kz]){
      if(localmesh->firstX()){
	// lower boundary (fine/coarse indices the same)
	for(int ix=0; ix<l.xs; ix++){
	  fine_error(kz,ix) = l.soln(kz,ix);
	  //output<<ix<<endl;
	}
      }
      else{
	// guard cells
	fine_error(kz,l.xs-1) = 0.5*(l.soln(kz,l.xs-1)+l.soln(kz,l.xs));
	  //output<<ix<<endl;
      }
      // interior points
      for(int ixc=l.xs; ixc<l.xe+1; ixc++){
	int ixf = 2*(ixc-l.xs)+l.xs;
	fine_error(kz,ixf) = l.soln(kz,ixc);
	fine_error(kz,ixf+1) = 0.5*(l.soln(kz,ixc)+l.soln(kz,ixc+1));
	//output<<ixc<<" "<<ixf<<endl;
      }
      // upper boundary
      for(int ixc=l.xe+1; ixc<l.ncx; ixc++){
	int ixf = l.xs + 2*(l.xe + 1 - l.xs)+ (ixc - l.xe - 1);
	fine_error(kz,ixf) = l.soln(kz,ixc);
	//output<<"ub "<<ixc<<" "<<ixf<<endl;
      }
    }
  }
  
  /*
  output<<"fine_error ";
  for(int ix=0; ix<2*(l.ncx-4)+4; ix++){
    output<<fine_error(0,ix)<<" ";
  }
  output<<endl;
  */
}

void LaplaceParallelTriMG::refine(Matrix<dcomplex> &xloc, Matrix<dcomplex> &xloclast){

  SCOREP0();
  // xloc[1] and xloc[3] don't change
  // xloc[0] unchanged if firstX, otherwise interpolated
  // xloc[2] always interpolated
  for(int kz=0; kz<nmode; kz++){
    if(!localmesh->firstX()){
      xloc(0,kz) = 0.5*(xloc(0,kz)+xloc(1,kz));
      xloclast(0,kz) = 0.5*(xloclast(0,kz)+xloclast(1,kz));
    }
    xloc(2,kz) = 0.5*(xloc(2,kz)+xloc(3,kz));

    xloclast(2,kz) = 0.5*(xloclast(2,kz)+xloclast(3,kz));
  }
}
