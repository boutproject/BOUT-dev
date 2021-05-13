/**************************************************************************
 * Perpendicular Laplacian inversion. Parallel code using FFTs in z
 * and parallel cyclic reduction in x.
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

#include "pcr.hxx"
#include "globals.hxx"

#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout/openmpwrap.hxx>
#include <bout/sys/timer.hxx>
#include <boutexception.hxx>
#include <cmath>
#include <fft.hxx>
#include <lapack_routines.hxx>
#include <utils.hxx>

#include "boutcomm.hxx"
#include <output.hxx>

#include <bout/scorepwrapper.hxx>

#include <cmath>
#include <mpi.h>
#include <cstdlib>
#include <algorithm>

using namespace std;


LaplacePCR::LaplacePCR(Options* opt, CELL_LOC loc, Mesh* mesh_in)
    : Laplacian(opt, loc, mesh_in),
      A(0.0, localmesh), C(1.0, localmesh), D(1.0, localmesh), nmode(maxmode + 1),
      ncx(localmesh->LocalNx), ny(localmesh->LocalNy), avec(ny, nmode, ncx),
      bvec(ny, nmode, ncx), cvec(ny, nmode, ncx) {

  A.setLocation(location);
  C.setLocation(location);
  D.setLocation(location);

  // Number of procs must be a power of 2
  const int n = localmesh->NXPE;
  if (!is_pow2(n)) {
    throw BoutException("LaplacePCR error: NXPE must be a power of 2");
  }

  // Cannot be run in serial
  if(localmesh->firstX() && localmesh->lastX()) {
    throw BoutException("Error: PCR method only works for NXPE > 1. Suggest using cyclic solver for NXPE = 1.\n");
  }

  // Number of x points must be a power of 2
  if (!is_pow2(localmesh->GlobalNx-4)) {
    throw BoutException("LaplacePCR error: GlobalNx must be a power of 2");
  }

  setup(localmesh->GlobalNx-4, localmesh->getNXPE(), localmesh->getXProcIndex());

}

/// Calculate the transpose of \p m in the pre-allocated \p m_t
namespace {
void transpose(Matrix<dcomplex>& m_t, const Matrix<dcomplex>& m) {
  SCOREP0();
  const auto n1 = std::get<1>(m.shape());
  const auto n2 = std::get<0>(m.shape());
  for (int i1 = 0; i1 < n1; i1++) {
    for (int i2 = 0; i2 < n2; i2++) {
      m_t(i1, i2) = m(i2, i1);
    }
  }
}
} // namespace

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
FieldPerp LaplacePCR::solve(const FieldPerp& b, const FieldPerp& x0) {

  SCOREP0();
  Timer timer("invert"); ///< Start timer

  /// SCOREP_USER_REGION_DEFINE(initvars);
  /// SCOREP_USER_REGION_BEGIN(initvars, "init vars",///SCOREP_USER_REGION_TYPE_COMMON);

  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  TRACE("LaplacePCR::solve(const, const)");

  FieldPerp x{emptyFrom(b)};

  // Info for halo swaps
  const int xproc = localmesh->getXProcIndex();
  const int yproc = localmesh->getYProcIndex();
  nproc = localmesh->getNXPE();
  myproc = yproc * nproc + xproc;
  proc_in = myproc - 1;
  proc_out = myproc + 1;

  jy = b.getIndex();

  const int ncz = localmesh->LocalNz; // Number of local z points

  xs = localmesh->xstart; // First interior point
  xe = localmesh->xend;   // Last interior point

  const BoutReal kwaveFactor = 2.0 * PI / coords->zlength();

  // Setting the width of the boundary.
  // NOTE: The default is a width of 2 guard cells
  const bool both_use_one_guard =
      isGlobalFlagSet(INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2);
  const int inbndry = (both_use_one_guard or isInnerBoundaryFlagSet(INVERT_BNDRY_ONE))
                          ? 1
                          : localmesh->xstart;
  const int outbndry = (both_use_one_guard or isOuterBoundaryFlagSet(INVERT_BNDRY_ONE))
                           ? 1
                           : localmesh->xstart;

  /* Allocation for
   * bk   = The fourier transformed of b, where b is one of the inputs in
   *        LaplacePCR::solve()
   * bk1d = The 1d array of bk
   * xk   = The fourier transformed of x, where x the output of
   *        LaplacePCR::solve()
   * xk1d = The 1d array of xk
   */
  auto bk = Matrix<dcomplex>(ncx, nmode);
  auto bk1d = Array<dcomplex>(ncx);
  auto xk = Matrix<dcomplex>(ncx, nmode);
  auto xk1d = Matrix<dcomplex>(nmode, ncx);

  /// SCOREP_USER_REGION_END(initvars);
  /// SCOREP_USER_REGION_DEFINE(fftloop);
  /// SCOREP_USER_REGION_BEGIN(fftloop, "init fft loop",SCOREP_USER_REGION_TYPE_COMMON);

  /* Coefficents in the tridiagonal solver matrix
   * Following the notation in "Numerical recipes"
   * avec is the lower diagonal of the matrix
   * bvec is the diagonal of the matrix
   * cvec is the upper diagonal of the matrix
   * NOTE: Do not confuse avec, bvec and cvec with the A, C, and D coefficients
   *       above
   */
  auto bcmplx = Matrix<dcomplex>(nmode, ncx);

  const bool invert_inner_boundary =
      isInnerBoundaryFlagSet(INVERT_SET) and localmesh->firstX();
  const bool invert_outer_boundary =
      isOuterBoundaryFlagSet(INVERT_SET) and localmesh->lastX();

  BOUT_OMP(parallel for)
  for (int ix = 0; ix < ncx; ix++) {
    /* This for loop will set the bk (initialized by the constructor)
     * bk is the z fourier modes of b in z
     * If the INVERT_SET flag is set (meaning that x0 will be used to set the
     * boundary values),
     */
    if ((invert_inner_boundary and (ix < inbndry))
        or (invert_outer_boundary and (ncx - ix - 1 < outbndry))) {
      // Use the values in x0 in the boundary
      rfft(x0[ix], ncz, &bk(ix, 0));
    } else {
      rfft(b[ix], ncz, &bk(ix, 0));
    }
  }
  /// SCOREP_USER_REGION_END(fftloop);
  /// SCOREP_USER_REGION_DEFINE(kzinit);
  /// SCOREP_USER_REGION_BEGIN(kzinit, "kz init",///SCOREP_USER_REGION_TYPE_COMMON);

  /* Solve differential equation in x for each fourier mode, so transpose to make x the
   * fastest moving index. Note that only the non-degenerate fourier modes are used (i.e.
   * the offset and all the modes up to the Nyquist frequency), so we only copy up to
   * `nmode` in the transpose.
   */
  //output.write("Before transpose\n");
  //output.write("ncx {}, nmode {}, ncz {}, ncz / 2 + 1 {}, xe-xs+1 {}  \n",ncx,nmode,ncz,ncz / 2 + 1,xe-xs+1);
  transpose(bcmplx, bk);

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
   *
   */
  //output.write("Before coefs\n");
  for (int kz = 0; kz < nmode; kz++) {
    // Note that this is called every time to deal with bcmplx and could mostly
    // be skipped when storing coefficients.
    tridagMatrix(&avec(jy, kz, 0), &bvec(jy, kz, 0), &cvec(jy, kz, 0), &bcmplx(kz, 0), jy,
                 // wave number index
                 kz,
                 // wave number (different from kz only if we are taking a part
                 // of the z-domain [and not from 0 to 2*pi])
                 kz * kwaveFactor, global_flags, inner_boundary_flags,
                 outer_boundary_flags, &A, &C, &D);

    // Patch up internal boundaries
    if (not localmesh->lastX()) {
      for (int ix = localmesh->xend + 1; ix < localmesh->LocalNx; ix++) {
        avec(jy, kz, ix) = 0;
        bvec(jy, kz, ix) = 1;
        cvec(jy, kz, ix) = 0;
        bcmplx(kz, ix) = 0;
      }
    }
    if (not localmesh->firstX()) {
      for (int ix = 0; ix < localmesh->xstart; ix++) {
        avec(jy, kz, ix) = 0;
        bvec(jy, kz, ix) = 1;
        cvec(jy, kz, ix) = 0;
        bcmplx(kz, ix) = 0;
      }
    }
  }
  /// SCOREP_USER_REGION_END(kzinit);
  /// SCOREP_USER_REGION_DEFINE(initlevels);
  /// SCOREP_USER_REGION_BEGIN(initlevels, "init
  /// levels",///SCOREP_USER_REGION_TYPE_COMMON);

///  output.write("before\n");
///  for(int kz=0;kz<nmode;kz++){
///    for(int ix=0;ix<localmesh->LocalNx;ix++){
///      output.write("{} ",avec(jy,kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<localmesh->LocalNx;ix++){
///      output.write("{} ",bvec(jy,kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<localmesh->LocalNx;ix++){
///      output.write("{} ",cvec(jy,kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<localmesh->LocalNx;ix++){
///      output.write("{} ",bcmplx(kz,ix).real());
///    }
///    output.write("\n");
///  }

  // eliminate boundary rows - this is necessary to ensure we solve a square
  // system of interior rows
  if (localmesh->firstX()) {
    for (int kz = 0; kz < nmode; kz++) {
      bvec(jy,kz,xs) = bvec(jy,kz,xs) - cvec(jy, kz, xs-1) * avec(jy,kz,xs) / bvec(jy, kz, xs-1);
    }
  }
  if (localmesh->lastX()) {
    for (int kz = 0; kz < nmode; kz++) {
      bvec(jy,kz,xe) = bvec(jy,kz,xe) - cvec(jy, kz, xe) * avec(jy,kz,xe+1) / bvec(jy, kz, xe+1);
    }
  }

  // Perform the parallel triadiagonal solver
  // Note the API switches sub and super diagonals
  //output.write("Before solve\n");
  cr_pcr_solver(cvec,bvec,avec,bcmplx,xk1d,jy);

  //output.write("Before bcs\n");
  // apply boundary conditions
  if (localmesh->firstX()) {
    for (int kz = 0; kz < nmode; kz++) {
      for(int ix = xs-1; ix >= 0; ix--){
        xk1d(kz,ix) = (bcmplx(kz, ix) 
                      - cvec(jy, kz, ix) * xk1d(kz,ix+1)) / bvec(jy, kz, ix);
      }
    }
  }
  if (localmesh->lastX()) {
    for (int kz = 0; kz < nmode; kz++) {
      for(int ix = xe+1; ix < localmesh->LocalNx; ix++){
        xk1d(kz,ix) = (bcmplx(kz,ix) 
		      - avec(jy, kz, ix) * xk1d(kz,ix-1)) / bvec(jy, kz, ix);
      }
    }
  }

///  output.write("after\n");
///  for(int kz=0;kz<nmode;kz++){
///    for(int ix=0;ix<localmesh->LocalNx;ix++){
///      output.write("{} ",avec(jy,kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<localmesh->LocalNx;ix++){
///      output.write("{} ",bvec(jy,kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<localmesh->LocalNx;ix++){
///      output.write("{} ",cvec(jy,kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<localmesh->LocalNx;ix++){
///      output.write("{} ",bcmplx(kz,ix).real());
///    }
///    output.write("\n");
///  }

  //output.write("After bcs\n");

#if CHECK > 2
  for (int ix = 0; ix < ncx; ix++) {
    for (int kz = 0; kz < nmode; kz++) {
      if (!finite(xk1d(kz, ix).real()) or !finite(xk1d(kz, ix).imag()))
        throw BoutException("Non-finite xloc at {:d}, {:d}, {:d}", ix, jy, kz);
    }
  }
#endif

  // If the global flag is set to INVERT_KX_ZERO
  if (isGlobalFlagSet(INVERT_KX_ZERO)) {
    dcomplex offset(0.0);
    for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
      offset += xk1d(0, ix);
    }
    offset /= static_cast<BoutReal>(localmesh->xend - localmesh->xstart + 1);
    for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
      xk1d(0, ix) -= offset;
    }
  }

  // Store the solution xk for the current fourier mode in a 2D array
  transpose(xk, xk1d);
  //output.write("After transpose 2\n");

  /// SCOREP_USER_REGION_END(afterloop);

  /// SCOREP_USER_REGION_DEFINE(fftback);
  /// SCOREP_USER_REGION_BEGIN(fftback, "fft back",///SCOREP_USER_REGION_TYPE_COMMON);
  // Done inversion, transform back
  for (int ix = 0; ix < ncx; ix++) {

    if (isGlobalFlagSet(INVERT_ZERO_DC)) {
      xk(ix, 0) = 0.0;
    }

    irfft(&xk(ix, 0), ncz, x[ix]);

#if CHECK > 2
    for (int kz = 0; kz < ncz; kz++)
      if (!finite(x(ix, kz)))
        throw BoutException("Non-finite at {:d}, {:d}, {:d}", ix, jy, kz);
#endif
  }

  //output.write("end\n");

  /// SCOREP_USER_REGION_END(fftback);
  return x; // Result of the inversion
}


/** 
 * @brief   Initialize local private variables from global input parameters.
 * @param   n Size of global array
 * @param   np_world Number of MPI process
 * @param   rank_world rank ID in MPI_COMM_WORLD
*/
void LaplacePCR :: setup(int n, int np_world, int rank_world)
{
    nprocs = np_world;
    myrank = rank_world;
    n_mpi = n / nprocs;
}
/** 
 * @brief   CR-PCR solver: cr_forward_multiple + pcr_forward_single + cr_backward_multiple
 * @param   a_mpi (input) Lower off-diagonal coeff., which is assigned to local private pointer a
 * @param   b_mpi (input) Diagonal coeff., which is assigned to local private pointer b
 * @param   c_mpi (input) Upper off-diagonal coeff.,, which is assigned to local private pointer c
 * @param   r_mpi (input) RHS vector, which is assigned to local private pointer r
 * @param   x_mpi (output) Solution vector, which is assigned to local private pointer x
*/
void LaplacePCR :: cr_pcr_solver(Tensor<dcomplex> &a_mpi, Tensor<dcomplex> &b_mpi, Tensor<dcomplex> &c_mpi, Matrix<dcomplex> &r_mpi, Matrix<dcomplex> &x_mpi, int jy)
{

  //output.write("nmode {}\n",nmode);
  int nxloc = localmesh->xend-localmesh->xstart+1;
  a.reallocate(nmode, nxloc+2);
  b.reallocate(nmode, nxloc+2);
  c.reallocate(nmode, nxloc+2);
  r.reallocate(nmode, nxloc+2);
  x.reallocate(nmode, nxloc+2);

  int xs = localmesh->xstart;
  for(int kz=0; kz<nmode; kz++){
    a(kz,0) = 0;
    b(kz,0) = 1;
    c(kz,0) = 0;
    r(kz,0) = 0;
    x(kz,0) = 0;
    for(int ix=xs; ix<localmesh->xend+1; ix++){
      a(kz,ix-xs+1) = a_mpi(jy,kz,ix);
      b(kz,ix-xs+1) = b_mpi(jy,kz,ix);
      c(kz,ix-xs+1) = c_mpi(jy,kz,ix);
      r(kz,ix-xs+1) = r_mpi(kz,ix);
      x(kz,ix-xs+1) = x_mpi(kz,ix);
    }
    a(kz,nxloc+1) = 0;
    b(kz,nxloc+1) = 1;
    c(kz,nxloc+1) = 0;
    r(kz,nxloc+1) = 0;
    x(kz,nxloc+1) = 0;
  }

  output.write("data\n");
///  for(int kz=0;kz<nmode;kz++){
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",a(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",b(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",c(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",r(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",x(kz,ix).real());
///    }
///    output.write("\n");
///  }

  cr_forward_multiple_row(a,b,c,r);

  output.write("after forward\n");
///  for(int kz=0;kz<nmode;kz++){
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",a(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",b(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",c(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",r(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",x(kz,ix).real());
///    }
///    output.write("\n");
///  }
//
    pcr_forward_single_row(a,b,c,r,x);     // Including 2x2 solver

    output.write("after forward single row\n");
///  for(int kz=0;kz<nmode;kz++){
///    output.write("kz {}\n",kz);
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",a(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",b(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",c(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",r(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",x(kz,ix).real());
///    }
///    output.write("\n");
///  }
/////
    cr_backward_multiple_row(a,b,c,r,x);

  output.write("after backward multiple row\n");
///  for(int kz=0;kz<nmode;kz++){
///    output.write("kz {}\n",kz);
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",a(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",b(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",c(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",r(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",x(kz,ix).real());
///    }
///    output.write("\n");
///  }

    for(int kz=0; kz<nmode; kz++){
      for(int ix=xs; ix<localmesh->xend+1; ix++){
        x_mpi(kz,ix) = x(kz,ix-xs+1);
      }
    }
    output.write("end\n");
}

/** 
 * @brief   Forward elimination of CR until a single row per MPI process remains.
 * @details After a single row per MPI process remains, PCR or CR between a single row is performed.
*/
void LaplacePCR :: cr_forward_multiple_row(Matrix<dcomplex> &a,Matrix<dcomplex> &b,Matrix<dcomplex> &c,Matrix<dcomplex> &r)
{
    MPI_Comm comm = BoutComm::get();
    int i, l;
    int nlevel;
    int ip, in, start, dist_row, dist2_row;
    Array<dcomplex> alpha(nmode);
    Array<dcomplex> gamma(nmode);
    Array<dcomplex> sbuf(4*nmode);
    Array<dcomplex> rbuf(4*nmode);

    MPI_Status status, status1;
    Array<MPI_Request> request(2);

   //int nxloc = localmesh->xend-localmesh->xstart+1;
  output.write("start\n");
///  for(int kz=0;kz<nmode;kz++){
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",a(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",b(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",c(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",r(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",x(kz,ix).real());
///    }
///    output.write("\n");
///  }

    /// Variable nlevel is used to indicates when single row remains.
    nlevel    = log2(n_mpi);
    dist_row  = 1;
    dist2_row = 2;
    
    for(l=0;l<nlevel;l++) {
        //output.write("level {}, n_mpi {}, nlevel {}\n",l,n_mpi,nlevel);
        //output.write("myrank {}, nprocs {}\n",myrank,nprocs);
        start = dist2_row;
        /// Data exchange is performed using MPI send/recv for each succesive reduction
        if(myrank<nprocs-1) {
          //output.write("before Irecv\n");
          MPI_Irecv(&rbuf[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank+1, 0, comm, &request[0]);
          //output.write("after Irecv\n");
        }
        if(myrank>0) {
            //output.write("filling sbuf\n");
            for(int kz=0; kz<nmode; kz++){
              //output.write("myrank {}, kz {}\n",myrank, kz);
              sbuf[0+4*kz] = a(kz,dist_row);
              sbuf[1+4*kz] = b(kz,dist_row);
              sbuf[2+4*kz] = c(kz,dist_row);
              sbuf[3+4*kz] = r(kz,dist_row);
	    }
            //output.write("before isend\n");
            MPI_Isend(&sbuf[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank-1, 0, comm, &request[1]);
        }
        if(myrank<nprocs-1) {
            //output.write("before wait\n");
            MPI_Wait(&request[0], &status1);
            //output.write("after wait\n");
            for(int kz=0; kz<nmode; kz++){
	      //output.write("myrank {}, kz {}\n",myrank, kz);
              //output.write("a {}\n",a(kz,n_mpi+1).real());
              //output.write("rbuf {}\n",rbuf[0+4*kz].real());
              a(kz,n_mpi+1) = rbuf[0+4*kz];
              b(kz,n_mpi+1) = rbuf[1+4*kz];
              c(kz,n_mpi+1) = rbuf[2+4*kz];
              r(kz,n_mpi+1) = rbuf[3+4*kz];
	    }
        }

        //output.write("after sends\n");


        /// Odd rows of remained rows are reduced to even rows of remained rows in each reduction step.
        /// Index in of global last row is out of range, but we treat it as a = c = r = 0 and b = 1 in main function.
        for(i=start;i<=n_mpi;i+=dist2_row) {
            ip = i - dist_row;
            in = min(i + dist_row, n_mpi + 1);
            for(int kz=0; kz<nmode; kz++){
              alpha[kz] = -a(kz,i) / b(kz,ip);
              gamma[kz] = -c(kz,i) / b(kz,in);

              b(kz,i) += (alpha[kz] * c(kz,ip) + gamma[kz] * a(kz,in));
              a(kz,i) = alpha[kz] * a(kz,ip);
              c(kz,i) = gamma[kz] * c(kz,in);
              r(kz,i) += (alpha[kz] * r(kz,ip) + gamma[kz] * r(kz,in));
	    }

        }
        //output.write("after loop\n");
        /// As reduction continues, the indices of required coefficients doubles.
        dist2_row *= 2;
        dist_row *= 2;
        
        if(myrank>0) {
            //MPI_Wait(request+1, &status);
            MPI_Wait(&request[1], &status);
        }
    }

///  output.write("after loops\n");
///  for(int kz=0;kz<nmode;kz++){
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",a(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",b(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",c(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",r(kz,ix).real());
///    }
///    output.write("\n");
///    for(int ix=0;ix<nxloc+2;ix++){
///      output.write("{} ",x(kz,ix).real());
///    }
///    output.write("\n");
///  }
}

/** 
 * @brief   Backward substitution of CR after single-row solution per MPI process is obtained.
*/
void LaplacePCR :: cr_backward_multiple_row(Matrix<dcomplex> &a,Matrix<dcomplex> &b,Matrix<dcomplex> &c,Matrix<dcomplex> &r,Matrix<dcomplex> &x)
{
    int i, l;
    int nlevel;
    int ip, in, dist_row, dist2_row;
    MPI_Comm comm = BoutComm::get();

    MPI_Status status;
    MPI_Request request[2];
    auto recvvec = Array<dcomplex>(nmode);
    auto sendvec = Array<dcomplex>(nmode);

    nlevel    = log2(n_mpi);
    dist_row = n_mpi/2;

    /// Each rank requires a solution on last row of previous rank.
    if(myrank>0) {
	//output.write("before irecv, myrank {}\n",myrank);
        MPI_Irecv(&recvvec[0], nmode, MPI_DOUBLE_COMPLEX, myrank-1, 100, comm, request);
    }
    if(myrank<nprocs-1) {
	//output.write("before isend\n");
        for(int kz=0; kz<nmode; kz++){
	  sendvec[kz] = x(kz,n_mpi);
	}
        MPI_Isend(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, myrank+1, 100, comm, request+1);
    }
    if(myrank>0) {
	//output.write("before wait\n");
        MPI_Wait(request, &status);
        for(int kz=0; kz<nmode; kz++){
	  x(kz,0) = recvvec[kz];
	}
	//output.write("after wait\n");
    }
    for(l=nlevel-1;l>=0;l--) {
        dist2_row = dist_row * 2;
        for(i=n_mpi-dist_row;i>=0;i-=dist2_row) {
            ip = i - dist_row;
            in = i + dist_row;
	    for(int kz=0;kz<nmode; kz++){
              x(kz,i) = r(kz,i)-c(kz,i)*x(kz,in)-a(kz,i)*x(kz,ip);
              x(kz,i) = x(kz,i)/b(kz,i);
	    }
        }
        dist_row = dist_row / 2;
    }
    if(myrank<nprocs-1) {
	//output.write("before wait\n");
        MPI_Wait(request+1, &status);
	//output.write("after wait\n");
    }
    //output.write("end part 3\n");
}

/** 
 * @brief   PCR between a single row per MPI process and 2x2 matrix solver between i and i+nprocs/2 rows. 
*/
void LaplacePCR :: pcr_forward_single_row(Matrix<dcomplex> &a,Matrix<dcomplex> &b,Matrix<dcomplex> &c,Matrix<dcomplex> &r,Matrix<dcomplex> &x)
{

    int i, l, nhprocs;
    int nlevel;
    int ip, in, dist_rank, dist2_rank;
    int myrank_level, nprocs_level;
    Array<dcomplex> alpha(nmode);
    Array<dcomplex> gamma(nmode);
    Array<dcomplex> sbuf(4*nmode);
    Array<dcomplex> rbuf0(4*nmode);
    Array<dcomplex> rbuf1(4*nmode);
    dcomplex det;

    MPI_Status status;
    Array<MPI_Request> request(4);
    MPI_Comm comm = BoutComm::get();

    nlevel      = log2(nprocs);
    nhprocs     = nprocs/2;
    dist_rank   = 1;
    dist2_rank  = 2;

    /// Parallel cyclic reduction continues until 2x2 matrix are made between a pair of rank, 
    /// (myrank, myrank+nhprocs).
    for(l=0;l<nlevel-1;l++) {

        /// Rank is newly calculated in each level to find communication pair.
        /// Nprocs is also newly calculated as myrank is changed.
        myrank_level = myrank / dist_rank;
        nprocs_level = nprocs / dist_rank;

        /// All rows exchange data for reduction and perform reduction successively.
        /// Coefficients are updated for every rows.
	for(int kz=0;kz<nmode;kz++){
          sbuf[0+4*kz] = a(kz,n_mpi);
          sbuf[1+4*kz] = b(kz,n_mpi);
          sbuf[2+4*kz] = c(kz,n_mpi);
          sbuf[3+4*kz] = r(kz,n_mpi);
	}

        if((myrank_level+1)%2 == 0) {
            if(myrank+dist_rank<nprocs) {
                MPI_Irecv(&rbuf1[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank+dist_rank, 202, comm, &request[0]);
                MPI_Isend(&sbuf[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank+dist_rank, 203, comm, &request[1]);
            }
            if(myrank-dist_rank>=0) {
                MPI_Irecv(&rbuf0[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank-dist_rank, 200, comm, &request[2]);
                MPI_Isend(&sbuf[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank-dist_rank, 201, comm, &request[3]);
            }
            if(myrank+dist_rank<nprocs) {
                MPI_Wait(&request[0], &status);
	        for(int kz=0;kz<nmode;kz++){
                  a(kz,n_mpi+1) = rbuf1[0+4*kz];
                  b(kz,n_mpi+1) = rbuf1[1+4*kz];
                  c(kz,n_mpi+1) = rbuf1[2+4*kz];
                  r(kz,n_mpi+1) = rbuf1[3+4*kz];
		}
                MPI_Wait(&request[1], &status);
            }
            if(myrank-dist_rank>=0) {
                MPI_Wait(&request[2], &status);
	        for(int kz=0;kz<nmode;kz++){
                  a(kz,0) = rbuf0[0+4*kz];
                  b(kz,0) = rbuf0[1+4*kz];
                  c(kz,0) = rbuf0[2+4*kz];
                  r(kz,0) = rbuf0[3+4*kz];
		}
                MPI_Wait(&request[3], &status);
            }
        }
        else if((myrank_level+1)%2 == 1) {
            if(myrank+dist_rank<nprocs) {
                MPI_Irecv(&rbuf1[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank+dist_rank, 201, comm, &request[0]);
                MPI_Isend(&sbuf[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank+dist_rank, 200, comm, &request[1]);
            }
            if(myrank-dist_rank>=0) {
                MPI_Irecv(&rbuf0[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank-dist_rank, 203, comm, &request[2]);
                MPI_Isend(&sbuf[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank-dist_rank, 202, comm, &request[3]);
            }
            if(myrank+dist_rank<nprocs) {
                MPI_Wait(&request[0], &status);
	        for(int kz=0;kz<nmode;kz++){
                  a(kz,n_mpi+1) = rbuf1[0+4*kz];
                  b(kz,n_mpi+1) = rbuf1[1+4*kz];
                  c(kz,n_mpi+1) = rbuf1[2+4*kz];
                  r(kz,n_mpi+1) = rbuf1[3+4*kz];
		}
                MPI_Wait(&request[1], &status);
            }
            if(myrank-dist_rank>=0) {
                MPI_Wait(&request[2], &status);
	        for(int kz=0;kz<nmode;kz++){
                  a(kz,0) = rbuf0[0+4*kz];
                  b(kz,0) = rbuf0[1+4*kz];
                  c(kz,0) = rbuf0[2+4*kz];
                  r(kz,0) = rbuf0[3+4*kz];
		}
                MPI_Wait(&request[3], &status);
            }
        }

        i = n_mpi;
        ip = 0;
        in = i + 1;
        if(myrank_level == 0) {
	  for(int kz=0;kz<nmode;kz++){
            alpha[kz] = 0.0;
	  }
        }
        else {
	  for(int kz=0;kz<nmode;kz++){
            alpha[kz] = -a(kz,i) / b(kz,ip);
	  }
        }
        if(myrank_level == nprocs_level-1) {
	  for(int kz=0;kz<nmode;kz++){
            gamma[kz] = 0.0;
	  }
        }
        else {
	  for(int kz=0;kz<nmode;kz++){
            gamma[kz] = -c(kz,i) / b(kz,in);
	  }
        }

	for(int kz=0;kz<nmode;kz++){
          b(kz,i) += (alpha[kz] * c(kz,ip) + gamma[kz] * a(kz,in));
          a(kz,i)  = alpha[kz] * a(kz,ip);
          c(kz,i)  = gamma[kz] * c(kz,in);
          r(kz,i) += (alpha[kz] * r(kz,ip) + gamma[kz] * r(kz,in));
	}

        dist_rank  *= 2;
        dist2_rank *= 2;
    }

    /// Solving 2x2 matrix. All pair of ranks, myrank and myrank+nhprocs, solves it simultaneously.
    for(int kz=0;kz<nmode;kz++){
      sbuf[0+4*kz] = a(kz,n_mpi);
      sbuf[1+4*kz] = b(kz,n_mpi);
      sbuf[2+4*kz] = c(kz,n_mpi);
      sbuf[3+4*kz] = r(kz,n_mpi);
    }
    if(myrank<nhprocs) {
        MPI_Irecv(&rbuf1[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank+nhprocs, 300, comm, &request[0]);
        MPI_Isend(&sbuf[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank+nhprocs, 301, comm, &request[1]);

        MPI_Wait(&request[0], &status);
        for(int kz=0;kz<nmode;kz++){
          a(kz,n_mpi+1) = rbuf1[0+4*kz];
          b(kz,n_mpi+1) = rbuf1[1+4*kz];
          c(kz,n_mpi+1) = rbuf1[2+4*kz];
          r(kz,n_mpi+1) = rbuf1[3+4*kz];
	}

        i = n_mpi;
        in = n_mpi+1;

        for(int kz=0;kz<nmode;kz++){
          det = b(kz,i)*b(kz,in) - c(kz,i)*a(kz,in);
          x(kz,i) = (r(kz,i)*b(kz,in) - r(kz,in)*c(kz,i))/det;
          x(kz,in) = (r(kz,in)*b(kz,i) - r(kz,i)*a(kz,in))/det;
	}
        MPI_Wait(&request[1], &status);

    }
    else if(myrank>=nhprocs) {
        MPI_Irecv(&rbuf0[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank-nhprocs, 301, comm, &request[2]);
        MPI_Isend(&sbuf[0], 4*nmode, MPI_DOUBLE_COMPLEX, myrank-nhprocs, 300, comm, &request[3]);

        MPI_Wait(&request[2], &status);
        for(int kz=0;kz<nmode;kz++){
          a(kz,0) = rbuf0[0+4*kz];
          b(kz,0) = rbuf0[1+4*kz];
          c(kz,0) = rbuf0[2+4*kz];
          r(kz,0) = rbuf0[3+4*kz];
	}

        ip = 0;
        i = n_mpi;

        for(int kz=0;kz<nmode;kz++){
          det = b(kz,ip)*b(kz,i) - c(kz,ip)*a(kz,i);
          x(kz,ip) = (r(kz,ip)*b(kz,i) - r(kz,i)*c(kz,ip))/det;
          x(kz,i) = (r(kz,i)*b(kz,ip) - r(kz,ip)*a(kz,i))/det;
	}
        MPI_Wait(&request[3], &status);
    }
}
