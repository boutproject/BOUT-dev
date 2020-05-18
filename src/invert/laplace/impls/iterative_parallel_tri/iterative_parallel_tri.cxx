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

#include "iterative_parallel_tri.hxx"
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

LaplaceIPT::LaplaceIPT(Options* opt, CELL_LOC loc, Mesh* mesh_in)
    : Laplacian(opt, loc, mesh_in), A(0.0), C(1.0), D(1.0), ipt_mean_its(0.), ncalls(0) {
  A.setLocation(location);
  C.setLocation(location);
  D.setLocation(location);

  OPTION(opt, rtol, 1.e-7);
  OPTION(opt, atol, 1.e-20);
  OPTION(opt, maxits, 100);
  OPTION(opt, max_level, 3);
  OPTION(opt, max_cycle, 3);
  OPTION(opt, predict_exit, false);

  // Number of procs must be a factor of 2
  const int n = localmesh->NXPE;
  if (!is_pow2(n)) {
    throw BoutException("LaplaceIPT error: NXPE must be a power of 2");
  }
  // Number of levels cannot must be such that nproc <= 2^(max_level-1)
  if (n > 1 and n < pow(2, max_level + 1)) {
    throw BoutException("LaplaceIPT error: number of levels and processors must satisfy "
                        "NXPE > 2^(max_levels+1).");
  }

  static int ipt_solver_count = 1;
  bout::globals::dump.addRepeat(
      ipt_mean_its, "ipt_solver" + std::to_string(ipt_solver_count) + "_mean_its");
  ++ipt_solver_count;

  first_call = Array<bool>(localmesh->LocalNy);

  x0saved = Tensor<dcomplex>(localmesh->LocalNy, 4, localmesh->LocalNz / 2 + 1);

  levels = std::vector<Level>(max_level + 1);

  converged = Array<bool>(localmesh->LocalNz / 2 + 1);

  fine_error = Matrix<dcomplex>(4, localmesh->LocalNz / 2 + 1);

  resetSolver();
}

/*
 * Reset the solver to its initial state
 */
void LaplaceIPT::resetSolver() {
  first_call = true;
  x0saved = 0.0;
  resetMeanIterations();
}

/*
 * Check whether the reduced matrix on the coarsest level is diagonally
 * dominant, i.e. whether for every row the absolute value of the diagonal
 * element is greater-or-equal-to the sum of the absolute values of the other
 * elements. Being diagonally dominant is sufficient (but not necessary) for
 * the Gauss-Seidel iteration to converge.
 */
bool LaplaceIPT::Level::is_diagonally_dominant(const LaplaceIPT& l) {

  for (int kz = 0; kz < l.nmode; kz++) {
    // Check index 1 on all procs, except: the last proc only has index 1 if the
    // max_level == 0.
    if (not l.localmesh->lastX() or l.max_level == 0) {
      if (std::fabs(ar(l.jy, 1, kz)) + std::fabs(cr(l.jy, 1, kz))
          > std::fabs(br(l.jy, 1, kz))) {
        output << BoutComm::rank() << " jy=" << l.jy << ", kz=" << kz
               << ", lower row not diagonally dominant" << endl;
        return false;
      }
    }
    // Check index 2 on final proc only.
    if (l.localmesh->lastX()) {
      if (std::fabs(ar(l.jy, 2, kz)) + std::fabs(cr(l.jy, 2, kz))
          > std::fabs(br(l.jy, 2, kz))) {
        output << BoutComm::rank() << " jy=" << l.jy << ", kz=" << kz
               << ", upper row not diagonally dominant" << endl;
        return false;
      }
    }
  }
  // Have checked all modes and all are diagonally dominant
  return true;
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
 * by rearranging:
 * my_xk1d(xs) = rl(xs) + al(xs)*lower_xk1d(xe) + bl(xs)*upper_xk1d(xs)
 *     xloc(1) = rl(xs) + al(xs)*xloc(0) + bl(xs)*xloc(3)
 */
void LaplaceIPT::Level::reconstruct_full_solution(const LaplaceIPT& l,
                                                  Matrix<dcomplex>& xk1d) {
  SCOREP0();

  Array<dcomplex> x_lower(l.nmode), x_upper(l.nmode);

  for (int kz = 0; kz < l.nmode; kz++) {

    x_lower[kz] = xloc(0, kz);
    x_upper[kz] = xloc(3, kz);

    if (not l.localmesh->firstX()) {
      x_lower[kz] = (xloc(1, kz) - rl[kz] - bl(l.jy, kz) * xloc(3, kz)) / al(l.jy, kz);
    }
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    for (int i = 0; i < l.ncx; i++) {
      xk1d(kz, i) = minvb(kz, i) + upperGuardVector(l.jy, kz, i) * x_upper[kz]
                    + lowerGuardVector(l.jy, kz, i) * x_lower[kz];
    }
  }
}

// TODO Move to Array
/*
 * Returns true if all values of bool array are true, otherwise returns false.
 */
bool LaplaceIPT::all(const Array<bool> a) {
  SCOREP0();
  return std::all_of(a.begin(), a.end(), [](bool v) { return v; });
}

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
// FieldPerp LaplaceIPT::solve(const FieldPerp& b, const FieldPerp& x0, const FieldPerp&
// b0) {
FieldPerp LaplaceIPT::solve(const FieldPerp& b, const FieldPerp& x0) {

  SCOREP0();
  Timer timer("invert"); ///< Start timer

  /// SCOREP_USER_REGION_DEFINE(initvars);
  /// SCOREP_USER_REGION_BEGIN(initvars, "init vars",///SCOREP_USER_REGION_TYPE_COMMON);

  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  TRACE("LaplaceIPT::solve(const, const)");

  FieldPerp x{emptyFrom(b)};

  // Info for halo swaps
  int xproc = localmesh->getXProcIndex();
  int yproc = localmesh->getYProcIndex();
  nproc = localmesh->getNXPE();
  myproc = yproc * nproc + xproc;
  proc_in = myproc - 1;
  proc_out = myproc + 1;

  nmode = maxmode + 1;
  jy = b.getIndex();

  int ncz = localmesh->LocalNz; // Number of local z points
  ncx = localmesh->LocalNx;     // Number of local x points

  int xs = localmesh->xstart; // First interior point
  int xe = localmesh->xend;   // Last interior point

  // Calculation variables
  // proc: |       p-1       |          p          |       p+1
  //    x: |            xs-1 | xs      ...      xe | xe+1
  // xloc: |xloc[0] ...      | xloc[1] ... xloc[2] | xloc[3]
  //
  // In this method, we write the original tridiagonal system as a reduced
  // tridiagonal system whose grid points are the first interior point on each
  // processor (xs), plus the last interior point on the final processor. Each
  // processor solves for its local points. Quantities are stored in length 4
  // arrays where index 1 corresponds to the value at a processor's xs on the
  // original grid. Indices 0 and 3 correspond to the lower and upper
  // neighbours' values at their xs respectively, except on physical
  // boundaries where these are the boundary points. Index 2 is used on the
  // final processor only to track its value at the last grid point xe. This
  // means we often have special cases of the equations for final processor.
  //
  auto xloc = Matrix<dcomplex>(4, nmode); // Values of x at the processor interface

  BoutReal kwaveFactor = 2.0 * PI / coords->zlength();

  // Should we store coefficients?
  store_coefficients = not(inner_boundary_flags & INVERT_AC_GRAD);
  store_coefficients = store_coefficients && not(outer_boundary_flags & INVERT_AC_GRAD);
  store_coefficients = store_coefficients && not(inner_boundary_flags & INVERT_SET);
  store_coefficients = store_coefficients && not(outer_boundary_flags & INVERT_SET);

  // Setting the width of the boundary.
  // NOTE: The default is a width of 2 guard cells
  int inbndry = localmesh->xstart, outbndry = localmesh->xstart;

  // If the flags to assign that only one guard cell should be used is set
  if ((global_flags & INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2)) {
    inbndry = outbndry = 1;
  }
  if (inner_boundary_flags & INVERT_BNDRY_ONE)
    inbndry = 1;
  if (outer_boundary_flags & INVERT_BNDRY_ONE)
    outbndry = 1;

  /* Allocation for
  * bk   = The fourier transformed of b, where b is one of the inputs in
  *        LaplaceIPT::solve()
  * bk1d = The 1d array of bk
  * xk   = The fourier transformed of x, where x the output of
  *        LaplaceIPT::solve()
  * xk1d = The 1d array of xk
  */
  auto bk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto bk1d = Array<dcomplex>(ncx);
  auto xk = Matrix<dcomplex>(ncx, ncz / 2 + 1);
  auto xk1d = Matrix<dcomplex>(ncz / 2 + 1, ncx);

  /// SCOREP_USER_REGION_END(initvars);
  /// SCOREP_USER_REGION_DEFINE(initloop);
  /// SCOREP_USER_REGION_BEGIN(initloop, "init xk loop",SCOREP_USER_REGION_TYPE_COMMON);

  // Initialise xk to 0 as we only visit 0<= kz <= maxmode in solve
  for (int ix = 0; ix < ncx; ix++) {
    // TODO This looks wrong: should it be kz<ncz?
    // This came from serial tri
    for (int kz = maxmode + 1; kz < ncz / 2 + 1; kz++) {
      xk(ix, kz) = 0.0;
    }
  }
  /// SCOREP_USER_REGION_END(initloop);
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
  auto avec = Matrix<dcomplex>(nmode, ncx);
  auto bvec = Matrix<dcomplex>(nmode, ncx);
  auto cvec = Matrix<dcomplex>(nmode, ncx);
  auto bcmplx = Matrix<dcomplex>(nmode, ncx);

  BOUT_OMP(parallel for)
  for (int ix = 0; ix < ncx; ix++) {
    /* This for loop will set the bk (initialized by the constructor)
    * bk is the z fourier modes of b in z
    * If the INVERT_SET flag is set (meaning that x0 will be used to set the
    * bounadry values),
    */
    if (((ix < inbndry) && (inner_boundary_flags & INVERT_SET) && localmesh->firstX())
        || ((ncx - ix - 1 < outbndry) && (outer_boundary_flags & INVERT_SET)
            && localmesh->lastX())) {
      // Use the values in x0 in the boundary

      // x0 is the input
      // bk is the output
      rfft(x0[ix], ncz, &bk(ix, 0));

    } else {
      // b is the input
      // bk is the output
      rfft(b[ix], ncz, &bk(ix, 0));
      // rfft(x0[ix], ncz, &xk(ix, 0));
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
  transpose(bcmplx, bk, nmode, ncx);

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
  for (int kz = 0; kz < nmode; kz++) {
    // Note that this is called every time to deal with bcmplx and could mostly
    // be skipped when storing coefficients.
    tridagMatrix(&avec(kz, 0), &bvec(kz, 0), &cvec(kz, 0), &bcmplx(kz, 0), jy,
                 // wave number index
                 kz,
                 // wave number (different from kz only if we are taking a part
                 // of the z-domain [and not from 0 to 2*pi])
                 kz * kwaveFactor, global_flags, inner_boundary_flags,
                 outer_boundary_flags, &A, &C, &D);

    // Patch up internal boundaries
    if (not localmesh->lastX()) {
      for (int ix = localmesh->xend + 1; ix < localmesh->LocalNx; ix++) {
        avec(kz, ix) = 0;
        bvec(kz, ix) = 1;
        cvec(kz, ix) = 0;
        bcmplx(kz, ix) = 0;
      }
    }
    if (not localmesh->firstX()) {
      for (int ix = 0; ix < localmesh->xstart; ix++) {
        avec(kz, ix) = 0;
        bvec(kz, ix) = 1;
        cvec(kz, ix) = 0;
        bcmplx(kz, ix) = 0;
      }
    }
  }
  /// SCOREP_USER_REGION_END(kzinit);
  /// SCOREP_USER_REGION_DEFINE(initlevels);
  /// SCOREP_USER_REGION_BEGIN(initlevels, "init
  /// levels",///SCOREP_USER_REGION_TYPE_COMMON);

  // Initialize levels. Note that the finest grid (level 0) has a different
  // routine to coarse grids (which generally depend on the grid one step
  // finer than itself).
  //
  // If the operator to invert doesn't change from one timestep to another,
  // much of the information for each level may be stored. Data that cannot
  // be cached (e.g. the changing right-hand sides) is calculated in init_rhs
  // below.
  if (first_call[jy] || not store_coefficients) {

    levels[0].init(*this, ncx, jy, avec, bvec, cvec, xs, xe);

    int ncx_coarse = ncx; //(xe-xs+1)/2 + xs + ncx - xe - 1;
    if (max_level > 0) {
      for (int l = 1; l <= max_level; l++) {
        ncx_coarse = (ncx_coarse - 4) / 2 + 4;
        levels[l].init(*this, levels[l - 1], ncx_coarse, xs, ncx_coarse - 3, l,
                       jy); // FIXME assumes mgy=2
      }
    }
  }

  // Compute coefficients that depend on the right-hand side and which
  // therefore change every time.
  levels[0].init_rhs(*this, jy, bcmplx);

  /// SCOREP_USER_REGION_END(initlevels);

  /// SCOREP_USER_REGION_DEFINE(setsoln);
  /// SCOREP_USER_REGION_BEGIN(setsoln, "set level 0
  /// solution",///SCOREP_USER_REGION_TYPE_COMMON);
  // Set initial values with cached values
  for (int ix = 0; ix < 4; ix++) {
    for (int kz = 0; kz < nmode; kz++) {
      levels[0].xloc(ix, kz) = x0saved(jy, ix, kz);
    }
  }
  /// SCOREP_USER_REGION_END(setsoln);
  /// SCOREP_USER_REGION_DEFINE(initwhileloop);
  /// SCOREP_USER_REGION_BEGIN(initwhileloop, "init while
  /// loop",///SCOREP_USER_REGION_TYPE_COMMON);

  int count = 0;      // Total iteration count
  int subcount = 0;   // Count of iterations on a level
  int cyclecount = 0; // Number of multigrid cycles
  int cycle_eta = 0;  // Predicted finishing cycle
  int current_level = 0;
  bool down = true;

  auto error_abs = Array<BoutReal>(nmode);
  auto error_abs_old = Array<BoutReal>(nmode);
  auto error_rel = Array<BoutReal>(nmode);
  auto error_rel_old = Array<BoutReal>(nmode);
  for (int kz = 0; kz < nmode; kz++) {
    converged[kz] = false;
    error_abs[kz] = 1e20;
    error_abs_old[kz] = 1e20;
    error_rel[kz] = 1e20;
    error_rel_old[kz] = 1e20;
  }

  /// SCOREP_USER_REGION_END(initwhileloop);
  /// SCOREP_USER_REGION_DEFINE(whileloop);
  /// SCOREP_USER_REGION_BEGIN(whileloop, "while loop",///SCOREP_USER_REGION_TYPE_COMMON);

  while (true) {

    levels[current_level].gauss_seidel_red_black(*this);

    /// SCOREP_USER_REGION_DEFINE(l0rescalc);
    /// SCOREP_USER_REGION_BEGIN(l0rescalc, "level 0 residual
    /// calculation",///SCOREP_USER_REGION_TYPE_COMMON);
    if (current_level == 0 and subcount == max_cycle - 1) {

      ++cyclecount;

      // The allreduce in calculate_total_residual is expensive at scale. To
      // minimize calls to this, we estimate when the algorithm will converge
      // and don't check for convergence until we get near this point.
      //
      // Need to do call calculate_residual every time, but everything else can
      // be called only on the first and second cycle (to set up) and after the
      // predicted number of iterations has elapsed.
      //
      // This approach is optional; whether it helps depends on the
      // particular problem. The approach saves a lot of time in allreduces,
      // but since it does not check for convergence, it does a lot of extra
      // work that is masked for modes that are known to have converged.

      // Keep the old error values if they are used in the
      // calculations in this cycle.
      if (cyclecount < 3 or cyclecount > cycle_eta) {
        for (int kz = 0; kz < nmode; kz++) {
          if (!converged[kz]) {
            error_abs_old[kz] = error_abs[kz];
            error_rel_old[kz] = error_rel[kz];
          }
        }
      }

      levels[current_level].calculate_residual(*this);

      if (cyclecount < 3 or cyclecount > cycle_eta - 5 or not predict_exit) {
        // Calculate the total residual. This also marks modes as converged, so the
        // algorithm cannot exit in cycles where this is not called.
        levels[current_level].calculate_total_residual(*this, error_abs, error_rel,
                                                       converged);

        // Based the error reduction per V-cycle, error_xxx/error_xxx_old,
        // predict when the slowest converging mode converges.
        if (cyclecount < 3 and predict_exit) {
          BoutReal ratio;
          int eta;
          cycle_eta = 0;
          for (int kz = 0; kz < nmode; kz++) {
            ratio = error_abs[kz] / error_abs_old[kz];
            eta = std::ceil(std::log(atol / error_abs[kz]) / std::log(ratio));
            cycle_eta = (cycle_eta > eta) ? cycle_eta : eta;

            ratio = error_rel[kz] / error_rel_old[kz];
            eta = std::ceil(std::log(rtol / error_rel[kz]) / std::log(ratio));
            cycle_eta = (cycle_eta > eta) ? cycle_eta : eta;
          }
        }
      }
    }

    /// SCOREP_USER_REGION_END(l0rescalc);
    /// SCOREP_USER_REGION_DEFINE(increment);
    /// SCOREP_USER_REGION_BEGIN(increment, "increment
    /// counters",///SCOREP_USER_REGION_TYPE_COMMON);
    ++count;
    ++subcount;
    /// SCOREP_USER_REGION_END(increment);

    // Force at least max_cycle iterations at each level
    // Do not skip with tolerence to minimize comms
    if (subcount < max_cycle) {
    } else if (all(converged) and current_level == 0) {
      break;
    } else if (not down) {
      refine(levels[current_level], levels[current_level - 1], fine_error, converged);
      --current_level;
      levels[current_level].update_solution(*this);
      levels[current_level].synchronize_reduced_field(*this, levels[current_level].xloc);

      subcount = 0;

      if (current_level == 0) {
        down = true;
      }
    } else if (down && max_level > 0) {

      if (current_level != 0) {
        // Prevents double call on level 0 - we just called this to check convergence
        levels[current_level].calculate_residual(*this);
      }
      levels[current_level].synchronize_reduced_field(*this,
                                                      levels[current_level].residual);
      ++current_level;
      levels[current_level].coarsen(*this, levels[current_level - 1].residual);
      subcount = 0;

      // If we are on the coarsest grid, stop trying to coarsen further
      if (current_level == max_level) {
        down = false;
      }
    } else {
      // When only using one level, need to ensure subcount < max_count
      subcount = 0;
    }

    // Throw error if we are performing too many iterations
    if (count > maxits) {
      // Maximum number of allowed iterations reached.
      // If the coarsest multigrid iteration matrix is diagonally-dominant,
      // then convergence is guaranteed, so maxits is set too low.
      // Otherwise, the method may or may not converge.
      bool is_dd = levels[max_level].is_diagonally_dominant(*this);

      bool global_is_dd;
      MPI_Allreduce(&is_dd, &global_is_dd, 1, MPI::BOOL, MPI_LAND, BoutComm::get());

      if (global_is_dd) {
        throw BoutException("LaplaceIPT error: Not converged within maxits=%i "
                            "iterations. The coarsest grained iteration matrix is "
                            "diagonally dominant and convergence is guaranteed. Please "
                            "increase maxits and retry.",
                            maxits);
      } else {
        throw BoutException(
            "LaplaceIPT error: Not converged within maxits=%i iterations. The coarsest "
            "iteration matrix is not diagonally dominant so there is no guarantee this "
            "method will converge. Consider (1) increasing maxits; or (2) increasing the "
            "number of levels (as grids become more diagonally dominant with "
            "coarsening). Using more grids may require larger NXPE.",
            maxits);
      }
    }
  }
/// SCOREP_USER_REGION_END(whileloop);
/// SCOREP_USER_REGION_DEFINE(afterloop);
/// SCOREP_USER_REGION_BEGIN(afterloop, "after faff",///SCOREP_USER_REGION_TYPE_COMMON);

#if CHECK > 2
  for (int ix = 0; ix < 4; ix++) {
    for (int kz = 0; kz < nmode; kz++) {
      if (!finite(levels[0].xloc(ix, kz).real())
          or !finite(levels[0].xloc(ix, kz).imag()))
        throw BoutException("Non-finite xloc at %d, %d, %d", ix, jy, kz);
    }
  }
#endif

  // Cache solution
  for (int ix = 0; ix < 4; ix++) {
    for (int kz = 0; kz < nmode; kz++) {
      x0saved(jy, ix, kz) = levels[0].xloc(ix, kz);
    }
  }

  levels[0].reconstruct_full_solution(*this, xk1d);

#if CHECK > 2
  for (int ix = 0; ix < ncx; ix++) {
    for (int kz = 0; kz < nmode; kz++) {
      if (!finite(xk1d(kz, ix).real()) or !finite(xk1d(kz, ix).imag()))
        throw BoutException("Non-finite xloc at %d, %d, %d", ix, jy, kz);
    }
  }
#endif

  ++ncalls;
  ipt_mean_its =
      (ipt_mean_its * BoutReal(ncalls - 1) + BoutReal(count)) / BoutReal(ncalls);

  // If the global flag is set to INVERT_KX_ZERO
  if (global_flags & INVERT_KX_ZERO) {
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
  transpose(xk, xk1d, ncx, nmode);

  /// SCOREP_USER_REGION_END(afterloop);

  /// SCOREP_USER_REGION_DEFINE(fftback);
  /// SCOREP_USER_REGION_BEGIN(fftback, "fft back",///SCOREP_USER_REGION_TYPE_COMMON);
  // Done inversion, transform back
  for (int ix = 0; ix < ncx; ix++) {

    if (global_flags & INVERT_ZERO_DC)
      xk(ix, 0) = 0.0;

    irfft(&xk(ix, 0), ncz, x[ix]);

#if CHECK > 2
    for (int kz = 0; kz < ncz; kz++)
      if (!finite(x(ix, kz)))
        throw BoutException("Non-finite at %d, %d, %d", ix, jy, kz);
#endif
  }

  first_call[jy] = false;

  /// SCOREP_USER_REGION_END(fftback);
  return x; // Result of the inversion
}

/*
 * Perform Gauss-Seidel with red-black colouring on the reduced system.
 * We don't attempt comm/comp overlap, as there is not sigificant work in the
 * x loop.
 */
void LaplaceIPT::Level::gauss_seidel_red_black(const LaplaceIPT& l) {

  SCOREP0();

  if (not included) {
    return;
  }

  Array<dcomplex> sendvec(l.nmode), recvecin(l.nmode), recvecout(l.nmode);
  MPI_Request rreqin, rreqout;

  // Processor colouring. There are p = 2^m processors, labelled 0 to p-1.
  // On level 0, even procs are coloured red, odd procs are coloured black.
  // The last proc, p-1, is coloured black, but also has a "red" point, the
  // final interior point. It does not need to receive during the black
  // sweep, as the red point only needs data from the first interior point
  // which is also local.

  // BLACK SWEEP
  //
  // Red processors communication only
  if (red) {
    // Post receives
    if (not l.localmesh->firstX()) {
      MPI_Irecv(&recvecin[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, BoutComm::get(),
                &rreqin);
    }
    if (not l.localmesh->lastX()) {
      MPI_Irecv(&recvecout[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, BoutComm::get(),
                &rreqout);
    }

    // Receive and put data in arrays
    if (!l.localmesh->firstX()) {
      MPI_Wait(&rreqin, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
          xloc(0, kz) = recvecin[kz];
        }
      }
    }
    if (!l.localmesh->lastX()) {
      MPI_Wait(&rreqout, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
          xloc(3, kz) = recvecout[kz];
        }
      }
    }
  }

  // Black processors: work and communication
  if (black) {
    // Black processors do work
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        // Due to extra point on final proc, indexing of last term is 2, not 3. To
        // remove branching, this is handled by l.index_end
        xloc(1, kz) = (rr(1, kz) - ar(l.jy, 1, kz) * xloc(0, kz)
                       - cr(l.jy, 1, kz) * xloc(index_end, kz))
                      * brinv(l.jy, 1, kz);
      }
    }
    // Send same data up and down
    MPI_Send(&xloc(1, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1,
             BoutComm::get()); // black never firstX
    if (not l.localmesh->lastX()) {
      MPI_Send(&xloc(1, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, BoutComm::get());
    }
  }

  // RED SWEEP
  //
  // Black processors only comms
  if (black) {
    // Post receives
    MPI_Irecv(&recvecin[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, BoutComm::get(),
              &rreqin); // black never first
    if (not l.localmesh
                ->lastX()) { // this is always be true is we force an even core count
      MPI_Irecv(&recvecout[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, BoutComm::get(),
                &rreqout);
    }

    // Receive and put data in arrays
    MPI_Wait(&rreqin, MPI_STATUS_IGNORE);
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        xloc(0, kz) = recvecin[kz];
      }
    }
    if (!l.localmesh->lastX()) {
      MPI_Wait(&rreqout, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
          xloc(3, kz) = recvecout[kz];
        }
      }
    }
  }

  // Red processors do work and comms
  if (red and not l.localmesh->lastX()) {
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        xloc(1, kz) =
            (rr(1, kz) - ar(l.jy, 1, kz) * xloc(0, kz) - cr(l.jy, 1, kz) * xloc(3, kz))
            * brinv(l.jy, 1, kz);
      }
    }
  }
  if (l.localmesh->lastX()) {
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        // index_start removes branches. On level 0, this is 1, otherwise 0
        xloc(2, kz) = (rr(2, kz) - ar(l.jy, 2, kz) * xloc(index_start, kz)
                       - cr(l.jy, 2, kz) * xloc(3, kz))
                      * brinv(l.jy, 2, kz);
      }
    }
  }

  if (red or l.localmesh->lastX()) { // red, or last proc when not on level zero
    // Send same data up and down
    if (not l.localmesh->firstX() and not l.localmesh->lastX()) { // excludes last proc
      MPI_Send(&xloc(1, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, BoutComm::get());
    } else if (l.localmesh->lastX() and current_level != 0) { // last proc on level > 0
      MPI_Send(&xloc(2, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, BoutComm::get());
    }
    if (not l.localmesh->lastX()) {
      MPI_Send(&xloc(1, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, BoutComm::get());
    }
  }

  if (current_level == 0) {
    // Update boundaries to match interior points
    // Do this after communication
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        if (l.localmesh->firstX()) {
          xloc(0, kz) = -cvec(l.jy, kz, xs - 1) * xloc(1, kz) / bvec(l.jy, kz, xs - 1);
        }
        if (l.localmesh->lastX()) {
          xloc(3, kz) = -avec(l.jy, kz, xe + 1) * xloc(2, kz) / bvec(l.jy, kz, xe + 1);
        }
      }
    }
  }
}

// Initialization routine for coarser grids. Initialization depends on the grid
// one step finer, lup.
void LaplaceIPT::Level::init(LaplaceIPT& l, const Level lup, int ncx_in, const int xs_in,
                             const int xe_in, const int current_level_in, const int jy) {

  SCOREP0();
  xs = xs_in;
  xe = xe_in;
  ncx = ncx_in;
  int ny = l.localmesh->LocalNy;
  current_level = current_level_in;

  auto sendvec = Array<dcomplex>(3 * l.nmode);
  auto recvecin = Array<dcomplex>(3 * l.nmode);
  auto recvecout = Array<dcomplex>(3 * l.nmode);

  // indexing to remove branches in tight loops
  if (l.localmesh->lastX()) {
    index_end = 2;
  } else {
    index_end = 3;
  }
  index_start = 0;

  myproc = lup.myproc;
  // Whether this proc is involved in the multigrid calculation
  included = (myproc % int((pow(2, current_level))) == 0) or l.localmesh->lastX();
  // Whether this proc is involved in the calculation on the grid one level more refined
  included_up = lup.included;

  if (not included) {
    return;
  }

  // Colouring of processor for Gauss-Seidel
  red = ((myproc / int((pow(2, current_level)))) % 2 == 0);
  black = ((myproc / int((pow(2, current_level)))) % 2 == 1);

  // The last processor is a special case. It is always included because of
  // the final grid point, which is treated explicitly. Otherwise it should
  // not be included in either the red or black work.
  if (l.localmesh->lastX()) {
    red = true;
    black = false;
  }

  // My neighbouring procs
  proc_in = myproc - int(pow(2, current_level));
  if (l.localmesh->lastX()) {
    proc_in += 1;
  }
  int p = myproc + int(pow(2, current_level));
  proc_out = (p < l.nproc - 1) ? p : l.nproc - 1;

  residual = Matrix<dcomplex>(4, l.nmode);
  xloc = Matrix<dcomplex>(4, l.nmode);
  rl = Array<dcomplex>(l.nmode);
  ru = Array<dcomplex>(l.nmode);

  // Coefficients for the reduced iterations
  ar = Tensor<dcomplex>(ny, 4, l.nmode);
  br = Tensor<dcomplex>(ny, 4, l.nmode);
  cr = Tensor<dcomplex>(ny, 4, l.nmode);
  rr = Matrix<dcomplex>(4, l.nmode);
  brinv = Tensor<dcomplex>(ny, 4, l.nmode);

  for (int kz = 0; kz < l.nmode; kz++) {
    if (l.localmesh->firstX()) {
      ar(l.jy, 1, kz) = 0.5 * lup.ar(l.jy, 1, kz);
      br(l.jy, 1, kz) = 0.5 * lup.br(l.jy, 1, kz) + 0.25 * lup.cr(l.jy, 1, kz)
                        + 0.25 * lup.ar(l.jy, 3, kz) + 0.125 * lup.br(l.jy, 3, kz);
      cr(l.jy, 1, kz) = 0.25 * lup.cr(l.jy, 1, kz) + 0.125 * lup.br(l.jy, 3, kz)
                        + 0.25 * lup.cr(l.jy, 3, kz);
    } else {
      ar(l.jy, 1, kz) = 0.25 * lup.ar(l.jy, 0, kz) + 0.125 * lup.br(l.jy, 0, kz)
                        + 0.25 * lup.ar(l.jy, 1, kz);
      br(l.jy, 1, kz) = 0.125 * lup.br(l.jy, 0, kz) + 0.25 * lup.cr(l.jy, 0, kz)
                        + 0.25 * lup.ar(l.jy, 1, kz) + 0.5 * lup.br(l.jy, 1, kz)
                        + 0.25 * lup.cr(l.jy, 1, kz) + 0.25 * lup.ar(l.jy, 3, kz)
                        + 0.125 * lup.br(l.jy, 3, kz);
      cr(l.jy, 1, kz) = 0.25 * lup.cr(l.jy, 1, kz) + 0.125 * lup.br(l.jy, 3, kz)
                        + 0.25 * lup.cr(l.jy, 3, kz);
    }

    // Last proc does calculation on index 2 as well as index 1.
    // If current_level=1, the point to my left on the level above it my
    // index 1. Otherwise, it is my index 0.
    if (l.localmesh->lastX()) {
      if (current_level == 1) {
        ar(l.jy, 2, kz) = 0.25 * lup.ar(l.jy, 1, kz) + 0.125 * lup.br(l.jy, 1, kz)
                          + 0.25 * lup.ar(l.jy, 2, kz);
        br(l.jy, 2, kz) = 0.125 * lup.br(l.jy, 1, kz) + 0.25 * lup.cr(l.jy, 1, kz)
                          + 0.25 * lup.ar(l.jy, 2, kz) + 0.5 * lup.br(l.jy, 2, kz);
        cr(l.jy, 2, kz) = 0.5 * lup.cr(l.jy, 2, kz);
      } else {
        ar(l.jy, 2, kz) = 0.25 * lup.ar(l.jy, 0, kz) + 0.125 * lup.br(l.jy, 0, kz)
                          + 0.25 * lup.ar(l.jy, 2, kz);
        br(l.jy, 2, kz) = 0.125 * lup.br(l.jy, 0, kz) + 0.25 * lup.cr(l.jy, 0, kz)
                          + 0.25 * lup.ar(l.jy, 2, kz) + 0.5 * lup.br(l.jy, 2, kz);
        cr(l.jy, 2, kz) = 0.5 * lup.cr(l.jy, 2, kz);
      }
    }
    brinv(l.jy, 1, kz) = 1.0 / br(l.jy, 1, kz);
    brinv(l.jy, 2, kz) = 1.0 / br(l.jy, 2, kz);

    // Need to communicate my index 1 to this level's neighbours
    // Index 2 if last proc.
    if (not l.localmesh->lastX()) {
      sendvec[kz] = ar(l.jy, 1, kz);
      sendvec[kz + l.nmode] = br(l.jy, 1, kz);
      sendvec[kz + 2 * l.nmode] = cr(l.jy, 1, kz);
    } else {
      sendvec[kz] = ar(l.jy, 2, kz);
      sendvec[kz + l.nmode] = br(l.jy, 2, kz);
      sendvec[kz + 2 * l.nmode] = cr(l.jy, 2, kz);
    }
  }

  MPI_Comm comm = BoutComm::get();

  // Communicate in
  if (not l.localmesh->firstX()) {
    MPI_Sendrecv(&sendvec[0], 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, &recvecin[0],
                 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, comm, MPI_STATUS_IGNORE);
  }

  // Communicate out
  if (not l.localmesh->lastX()) {
    MPI_Sendrecv(&sendvec[0], 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, &recvecout[0],
                 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, comm, MPI_STATUS_IGNORE);
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.localmesh->firstX()) {
      ar(l.jy, 0, kz) = recvecin[kz];
      br(l.jy, 0, kz) = recvecin[kz + l.nmode];
      cr(l.jy, 0, kz) = recvecin[kz + 2 * l.nmode];
    }
    if (not l.localmesh->lastX()) {
      ar(l.jy, 3, kz) = recvecout[kz];
      br(l.jy, 3, kz) = recvecout[kz + l.nmode];
      cr(l.jy, 3, kz) = recvecout[kz + 2 * l.nmode];
    }
  }
}

// Init routine for finest level
void LaplaceIPT::Level::init(LaplaceIPT& l, const int ncx_in, const int jy,
                             const Matrix<dcomplex> avec_in,
                             const Matrix<dcomplex> bvec_in,
                             const Matrix<dcomplex> cvec_in, const int xs_in,
                             const int xe_in) {

  // Basic definitions for conventional multigrid
  SCOREP0();
  xs = xs_in;
  xe = xe_in;
  ncx = ncx_in;
  int ny = l.localmesh->LocalNy;
  current_level = 0;

  // Processor information
  myproc = l.myproc;     // unique id
  proc_in = myproc - 1;  // in-neighbour
  proc_out = myproc + 1; // out-neighbour
  included = true;       // whether processor is included in this level's calculation
  // Colouring of processor for Gauss-Seidel
  red = (myproc % 2 == 0);
  black = (myproc % 2 == 1);

  // indexing to remove branching in tight loops
  if (l.localmesh->lastX()) {
    index_end = 2;
  } else {
    index_end = 3;
  }
  index_start = 1;

  avec = Tensor<dcomplex>(ny, l.nmode, ncx);
  bvec = Tensor<dcomplex>(ny, l.nmode, ncx);
  cvec = Tensor<dcomplex>(ny, l.nmode, ncx);

  // Coefficients for the reduced iterations
  ar = Tensor<dcomplex>(ny, 4, l.nmode);
  br = Tensor<dcomplex>(ny, 4, l.nmode);
  cr = Tensor<dcomplex>(ny, 4, l.nmode);
  rr = Matrix<dcomplex>(4, l.nmode);
  brinv = Tensor<dcomplex>(ny, 4, l.nmode);

  residual = Matrix<dcomplex>(4, l.nmode);

  for (int kz = 0; kz < l.nmode; kz++) {
    for (int ix = 0; ix < ncx; ix++) {
      avec(jy, kz, ix) = avec_in(kz, ix);
      bvec(jy, kz, ix) = bvec_in(kz, ix);
      cvec(jy, kz, ix) = cvec_in(kz, ix);
    }
    for (int ix = 0; ix < 4; ix++) {
      residual(ix, kz) = 0.0;
    }
  }
  // end basic definitions

  // Define sizes of local coefficients
  xloc = Matrix<dcomplex>(4, l.nmode); // Reduced grid x values

  // Arrays to construct global solution from halo values
  minvb = Matrix<dcomplex>(l.nmode, ncx);                // Local M^{-1} f
  lowerGuardVector = Tensor<dcomplex>(ny, l.nmode, ncx); // alpha
  upperGuardVector = Tensor<dcomplex>(ny, l.nmode, ncx); // beta

  // Coefficients of first and last interior rows
  al = Matrix<dcomplex>(ny, l.nmode); // alpha^l
  bl = Matrix<dcomplex>(ny, l.nmode); // beta^l
  au = Matrix<dcomplex>(ny, l.nmode); // alpha^u
  bu = Matrix<dcomplex>(ny, l.nmode); // beta^u
  rl = Array<dcomplex>(l.nmode);      // r^l
  ru = Array<dcomplex>(l.nmode);      // r^u

  // Coefs used to compute rl from domain below
  r1 = Matrix<dcomplex>(ny, l.nmode);
  r2 = Matrix<dcomplex>(ny, l.nmode);

  // Work arrays
  auto evec = Array<dcomplex>(ncx);
  auto tmp = Array<dcomplex>(ncx);

  // Communication arrays
  auto sendvec = Array<dcomplex>(3 * l.nmode);
  auto recvecin = Array<dcomplex>(3 * l.nmode);
  auto recvecout = Array<dcomplex>(3 * l.nmode);

  for (int kz = 0; kz < l.nmode; kz++) {

    /// SCOREP_USER_REGION_DEFINE(invert);
    /// SCOREP_USER_REGION_BEGIN(invert, "invert local
    /// matrices",///SCOREP_USER_REGION_TYPE_COMMON);

    // Invert local matrices to find upper/lower guard vectos.
    // Note Minv*b is calculated in init_rhs.
    //
    // Upper interface
    if (not l.localmesh->lastX()) {
      // Need the xend-th element
      for (int i = 0; i < ncx; i++) {
        evec[i] = 0.0;
      }
      evec[xe + 1] = 1.0;
      tridag(&avec_in(kz, 0), &bvec_in(kz, 0), &cvec_in(kz, 0), std::begin(evec),
             std::begin(tmp), ncx);
      for (int i = 0; i < ncx; i++) {
        upperGuardVector(jy, kz, i) = tmp[i];
      }
    } else {
      for (int i = 0; i < ncx; i++) {
        upperGuardVector(jy, kz, i) = 0.0;
      }
    }

    // Lower interface
    if (not l.localmesh->firstX()) {
      for (int i = 0; i < ncx; i++) {
        evec[i] = 0.0;
      }
      evec[xs - 1] = 1.0;
      tridag(&avec_in(kz, 0), &bvec_in(kz, 0), &cvec_in(kz, 0), std::begin(evec),
             std::begin(tmp), ncx);
      for (int i = 0; i < ncx; i++) {
        lowerGuardVector(jy, kz, i) = tmp[i];
      }
    } else {
      for (int i = 0; i < ncx; i++) {
        lowerGuardVector(jy, kz, i) = 0.0;
      }
    }

    /// SCOREP_USER_REGION_END(invert);
    /// SCOREP_USER_REGION_DEFINE(coefs);
    /// SCOREP_USER_REGION_BEGIN(coefs, "calculate
    /// coefs",///SCOREP_USER_REGION_TYPE_COMMON);

    bl(jy, kz) = upperGuardVector(jy, kz, xs);
    al(jy, kz) = lowerGuardVector(jy, kz, xs);

    bu(jy, kz) = upperGuardVector(jy, kz, xe);
    au(jy, kz) = lowerGuardVector(jy, kz, xe);

    // First compute coefficients that depend on the matrix to be inverted
    // and which therefore might be constant throughout a run.

    // Boundary processor values to be overwritten when relevant
    MPI_Request req;
    auto AdBd = Array<dcomplex>(2);  // A and B coefficients from proc down
    auto ABtmp = Array<dcomplex>(2); // Send array for A and B coefs from proc down
    AdBd[0] = 1.0;
    AdBd[1] = 0.0;
    if (not l.localmesh->firstX()) {
      MPI_Irecv(&AdBd[0], 2, MPI_DOUBLE_COMPLEX, proc_in, 0, BoutComm::get(), &req);
    }
    if (not l.localmesh->lastX()) {
      // Send coefficients up
      ABtmp[0] = 0.0;
      ABtmp[1] = bu(jy, kz);
      if (std::fabs(al(jy, kz)) > 1e-14) {
        ABtmp[0] = au(jy, kz) / al(jy, kz);
        ABtmp[1] -= ABtmp[0] * bl(jy, kz);
      }
      // Send these
      MPI_Send(&ABtmp[0], 2, MPI_DOUBLE_COMPLEX, proc_out, 0, BoutComm::get());
    }

    if (not l.localmesh->firstX()) {
      MPI_Wait(&req, MPI_STATUS_IGNORE);
    }

    const dcomplex Delta = 1.0 / (1.0 - al(jy, kz) * AdBd[1]);
    ar(jy, 1, kz) = -Delta * al(jy, kz) * AdBd[0];
    cr(jy, 1, kz) = -Delta * bl(jy, kz);

    r1(jy, kz) = Delta * al(jy, kz);
    r2(jy, kz) = Delta;

    // lastX is a special case having two points on the level 0 grid
    if (l.localmesh->lastX()) {
      // Note that if the first proc is also the last proc, then both alold and
      // auold are zero, and l.au = l.auold is already correct.
      if (not l.localmesh->lastX()) {
        ar(jy, 2, kz) = -au(jy, kz) / al(jy, kz);
        cr(jy, 2, kz) =
            -(bu(jy, kz) + ar(jy, 2, kz) * bl(jy, kz)); // NB depends on previous line
      }

      // Use BCs to replace x(xe+1) = -avec(xe+1) x(xe) / bvec(xe+1)
      //  => only bl changes
      cr(jy, 1, kz) = avec(jy, kz, xe + 1) * bl(jy, kz) / bvec(jy, kz, xe + 1);
    }

    // Now set coefficients for reduced iterations (shared by all levels)
    br(jy, 1, kz) = 1.0;
    br(jy, 2, kz) = 1.0;
    brinv(jy, 1, kz) = 1.0;
    brinv(jy, 2, kz) = 1.0;

    sendvec[kz] = ar(jy, 1, kz);
    sendvec[kz + l.nmode] = br(jy, 1, kz);
    sendvec[kz + 2 * l.nmode] = cr(jy, 1, kz);

    /// SCOREP_USER_REGION_END(coefs);
  } // end of kz loop

  // Synchonize reduced coefficients with neighbours
  MPI_Comm comm = BoutComm::get();

  // Communicate in
  if (not l.localmesh->firstX()) {
    MPI_Sendrecv(&sendvec[0], 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, &recvecin[0],
                 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, comm, MPI_STATUS_IGNORE);
  }

  // Communicate out
  if (not l.localmesh->lastX()) {
    MPI_Sendrecv(&sendvec[0], 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, &recvecout[0],
                 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, comm, MPI_STATUS_IGNORE);
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.localmesh->firstX()) {
      ar(jy, 0, kz) = recvecin[kz];
      br(jy, 0, kz) = recvecin[kz + l.nmode];
      cr(jy, 0, kz) = recvecin[kz + 2 * l.nmode];
    }
    if (not l.localmesh->lastX()) {
      ar(jy, 3, kz) = recvecout[kz];
      br(jy, 3, kz) = recvecout[kz + l.nmode];
      cr(jy, 3, kz) = recvecout[kz + 2 * l.nmode];
    }
  }
}

// Init routine for finest level information that cannot be cached
void LaplaceIPT::Level::init_rhs(LaplaceIPT& l, const int jy,
                                 const Matrix<dcomplex> bcmplx) {

  SCOREP0();

  auto Rd = Array<dcomplex>(l.nmode);
  auto Rsendup = Array<dcomplex>(l.nmode);
  MPI_Request req;

  for (int kz = 0; kz < l.nmode; kz++) {

    /// SCOREP_USER_REGION_DEFINE(invertforrhs);
    /// SCOREP_USER_REGION_BEGIN(invertforrhs, "invert local matrices for
    /// rhs",///SCOREP_USER_REGION_TYPE_COMMON);

    // Invert local matrices
    // Calculate Minv*b
    tridag(&avec(l.jy, kz, 0), &bvec(l.jy, kz, 0), &cvec(l.jy, kz, 0), &bcmplx(kz, 0),
           &minvb(kz, 0), l.ncx);
    // Now minvb is a constant vector throughout the iterations

    /// SCOREP_USER_REGION_END(invertforrhs);
    /// SCOREP_USER_REGION_DEFINE(coefsforrhs);
    /// SCOREP_USER_REGION_BEGIN(coefsforrhs, "calculate coefs for
    /// rhs",///SCOREP_USER_REGION_TYPE_COMMON);

    rl[kz] = minvb(kz, xs);
    ru[kz] = minvb(kz, xe);

    // Boundary processor value to be overwritten when relevant
    Rd[kz] = 0.0;

    if (not l.localmesh->lastX()) {
      // Send coefficients up
      Rsendup[kz] = ru[kz];
      if (std::fabs(al(l.jy, kz)) > 1e-14) {
        Rsendup[kz] -= rl[kz] * au(jy, kz) / al(jy, kz);
      }
    }
    /// SCOREP_USER_REGION_END(coefsforrhs);
  } // end of kz loop

  if (not l.localmesh->firstX()) {
    MPI_Irecv(&Rd[0], l.nmode, MPI_DOUBLE_COMPLEX, l.proc_in, 0, BoutComm::get(), &req);
  }
  if (not l.localmesh->lastX()) {
    MPI_Send(&Rsendup[0], l.nmode, MPI_DOUBLE_COMPLEX, l.proc_out, 0, BoutComm::get());
  }
  if (not l.localmesh->firstX()) {
    MPI_Wait(&req, MPI_STATUS_IGNORE);
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    rr(1, kz) = r1(l.jy, kz) * Rd[kz] + r2(l.jy, kz) * rl[kz];

    // Special case for multiple points on last proc
    // Note that if the first proc is also the last proc, then both al and
    // au are zero, and l.ru is already correct.
    if (l.localmesh->lastX() and not l.localmesh->firstX()) {
      rr(2, kz) = ru[kz] - au(l.jy, kz) * rl[kz] / al(l.jy, kz);
    } else {
      rr(2, kz) = ru[kz];
    }
  }
}

/*
 * Sum and communicate total residual for the reduced system
 */
void LaplaceIPT::Level::calculate_total_residual(LaplaceIPT& l,
                                                 Array<BoutReal>& error_abs,
                                                 Array<BoutReal>& error_rel,
                                                 Array<bool>& converged) {

  SCOREP0();
  // Communication arrays:
  // residual in (0 .. nmode-1)
  // solution in (nmode .. 2*nmode-1)
  auto total = Array<BoutReal>(2 * l.nmode);    // global summed residual
  auto subtotal = Array<BoutReal>(2 * l.nmode); // local contribution to residual

  for (int kz = 0; kz < l.nmode; kz++) {
    if (!converged[kz]) {
      total[kz] = 0.0;
      total[kz + l.nmode] = 0.0;

      // TODO add special case for last proc
      // Only xs and xe have nonzero residuals
      subtotal[kz] = pow(residual(1, kz).real(), 2) + pow(residual(1, kz).imag(), 2)
                     + pow(residual(2, kz).real(), 2)
                     + pow(residual(2, kz).imag(), 2);

      // TODO This approximation will increase iteration count. The alternatives are:
      // + reconstructing the solution and calculating properly
      // + multiply approximation by (interior points/2) - this can be done
      //   at runtime by changing rtol
      // Strictly this should be all contributions to the solution, but this
      // under-approximation saves work.
      subtotal[kz + l.nmode] = pow(xloc(1, kz).real(), 2) + pow(xloc(1, kz).imag(), 2)
                               + pow(xloc(2, kz).real(), 2) + pow(xloc(2, kz).imag(), 2);
    }
  }

  // Communication needed to ensure processors break on same iteration
  MPI_Allreduce(&subtotal[0], &total[0], 2 * l.nmode, MPI_DOUBLE, MPI_SUM,
                BoutComm::get());

  for (int kz = 0; kz < l.nmode; kz++) {
    if (!converged[kz]) {
      error_abs[kz] = sqrt(total[kz]);
      error_rel[kz] = error_abs[kz] / sqrt(total[kz + l.nmode]);
      if (error_abs[kz] < l.atol or error_rel[kz] < l.rtol) {
        converged[kz] = true;
      }
    }
  }
}

/*
 * Calculate residual on a reduced x grid. By construction, the residual is
 * zero, except at the points on the reduced grid.
 * Note: this does not synchronize the residuals between processors, as we can
 * calculate the total residual without guard cells. Coarsening requires the
 * guard cells, and an explicit synchronization is called before coarsening.
 */
void LaplaceIPT::Level::calculate_residual(const LaplaceIPT& l) {

  SCOREP0();
  if (not included) {
    return;
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.converged[kz]) {
      residual(1, kz) = rr(1, kz) - ar(l.jy, 1, kz) * xloc(0, kz)
                        - br(l.jy, 1, kz) * xloc(1, kz)
                        - cr(l.jy, 1, kz) * xloc(index_end, kz);
      if (not l.localmesh->lastX()) {
        residual(2, kz) = 0.0; // Need to ensure this, as this point is included in
                                 // residual calculations
      } else {
        residual(2, kz) = rr(2, kz) - ar(l.jy, 2, kz) * xloc(index_start, kz)
                          - br(l.jy, 2, kz) * xloc(2, kz) - cr(l.jy, 2, kz) * xloc(3, kz);
      }
    }
  }
}

/*
 * Coarsen the fine residual
 */
void LaplaceIPT::Level::coarsen(const LaplaceIPT& l,
                                const Matrix<dcomplex>& fine_residual) {

  SCOREP0();
  if (not included) {
    return;
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.converged[kz]) {
      if (not l.localmesh->lastX()) {
        residual(1, kz) = 0.25 * fine_residual(0, kz) + 0.5 * fine_residual(1, kz)
                            + 0.25 * fine_residual(3, kz);
      } else {
        // NB point(1,kz) on last proc only used on level=0
        residual(2, kz) = 0.25 * fine_residual(1, kz) + 0.5 * fine_residual(2, kz)
                            + 0.25 * fine_residual(3, kz);
      }

      // Set initial guess for coarse grid levels to zero
      for (int ix = 0; ix < 4; ix++) {
        xloc(ix, kz) = 0.0;
      }

      // Set RHS equal to residual
      rr(1, kz) = residual(1, kz);
      if (l.localmesh->lastX()) {
        rr(2, kz) = residual(2, kz);
      }
    }
  }
}

/*
 * Update the solution on the refined grid by adding the error calculated on the coarser
 * grid.
 * Note that this does not update guard cells, so we must synchronize xloc before using
 * it.
 */
void LaplaceIPT::Level::update_solution(const LaplaceIPT& l) {

  SCOREP0();
  if (not included) {
    return;
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.converged[kz]) {
      for (int ix = 1; ix < 3; ix++) {
        xloc(ix, kz) += l.fine_error(ix, kz);
      }
    }
  }
}

/*
 * Refine the reduced system.
 * There are three types of proc to cover:
 *  + procs included at this level. Calculate error and send contributions to neighbours
 *  + procs not included at this level but included at the refined level. Receive
 * contributions
 *  + procs included neither on this level or the level above. Do nothing
 *  Special case: last processor is always included, but must also receives if refining
 * from
 *  level 1 to level 0. It only sends if refining from level 0 to level 1.
 */
void LaplaceIPT::refine(const Level& l, const Level& lup, Matrix<dcomplex>& fine_error,
                        const Array<bool>& converged) {

  SCOREP0();
  Array<dcomplex> sendvec(nmode), recvecin(nmode), recvecout(nmode);
  MPI_Request rreqin, rreqout;

  // Included processors send their contribution to procs that are included on
  // the level above.
  // Special case: last proc sends if on level > 1, but NOT on level 1
  if (l.included and (not localmesh->lastX() or l.current_level > 1)) {
    for (int kz = 0; kz < nmode; kz++) {
      if (!converged[kz]) {
        fine_error(1, kz) = l.xloc(1, kz);
        sendvec[kz] = l.xloc(1, kz);
        if (localmesh->lastX()) {
          fine_error(2, kz) = l.xloc(2, kz);
        }
      }
    }

    if (!localmesh->lastX()) {
      MPI_Send(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, lup.proc_out, 0, BoutComm::get());
    }
    if (!localmesh->firstX()) {
      MPI_Send(&sendvec[0], nmode, MPI_DOUBLE_COMPLEX, lup.proc_in, 1, BoutComm::get());
    }
  }

  // Receive if proc is included on the level above, but not this level.
  // Special case: last proc receives if on level 1
  if ((l.included_up and not l.included)
      or (localmesh->lastX() and l.current_level == 1)) {
    if (!localmesh->firstX()) {
      MPI_Irecv(&recvecin[0], nmode, MPI_DOUBLE_COMPLEX, lup.proc_in, 0, BoutComm::get(),
                &rreqin);
    }
    if (!localmesh->lastX()) {
      MPI_Irecv(&recvecout[0], nmode, MPI_DOUBLE_COMPLEX, lup.proc_out, 1,
                BoutComm::get(), &rreqout);
    }

    for (int kz = 0; kz < nmode; kz++) {
      fine_error(1, kz) = 0.0;
    }

    if (!localmesh->firstX()) {
      MPI_Wait(&rreqin, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < nmode; kz++) {
        if (!converged[kz]) {
          fine_error(1, kz) += 0.5 * recvecin[kz];
        }
      }
    }
    if (!localmesh->lastX()) {
      MPI_Wait(&rreqout, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < nmode; kz++) {
        if (!converged[kz]) {
          fine_error(1, kz) += 0.5 * recvecout[kz];
        }
      }
    }
  }
  // Special case where we need to fill (1,kz) on final proc
  if (localmesh->lastX() and l.current_level == 1) {
    for (int kz = 0; kz < nmode; kz++) {
      if (!converged[kz]) {
        fine_error(1, kz) += 0.5 * l.xloc(2, kz);
      }
    }
  }
}

/*
 * Synchronize the values of a reduced field(4,nmode) between processors that
 * are neighbours on level l. This assumes each processor's value of
 * field(1,:) is correct, and puts the in-neighbour's value into field(0,:)
 * and out-neighbour's value into field(3,:).
 */
void LaplaceIPT::Level::synchronize_reduced_field(const LaplaceIPT& l,
                                                  Matrix<dcomplex>& field) {

  SCOREP0();
  if (not included) {
    return;
  }

  MPI_Comm comm = BoutComm::get();

  // Communicate in
  if (not l.localmesh->firstX()) {
    MPI_Sendrecv(&field(1, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, &field(0, 0),
                 l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, comm, MPI_STATUS_IGNORE);
  }

  // Communicate out
  if (not l.localmesh->lastX()) {
    MPI_Sendrecv(&field(1, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, &field(3, 0),
                 l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, comm, MPI_STATUS_IGNORE);
  }
}

/*
 * Returns the transpose of a matrix
 */
void LaplaceIPT::transpose(Matrix<dcomplex>& m_t, const Matrix<dcomplex>& m, const int n1,
                           const int n2) {
  SCOREP0();
  for (int i1 = 0; i1 < n1; i1++) {
    for (int i2 = 0; i2 < n2; i2++) {
      m_t(i1, i2) = m(i2, i1);
    }
  }
}
