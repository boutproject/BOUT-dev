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

#include "bout/build_config.hxx"
#include "iterative_parallel_tri.hxx"

#if not BOUT_USE_METRIC_3D

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
    : Laplacian(opt, loc, mesh_in),
      rtol((*opt)["rtol"].doc("Relative tolerance").withDefault(1.e-7)),
      atol((*opt)["atol"].doc("Absolute tolerance").withDefault(1.e-20)),
      maxits((*opt)["maxits"].doc("Maximum number of iterations").withDefault(100)),
      max_level((*opt)["max_level"].doc("Maximum number of coarse grids").withDefault(0)),
      max_cycle((*opt)["max_cycle"]
                    .doc("Maximum number of iterations per coarse grid")
                    .withDefault(1)),
      predict_exit((*opt)["predict_exit"]
                       .doc("Predict when convergence will be reached, and skip "
                            "expensive convergence checks at earlier iterations")
                       .withDefault(false)),
      A(0.0, localmesh), C(1.0, localmesh), D(1.0, localmesh), nmode(maxmode + 1),
      ncx(localmesh->LocalNx), ny(localmesh->LocalNy), avec(ny, nmode, ncx),
      bvec(ny, nmode, ncx), cvec(ny, nmode, ncx), upperGuardVector(ny, nmode, ncx),
      lowerGuardVector(ny, nmode, ncx), minvb(nmode, ncx), al(ny, nmode), bl(ny, nmode),
      au(ny, nmode), bu(ny, nmode), rl(nmode), ru(nmode), r1(ny, nmode), r2(ny, nmode),
      first_call(ny), x0saved(ny, 4, nmode), converged(nmode), fine_error(4, nmode) {

  A.setLocation(location);
  C.setLocation(location);
  D.setLocation(location);

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
  // Cannot use multigrid on 1 core
  if (n == 1 and max_level != 0) {
    throw BoutException("LaplaceIPT error: must have max_level=0 if using one processor. ");
  }

  static int ipt_solver_count = 1;
  bout::globals::dump.addRepeat(
      ipt_mean_its, "ipt_solver" + std::to_string(ipt_solver_count) + "_mean_its");
  ++ipt_solver_count;

  resetSolver();
}

/*
 * Reset the solver to its initial state
 */
void LaplaceIPT::resetSolver() {
  std::fill(std::begin(first_call), std::end(first_call), true);
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
bool LaplaceIPT::Level::is_diagonally_dominant(const LaplaceIPT& l) const {

  if (not included) {
    return true;
  }

  if (not included) {
    // Return true, as this contributes to an all_reduce over AND.  True here means that
    // skipped procs do not affect the result.
    return true;
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    // Check index 1 on all procs, except: the last proc only has index 1 if the
    // max_level == 0.
    if (not l.localmesh->lastX() or l.max_level == 0) {
      if (std::fabs(ar(l.jy, 1, kz)) + std::fabs(cr(l.jy, 1, kz))
          > std::fabs(br(l.jy, 1, kz))) {
        output_error.write("Rank {}, jy={}, kz={}, lower row not diagonally dominant\n",
                           BoutComm::rank(), l.jy, kz);
        output_error.flush();
        return false;
      }
    }
    // Check index 2 on final proc only.
    if (l.localmesh->lastX()) {
      if (std::fabs(ar(l.jy, 2, kz)) + std::fabs(cr(l.jy, 2, kz))
          > std::fabs(br(l.jy, 2, kz))) {
        output_error.write("Rank {}, jy={}, kz={}, upper row not diagonally dominant\n",
                           BoutComm::rank(), l.jy, kz);
        output_error.flush();
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
                                                  Matrix<dcomplex>& xk1d) const {
  SCOREP0();

  Array<dcomplex> x_lower(l.nmode), x_upper(l.nmode);

  for (int kz = 0; kz < l.nmode; kz++) {

    x_upper[kz] = xloc(3, kz);

    if (l.localmesh->firstX()) {
      x_lower[kz] = xloc(0, kz);
    } else {
      x_lower[kz] =
          (xloc(1, kz) - l.rl[kz] - l.bl(l.jy, kz) * xloc(3, kz)) / l.al(l.jy, kz);
    }
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    for (int i = 0; i < l.ncx; i++) {
      xk1d(kz, i) = l.minvb(kz, i) + l.upperGuardVector(l.jy, kz, i) * x_upper[kz]
                    + l.lowerGuardVector(l.jy, kz, i) * x_lower[kz];
    }
  }
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

  const BoutReal kwaveFactor = 2.0 * PI / getUniform(coords->zlength());

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

  // Should we store coefficients? True when matrix to be inverted is
  // constant, allowing results to be cached and work skipped
  const bool store_coefficients = not isInnerBoundaryFlagSet(INVERT_AC_GRAD)
                                  and not isOuterBoundaryFlagSet(INVERT_AC_GRAD)
                                  and not isInnerBoundaryFlagSet(INVERT_SET)
                                  and not isOuterBoundaryFlagSet(INVERT_SET);

  // Initialize levels. Note that the finest grid (level 0) has a different
  // routine to coarse grids (which generally depend on the grid one step
  // finer than itself).
  //
  // If the operator to invert doesn't change from one timestep to another,
  // much of the information for each level may be stored. Data that cannot
  // be cached (e.g. the changing right-hand sides) is calculated in init_rhs
  // below.
  levels.reserve(max_level + 1);
  if (first_call[jy] || not store_coefficients) {

    levels.emplace_back(*this);

    if (max_level > 0) {
      for (std::size_t l = 1; l < (static_cast<std::size_t>(max_level) + 1); ++l) {
        levels.emplace_back(*this, levels[l - 1], l);
      }
    }
  }

  // Compute coefficients that depend on the right-hand side and which
  // therefore change every time.
  levels[0].init_rhs(*this, bcmplx);

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
  std::size_t current_level = 0;
  bool down = true;

  auto errornorm = Array<BoutReal>(nmode);
  auto errornorm_old = Array<BoutReal>(nmode);
  constexpr BoutReal initial_error = 1e6;
  errornorm = initial_error;
  errornorm_old = initial_error;
  std::fill(std::begin(converged), std::end(converged), false);

  /// SCOREP_USER_REGION_END(initwhileloop);
  /// SCOREP_USER_REGION_DEFINE(whileloop);
  /// SCOREP_USER_REGION_BEGIN(whileloop, "while loop",///SCOREP_USER_REGION_TYPE_COMMON);

  const auto all = [](const Array<bool>& a) {
    return std::all_of(a.begin(), a.end(), [](bool v) { return v; });
  };

  // Check for convergence before loop to skip work with cvode
  levels[0].calculate_residual(*this);
  levels[0].calculate_total_residual(*this, errornorm, converged);
  bool execute_loop = not all(converged);

  while (execute_loop) {

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
            errornorm_old[kz] = errornorm[kz];
          }
        }
      }

      levels[0].calculate_residual(*this);

      if (cyclecount < 3 or cyclecount > cycle_eta - 5 or not predict_exit) {
        // Calculate the total residual. This also marks modes as converged, so the
        // algorithm cannot exit in cycles where this is not called.
        levels[0].calculate_total_residual(*this, errornorm, converged);

        // Based the error reduction per V-cycle, errornorm/errornorm_old,
        // predict when the slowest converging mode converges.
        if (cyclecount < 3 and predict_exit) {
          cycle_eta = 0;
          for (int kz = 0; kz < nmode; kz++) {
            const BoutReal ratio = errornorm[kz] / errornorm_old[kz];
            const int eta =
                std::ceil(std::log(1.0 / errornorm[kz]) / std::log(ratio));
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
      levels[current_level].refine(*this, fine_error);
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
      if (current_level == static_cast<std::size_t>(max_level)) {
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
      const bool is_dd = levels[max_level].is_diagonally_dominant(*this);

      bool global_is_dd;
      MPI_Allreduce(&is_dd, &global_is_dd, 1, MPI_C_BOOL, MPI_LAND, BoutComm::get());

      if (global_is_dd) {
        throw BoutException("LaplaceIPT error: Not converged within maxits={:d} "
                            "iterations. The coarsest grained iteration matrix is "
                            "diagonally dominant and convergence is guaranteed. Please "
                            "increase maxits, rtol, or atol and retry.",
                            maxits);
      }
      throw BoutException(
          "LaplaceIPT error: Not converged within maxits={:d} iterations. The coarsest "
          "iteration matrix is not diagonally dominant so there is no guarantee this "
          "method will converge. Consider (1) increasing maxits; or (2) increasing the "
          "number of levels (as grids become more diagonally dominant with "
          "coarsening). Using more grids may require larger NXPE.",
          maxits);
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
        throw BoutException("Non-finite xloc at {:d}, {:d}, {:d}", ix, jy, kz);
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
        throw BoutException("Non-finite xloc at {:d}, {:d}, {:d}", ix, jy, kz);
    }
  }
#endif

  ++ncalls;
  ipt_mean_its =
      (ipt_mean_its * BoutReal(ncalls - 1) + BoutReal(count)) / BoutReal(ncalls);

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
          xloc(0, kz) =
              -l.cvec(l.jy, kz, l.xs - 1) * xloc(1, kz) / l.bvec(l.jy, kz, l.xs - 1);
        }
        if (l.localmesh->lastX()) {
          xloc(3, kz) =
              -l.avec(l.jy, kz, l.xe + 1) * xloc(2, kz) / l.bvec(l.jy, kz, l.xe + 1);
        }
      }
    }
  }
}

// Initialization routine for coarser grids. Initialization depends on the grid
// one step finer, lup.
LaplaceIPT::Level::Level(const LaplaceIPT& l, const Level& lup,
                         const std::size_t current_level_in)
    : myproc(lup.myproc), current_level(current_level_in), index_start(0),
      index_end(l.localmesh->lastX() ? 2 : 3), included_up(lup.included),
      proc_in_up(lup.proc_in), proc_out_up(lup.proc_out) {

  SCOREP0();

  // 2^current_level
  const auto scale = 1 << current_level;

  // Whether this proc is involved in the multigrid calculation
  included = (myproc % scale == 0) or l.localmesh->lastX();

  if (not included) {
    return;
  }

  // Colouring of processor for Gauss-Seidel
  red = ((myproc / scale) % 2 == 0);
  black = ((myproc / scale) % 2 == 1);

  // The last processor is a special case. It is always included because of
  // the final grid point, which is treated explicitly. Otherwise it should
  // not be included in either the red or black work.
  if (l.localmesh->lastX()) {
    red = true;
    black = false;
  }

  // My neighbouring procs
  proc_in = myproc - scale;
  if (l.localmesh->lastX()) {
    proc_in += 1;
  }
  const int p = myproc + scale;
  proc_out = (p < l.nproc - 1) ? p : l.nproc - 1;

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
  xloc.reallocate(4, l.nmode);
  residual.reallocate(4, l.nmode);

  // Coefficients for the reduced iterations
  ar.reallocate(l.ny, 4, l.nmode);
  br.reallocate(l.ny, 4, l.nmode);
  cr.reallocate(l.ny, 4, l.nmode);
  rr.reallocate(4, l.nmode);
  brinv.reallocate(l.ny, 4, l.nmode);

  auto sendvec = Array<dcomplex>(3 * l.nmode);
  auto recvecin = Array<dcomplex>(3 * l.nmode);
  auto recvecout = Array<dcomplex>(3 * l.nmode);

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
LaplaceIPT::Level::Level(LaplaceIPT& l)
    : myproc(l.myproc), proc_in(myproc - 1), proc_out(myproc + 1), included(true),
      red(myproc % 2 == 0), black(myproc % 2 == 1), current_level(0), index_start(1),
      index_end(l.localmesh->lastX() ? 2 : 3) {

  // Basic definitions for conventional multigrid
  SCOREP0();

  const int ny = l.localmesh->LocalNy;

  // Coefficients for the reduced iterations
  ar.reallocate(ny, 4, l.nmode);
  br.reallocate(ny, 4, l.nmode);
  cr.reallocate(ny, 4, l.nmode);
  rr.reallocate(4, l.nmode);
  brinv.reallocate(ny, 4, l.nmode);

  residual.reallocate(4, l.nmode);
  residual = 0.0;

  // Define sizes of local coefficients
  xloc.reallocate(4, l.nmode); // Reduced grid x values

  // Work arrays
  auto evec = Array<dcomplex>(l.ncx);
  auto tmp = Array<dcomplex>(l.ncx);

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
      for (int i = 0; i < l.ncx; i++) {
        evec[i] = 0.0;
      }
      evec[l.xe + 1] = 1.0;
      tridag(&l.avec(l.jy, kz, 0), &l.bvec(l.jy, kz, 0), &l.cvec(l.jy, kz, 0),
             std::begin(evec), std::begin(tmp), l.ncx);
      for (int i = 0; i < l.ncx; i++) {
        l.upperGuardVector(l.jy, kz, i) = tmp[i];
      }
    } else {
      for (int i = 0; i < l.ncx; i++) {
        l.upperGuardVector(l.jy, kz, i) = 0.0;
      }
    }

    // Lower interface
    if (not l.localmesh->firstX()) {
      for (int i = 0; i < l.ncx; i++) {
        evec[i] = 0.0;
      }
      evec[l.xs - 1] = 1.0;
      tridag(&l.avec(l.jy, kz, 0), &l.bvec(l.jy, kz, 0), &l.cvec(l.jy, kz, 0),
             std::begin(evec), std::begin(tmp), l.ncx);
      for (int i = 0; i < l.ncx; i++) {
        l.lowerGuardVector(l.jy, kz, i) = tmp[i];
      }
    } else {
      for (int i = 0; i < l.ncx; i++) {
        l.lowerGuardVector(l.jy, kz, i) = 0.0;
      }
    }

    /// SCOREP_USER_REGION_END(invert);
    /// SCOREP_USER_REGION_DEFINE(coefs);
    /// SCOREP_USER_REGION_BEGIN(coefs, "calculate
    /// coefs",///SCOREP_USER_REGION_TYPE_COMMON);

    l.bl(l.jy, kz) = l.upperGuardVector(l.jy, kz, l.xs);
    l.al(l.jy, kz) = l.lowerGuardVector(l.jy, kz, l.xs);

    l.bu(l.jy, kz) = l.upperGuardVector(l.jy, kz, l.xe);
    l.au(l.jy, kz) = l.lowerGuardVector(l.jy, kz, l.xe);

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
      ABtmp[1] = l.bu(l.jy, kz);
      if (std::fabs(l.al(l.jy, kz)) > 1e-14) {
        ABtmp[0] = l.au(l.jy, kz) / l.al(l.jy, kz);
        ABtmp[1] -= ABtmp[0] * l.bl(l.jy, kz);
      }
      // Send these
      MPI_Send(&ABtmp[0], 2, MPI_DOUBLE_COMPLEX, proc_out, 0, BoutComm::get());
    }

    if (not l.localmesh->firstX()) {
      MPI_Wait(&req, MPI_STATUS_IGNORE);
    }

    const dcomplex Delta = 1.0 / (1.0 - l.al(l.jy, kz) * AdBd[1]);
    ar(l.jy, 1, kz) = -Delta * l.al(l.jy, kz) * AdBd[0];
    cr(l.jy, 1, kz) = -Delta * l.bl(l.jy, kz);

    l.r1(l.jy, kz) = Delta * l.al(l.jy, kz);
    l.r2(l.jy, kz) = Delta;

    // lastX is a special case having two points on the level 0 grid
    if (l.localmesh->lastX()) {
      // Note that if the first proc is also the last proc, then both alold and
      // auold are zero, and l.au = l.auold is already correct.
      if (not l.localmesh->lastX()) {
        ar(l.jy, 2, kz) = -l.au(l.jy, kz) / l.al(l.jy, kz);
        // NB this depends on previous line
        cr(l.jy, 2, kz) = -(l.bu(l.jy, kz) + ar(l.jy, 2, kz) * l.bl(l.jy, kz));
      }

      // Use BCs to replace x(xe+1) = -avec(xe+1) x(xe) / bvec(xe+1)
      //  => only bl changes
      cr(l.jy, 1, kz) =
          l.avec(l.jy, kz, l.xe + 1) * l.bl(l.jy, kz) / l.bvec(l.jy, kz, l.xe + 1);
    }

    // Now set coefficients for reduced iterations (shared by all levels)
    br(l.jy, 1, kz) = 1.0;
    br(l.jy, 2, kz) = 1.0;
    brinv(l.jy, 1, kz) = 1.0;
    brinv(l.jy, 2, kz) = 1.0;

    sendvec[kz] = ar(l.jy, 1, kz);
    sendvec[kz + l.nmode] = br(l.jy, 1, kz);
    sendvec[kz + 2 * l.nmode] = cr(l.jy, 1, kz);

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

// Init routine for finest level information that cannot be cached
void LaplaceIPT::Level::init_rhs(LaplaceIPT& l, const Matrix<dcomplex>& bcmplx) {

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
    tridag(&l.avec(l.jy, kz, 0), &l.bvec(l.jy, kz, 0), &l.cvec(l.jy, kz, 0),
           &bcmplx(kz, 0), &l.minvb(kz, 0), l.ncx);
    // Now minvb is a constant vector throughout the iterations

    /// SCOREP_USER_REGION_END(invertforrhs);
    /// SCOREP_USER_REGION_DEFINE(coefsforrhs);
    /// SCOREP_USER_REGION_BEGIN(coefsforrhs, "calculate coefs for
    /// rhs",///SCOREP_USER_REGION_TYPE_COMMON);

    l.rl[kz] = l.minvb(kz, l.xs);
    l.ru[kz] = l.minvb(kz, l.xe);

    // Boundary processor value to be overwritten when relevant
    Rd[kz] = 0.0;

    if (not l.localmesh->lastX()) {
      // Send coefficients up
      Rsendup[kz] = l.ru[kz];
      if (std::fabs(l.al(l.jy, kz)) > 1e-14) {
        Rsendup[kz] -= l.rl[kz] * l.au(l.jy, kz) / l.al(l.jy, kz);
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
    rr(1, kz) = l.r1(l.jy, kz) * Rd[kz] + l.r2(l.jy, kz) * l.rl[kz];

    // Special case for multiple points on last proc
    // Note that if the first proc is also the last proc, then both al and
    // au are zero, and l.ru is already correct.
    if (l.localmesh->lastX() and not l.localmesh->firstX()) {
      rr(2, kz) = l.ru[kz] - l.au(l.jy, kz) * l.rl[kz] / l.al(l.jy, kz);
    } else {
      rr(2, kz) = l.ru[kz];
    }
  }
}

/*
 * Sum and communicate total residual for the reduced system
 * NB This calculation assumes we are using the finest grid, level 0.
 */
void LaplaceIPT::Level::calculate_total_residual(const LaplaceIPT& l,
                                                 Array<BoutReal>& errornorm,
                                                 Array<bool>& converged) {

  SCOREP0();

  if (current_level != 0) {
    throw BoutException(
        "LaplaceIPT error: calculate_total_residual can only be called on level 0");
  }

  // Communication arrays
  auto subtotal = Array<BoutReal>(l.nmode); // local contribution to residual

  for (int kz = 0; kz < l.nmode; kz++) {
    if (!converged[kz]) {
      errornorm[kz] = 0.0;

      BoutReal w = pow( l.rtol*sqrt(pow(xloc(1, kz).real(), 2) + pow(xloc(1, kz).imag(), 2)) + l.atol , 2);
      subtotal[kz] = ( pow(residual(1, kz).real(), 2) + pow(residual(1, kz).imag(), 2) ) / w;
      if (l.localmesh->lastX()) {
        w = pow( l.rtol*sqrt(pow(xloc(2, kz).real(), 2) + pow(xloc(2, kz).imag(), 2)) + l.atol , 2);
        subtotal[kz] +=
            ( pow(residual(2, kz).real(), 2) + pow(residual(2, kz).imag(), 2) ) / w;
      }
    }
  }

  // Communication needed to ensure processors break on same iteration
  MPI_Allreduce(subtotal.begin(), errornorm.begin(), l.nmode, MPI_DOUBLE, MPI_SUM,
                BoutComm::get());

  for (int kz = 0; kz < l.nmode; kz++) {
    if (!converged[kz]) {
      errornorm[kz] = sqrt(errornorm[kz]/BoutReal(l.ncx));
      if (errornorm[kz] < 1.0) {
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

  // The residual should be calculated at index 1
  //   + for all procs on level 0
  //   + for all procs, except the last proc on other levels
  // The residual should be calculated at index 2 on the last proc only on all levels
  // To remove branching:
  //   + we calculate something for indices that should be skipped - having a non-zero
  //   value is "wrong", but the value is never used.
  //   + we use index_start and index_end variable to index correctly depending on level
  //   and whether the last processor
  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.converged[kz]) {
      residual(1, kz) = rr(1, kz) - ar(l.jy, 1, kz) * xloc(0, kz)
                        - br(l.jy, 1, kz) * xloc(1, kz)
                        - cr(l.jy, 1, kz) * xloc(index_end, kz);
      residual(2, kz) = rr(2, kz) - ar(l.jy, 2, kz) * xloc(index_start, kz)
                        - br(l.jy, 2, kz) * xloc(2, kz) - cr(l.jy, 2, kz) * xloc(3, kz);
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
void LaplaceIPT::Level::refine(const LaplaceIPT& l, Matrix<dcomplex>& fine_error) {

  SCOREP0();
  Array<dcomplex> sendvec(l.nmode), recvecin(l.nmode), recvecout(l.nmode);
  MPI_Request rreqin, rreqout;

  // Included processors send their contribution to procs that are included on
  // the level above.
  // Special case: last proc sends if on level > 1, but NOT on level 1
  if (included and (not l.localmesh->lastX() or current_level > 1)) {
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        fine_error(1, kz) = xloc(1, kz);
        sendvec[kz] = xloc(1, kz);
        if (l.localmesh->lastX()) {
          fine_error(2, kz) = xloc(2, kz);
        }
      }
    }

    if (not l.localmesh->lastX()) {
      MPI_Send(&sendvec[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_out_up, 0, BoutComm::get());
    }
    if (not l.localmesh->firstX()) {
      MPI_Send(&sendvec[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_in_up, 1, BoutComm::get());
    }
  }

  // Receive if proc is included on the level above, but not this level.
  // Special case: last proc receives if on level 1
  if ((included_up and not included) or (l.localmesh->lastX() and current_level == 1)) {
    if (not l.localmesh->firstX()) {
      MPI_Irecv(&recvecin[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_in_up, 0, BoutComm::get(),
                &rreqin);
    }
    if (not l.localmesh->lastX()) {
      MPI_Irecv(&recvecout[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_out_up, 1,
                BoutComm::get(), &rreqout);
    }

    for (int kz = 0; kz < l.nmode; kz++) {
      fine_error(1, kz) = 0.0;
    }

    if (not l.localmesh->firstX()) {
      MPI_Wait(&rreqin, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
          fine_error(1, kz) += 0.5 * recvecin[kz];
        }
      }
    }
    if (not l.localmesh->lastX()) {
      MPI_Wait(&rreqout, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
          fine_error(1, kz) += 0.5 * recvecout[kz];
        }
      }
    }
  }
  // Special case where we need to fill (1,kz) on final proc
  if (l.localmesh->lastX() and current_level == 1) {
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        fine_error(1, kz) += 0.5 * xloc(2, kz);
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
  // Send index 1 to the proc below, unless last proc and not level zero, then send 2
  const int send_in_index = (l.localmesh->lastX() and current_level != 0) ? 2 : 1;

  // Communicate in
  if (not l.localmesh->firstX()) {
    MPI_Sendrecv(&field(send_in_index, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1,
                 &field(0, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, comm,
                 MPI_STATUS_IGNORE);
  }

  // Communicate out
  if (not l.localmesh->lastX()) {
    MPI_Sendrecv(&field(1, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, &field(3, 0),
                 l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, comm, MPI_STATUS_IGNORE);
  }
}

#endif // BOUT_USE_METRIC_3D
