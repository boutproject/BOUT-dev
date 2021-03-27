/**************************************************************************
 * Perpendicular Laplacian inversion. Parallel code using FFTs in z
 * and multigrid in x.
 *
 **************************************************************************
 * Copyright 2021 Joseph Parker
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

#include "1d_multigrid.hxx"
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

Laplace1DMG::Laplace1DMG(Options* opt, CELL_LOC loc, Mesh* mesh_in)
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
      bvec(ny, nmode, ncx), cvec(ny, nmode, ncx),
      first_call(ny), x0saved(ny, ncx, nmode), converged(nmode), fine_error(ncx, nmode) {

  A.setLocation(location);
  C.setLocation(location);
  D.setLocation(location);

  // Number of x grid points must be a power of 2
  const int ngx = localmesh->GlobalNx;
  const int mxg = ngx - (localmesh->xend-localmesh->xstart+1);
	  std::cout<<mxg<<"\n";
  if (!is_pow2(ngx-mxg)) {
    throw BoutException("Laplace1DMG error: internal grid points ({:d}) must be a power of 2", ngx-mxg);
  }
  // Number of procs must be a power of 2
  const int n = localmesh->NXPE;
  if (!is_pow2(n)) {
    throw BoutException("Laplace1DMG error: NXPE must be a power of 2");
  }
  // Number of levels cannot must be such that nproc <= 2^(max_level-1)
  BoutReal lognx = log2(ngx-4);
  if ( lognx <= max_level ) {
    std::cout << "WARNING: Specified max_level "<<max_level<<" is too large. Setting max_level to largest allowable level "<<lognx-1<< ", which gives one interior point and two bounarary points on finest grid. This choice is sensible and often optimal.\n";
    max_level = lognx - 1;
  }

  static int ipt_solver_count = 1;
  bout::globals::dump.addRepeat(
      ipt_mean_its, "1dmg_solver" + std::to_string(ipt_solver_count) + "_mean_its");
  ++ipt_solver_count;

  resetSolver();
}

/*
 * Reset the solver to its initial state
 */
void Laplace1DMG::resetSolver() {
  std::fill(std::begin(first_call), std::end(first_call), true);
  x0saved = 0.0;
  resetMeanIterations();
}

// TODO Move to Array
/*
 * Returns true if all values of bool array are true, otherwise returns false.
 */
bool Laplace1DMG::all(const Array<bool> a) {
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
// FieldPerp Laplace1DMG::solve(const FieldPerp& b, const FieldPerp& x0, const FieldPerp&
// b0) {
FieldPerp Laplace1DMG::solve(const FieldPerp& b, const FieldPerp& x0) {

  SCOREP0();
  Timer timer("invert"); ///< Start timer

  /// SCOREP_USER_REGION_DEFINE(initvars);
  /// SCOREP_USER_REGION_BEGIN(initvars, "init vars",///SCOREP_USER_REGION_TYPE_COMMON);

  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  TRACE("Laplace1DMG::solve(const, const)");

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
  *        Laplace1DMG::solve()
  * bk1d = The 1d array of bk
  * xk   = The fourier transformed of x, where x the output of
  *        Laplace1DMG::solve()
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
    * bounadry values),
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

//    // Patch up internal boundaries
//    if (not localmesh->lastX()) {
//      for (int ix = localmesh->xend + 1; ix < localmesh->LocalNx; ix++) {
//        avec(jy, kz, ix) = 0;
//        bvec(jy, kz, ix) = 1;
//        cvec(jy, kz, ix) = 0;
//        bcmplx(kz, ix) = 0;
//      }
//    }
//    if (not localmesh->firstX()) {
//      for (int ix = 0; ix < localmesh->xstart; ix++) {
//        avec(jy, kz, ix) = 0;
//        bvec(jy, kz, ix) = 1;
//        cvec(jy, kz, ix) = 0;
//        bcmplx(kz, ix) = 0;
//      }
//    }
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
  //std::cout<<"jy "<<jy<<" "<<first_call[jy]<<" "<<store_coefficients<<"\n";
  levels.reserve(max_level + 1);
  if (first_call[jy] || not store_coefficients) {

    levels.emplace_back(*this);

    if (max_level > 0) {
      for (std::size_t l = 1; l < (static_cast<std::size_t>(max_level) + 1); ++l) {
        levels.emplace_back(*this, levels[l - 1], l);
      }
    }
  }

  for(int lev=0; lev<=max_level; lev++){
    if(levels[lev].included){
	  std::cout<<"\n";
	  output.write("{} ar\n",lev);
	  for(int ix = 0; ix<levels[lev].nxloc; ix++){
		  output.write("{} ",levels[lev].ar(jy,ix,0).real());
	  }
	  output.write("\n{} br\n",lev);
	  for(int ix = 0; ix<levels[lev].nxloc; ix++){
		  output.write("{} ",levels[lev].br(jy,ix,0).real());
	  }
	  output.write("\n{} cr\n",lev);
	  for(int ix = 0; ix<levels[lev].nxloc; ix++){
		  output.write("{} ",levels[lev].cr(jy,ix,0).real());
	  }
	  output.write("\n");
    }
  }

  // Compute coefficients that depend on the right-hand side and which
  // therefore change every time.
  levels[0].init_rhs(*this, bcmplx);
//  std::cout<<"\nrr\n";
//  for(int ix = 0; ix<ncx; ix++){
//	  std::cout<<levels[0].rr(ix,0)<<" ";
//  }
//  std::cout<<"\n\n";

//  std::cout<<"x0saved\n";
//  for(int ix = 0; ix<ncx; ix++){
//	  std::cout<<x0saved(jy,ix,0)<<" ";
//  }
//  std::cout<<"\n";
  /// SCOREP_USER_REGION_END(initlevels);

  /// SCOREP_USER_REGION_DEFINE(setsoln);
  /// SCOREP_USER_REGION_BEGIN(setsoln, "set level 0
  /// solution",///SCOREP_USER_REGION_TYPE_COMMON);
  // Set initial values with cached values
  for (int ix = 0; ix < ncx; ix++) {
    for (int kz = 0; kz < nmode; kz++) {
      levels[0].xloc(ix, kz) = x0saved(jy, ix, kz);
    }
  }
//  levels[0].xloc(0,0) = 0.0;
//  levels[0].xloc(1,0) = -2.1414652409701063e-06;
//  //levels[0].xloc(2,0) = 2.1414652409701063e-06;
//  levels[0].xloc(3,0) =  1.1421489915894187e-05;
//  levels[0].xloc(4,0) =  1.542640865124478e-05;
//  levels[0].xloc(5,0) = 1.1242862372475545e-05;
//  levels[0].xloc(6,0) = 1.7842101647995303e-06;
//  levels[0].xloc(7,0) = -6.960278331832829e-06;
//  levels[0].xloc(8,0) = -1.0304009426179799e-05;
//  levels[0].xloc(9,0) = -6.137204678369813e-06;
//  levels[0].xloc(10,0) = 6.137204678369813e-06;
//  levels[0].xloc(11,0) = 0.0;

  output.write("xloc initial\n");
  for(int ix = 0; ix<ncx; ix++){
    output.write("{} ",levels[0].xloc(ix,0).real());
  }
  output.write("\n");
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

    output.write("before loop {}",count);
    output.write("Current level {}\n",current_level);
    if(levels[current_level].included){
      output.write("\nxloc\n");
      for(int ix = 0; ix<levels[current_level].nxloc; ix++){
        output.write("{} {}\n",levels[current_level].xloc(ix,0).real(),ix);
      }
      output.write("\n");
//	  std::cout<<"\nresidual\n";
//	  for(int ix = 0; ix<levels[current_level].nxloc; ix++){
//		  std::cout<<levels[current_level].residual(ix,0)<<" ";
//	  }
      output.write("before smoothing loop {}",count);
      output.write("\nrr\n");
      for(int ix = 0; ix<levels[current_level].nxloc; ix++){
        output.write("{} {}\n",levels[current_level].rr(ix,0).real(),ix);
      }
      output.write("\n");
    } else {
      output.write("not included\n");
    }

    if( levels[current_level].ninternal > 1 ){
      levels[current_level].gauss_seidel_red_black_local(*this);
    } else {
      output.write("\nsmoothing with ipt routine\n");
      levels[current_level].gauss_seidel_red_black(*this);
    }

    if(levels[current_level].included){
      output.write("after smoothing loop {}",count);
      output.write("\nxloc\n");
      for(int ix = 0; ix<levels[current_level].nxloc; ix++){
        output.write("{} {}\n",levels[current_level].xloc(ix,0).real(),ix);
      }
      output.write("\n");
    }

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

        output.write("cycle {}\t iteration {}\t total weighted residual {}\t reduction factor {}\n",cyclecount, count, errornorm[0], errornorm[0]/errornorm_old[0]);

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
      if( levels[current_level].proc_level < 1 ){
        levels[current_level].refine_local(*this, fine_error);
      } else {
        output.write("\nrefining with ipt routine\n");
        levels[current_level].refine(*this, fine_error);
      }

      output.write("after refine (get fine error), loop {}",count);
      output.write("\nfine_error\n");
      for(int ix = 0; ix<levels[current_level-1].nxloc; ix++){
        output.write("{} {}\n",fine_error(ix,0).real(),ix);
      }
      output.write("\n");

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

      if(levels[current_level].included){
        output.write("residul before coarsening, loop {}",count);
        output.write("\nresidual\n");
        for(int ix = 0; ix<levels[current_level].nxloc; ix++){
          output.write("{} {}\n",levels[current_level].residual(ix,0).real(),ix);
        }
        output.write("\n");
      }

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
      throw BoutException(
          "Laplace1DMG error: Not converged within maxits={:d} iterations. Sometimes multigrid does not converge. "
          "Check the error norm reduction factor in this case. If it it is less than one, multigrid will converge "
	 "given enough iterations. If so, increase maxits and rerun. "
	 "Changing solver parameters can also increase convergence speed. "
	 "Consider (1) increasing max_levels the number of levels multigrid uses. Usually, the more levels, "
	 "the faster the convergence, though one or two levels fewer than the maximum allowed is sometimes faster; "
	 "(2) changing max_cycle, the number of smoothing cycles per multigrid level. "
	 "Usually values between 1 and 5 are good. This value is a compromise: "
	 "higher values mean that a V-cycle reduces the error by a larger factor, "
	 "but requires more iterations."
	 "\n"
	 "NB: in parallel, we are assuming the cost of the algorithm is determined by communication costs rather than work required. "
	 "Therefore we assume cost is proportional to the iteration count, rather than the V-cycle count.",
          maxits);
    }
  }
/// SCOREP_USER_REGION_END(whileloop);
/// SCOREP_USER_REGION_DEFINE(afterloop);
/// SCOREP_USER_REGION_BEGIN(afterloop, "after faff",///SCOREP_USER_REGION_TYPE_COMMON);

#if CHECK > 2
  for (int ix = 0; ix < ncx; ix++) {
    for (int kz = 0; kz < nmode; kz++) {
      if (!finite(levels[0].xloc(ix, kz).real())
          or !finite(levels[0].xloc(ix, kz).imag()))
        throw BoutException("Non-finite xloc at {:d}, {:d}, {:d}", ix, jy, kz);
    }
  }
#endif

  // Cache solution
  for (int ix = 0; ix < ncx; ix++) {
    for (int kz = 0; kz < nmode; kz++) {
      x0saved(jy, ix, kz) = levels[0].xloc(ix, kz);
    }
  }

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

void Laplace1DMG::Level::gauss_seidel_red_black_local(const Laplace1DMG& l) {

  SCOREP0();

  if (not included) {
    return;
  }

//  output.write("xloc before\n");
//  for(int ix = 0; ix<nxloc; ix++){
//    output.write("{} ",xloc(ix,0).real());
//  }
//  output.write("\n");

  Array<dcomplex> sendvec(l.nmode), recvecin(l.nmode), recvecout(l.nmode);
  MPI_Request rreqin, rreqout;

  // Even pass: xloc starts in sync, so no comm before first loop,
  // but post expect receive for after first loop: xs from above
  if (not l.localmesh->firstX()) {
      MPI_Irecv(&recvecin[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, BoutComm::get(), &rreqin);
  }
  if (not l.localmesh->lastX()) {
    MPI_Irecv(&recvecout[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, BoutComm::get(), &rreqout);
  }

  int ixend = xe + 2;
  //if( l.localmesh->lastX() ) ixend = xe + 2;

  // Sweep over even x points
  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.converged[kz]) {
      for (int ix = l.xs+1; ix < ixend; ix+=2) {
	 if(kz==0){
	   output.write("ix,xl,rr,ar,xloc-,cr,xloc+,brinv: {} {} {} {} {} {} {} {}\n",ix,xloc(ix,0).real(),rr(ix,0).real(),ar(l.jy,ix,0).real(),xloc(ix-1,0).real(),cr(l.jy,ix,0).real(),xloc(ix+1,0).real(),brinv(l.jy,ix,0).real());
	 }
         xloc(ix, kz) = (rr(ix, kz) 
		      - ar(l.jy, ix, kz) * xloc(ix-1, kz)
                      - cr(l.jy, ix, kz) * xloc(ix+1, kz))*brinv(l.jy, ix, kz);
	 if(kz==0){
	   output.write("xl after: {}\n",xloc(ix,0).real());
	 }
      }
    }
  }

  if (not l.localmesh->lastX()) {
    MPI_Send(&xloc(xe, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, BoutComm::get());
  }
  if (not l.localmesh->firstX()) {
    // Receive and put data in arrays
    MPI_Wait(&rreqin, MPI_STATUS_IGNORE);
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        xloc(l.xs-1, kz) = recvecin[kz];
      }
    }
  }

//  output.write("xloc middle\n");
//  for(int ix = 0; ix<nxloc; ix++){
//    output.write("{} ",xloc(ix,0).real());
//  }
//  output.write("\n");

  // Sweep over odd x points
  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.converged[kz]) {
      for (int ix = l.xs; ix < ixend; ix+=2) {
	 if(kz==0){
	   output.write("B ix,xl,rr,ar,xloc-,cr,xloc+,brinv: {} {} {} {} {} {} {} {}\n",ix, xloc(ix,0).real(),rr(ix,0).real(),ar(l.jy,ix,0).real(),xloc(ix-1,0).real(),cr(l.jy,ix,0).real(),xloc(ix+1,0).real(),brinv(l.jy,ix,0).real());
	 }
         xloc(ix, kz) = (rr(ix, kz) 
		      - ar(l.jy, ix, kz) * xloc(ix-1, kz)
                      - cr(l.jy, ix, kz) * xloc(ix+1, kz))*brinv(l.jy, ix, kz);
	 if(kz==0){
	   output.write("xl after: {}\n",xloc(ix,0).real());
	 }
      }
    }
  }

  if (not l.localmesh->firstX()) {
    // Send my xs down
    MPI_Send(&xloc(l.xs, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, BoutComm::get());
  }
  if (not l.localmesh->lastX()) {
    MPI_Wait(&rreqout, MPI_STATUS_IGNORE);
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        xloc(xe+1, kz) = recvecout[kz];
      }
    }
  }

//  output.write("xloc before  bcs\n");
//  for(int ix = 0; ix<nxloc; ix++){
//    output.write("{} ",xloc(ix,0).real());
//  }
//  output.write("\n");

  if(current_level==0){
    // apply boundary conditions
    if (l.localmesh->firstX()) {
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
	  int ix = l.xs-1;
          xloc(ix, kz) = (rr(ix, kz) 
                      - l.cvec(l.jy, kz, ix) * xloc(ix+1, kz)) / l.bvec(l.jy, kz, ix);
        }
      }
    }
    if (l.localmesh->lastX()) {
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
	  int ix = xe+1;
          xloc(ix, kz) = (rr(ix, kz) 
		      - l.avec(l.jy, kz, ix) * xloc(ix-1, kz)) / l.bvec(l.jy, kz, ix);
	}
      }
    }
  }

//  output.write("xloc after\n");
//  for(int ix = 0; ix<nxloc; ix++){
//    output.write("{} ",xloc(ix,0).real());
//  }
//  output.write("\n");
}

/*
 * Perform Gauss-Seidel with red-black colouring on the reduced system.
 * We don't attempt comm/comp overlap, as there is not sigificant work in the
 * x loop.
 */
void Laplace1DMG::Level::gauss_seidel_red_black(const Laplace1DMG& l) {

  SCOREP0();

  output.write("included {} red {} black {}\n", included, red, black);

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
          xloc(l.xs-1, kz) = recvecin[kz];
        }
      }
    }
    if (!l.localmesh->lastX()) {
      MPI_Wait(&rreqout, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
          xloc(xe+1, kz) = recvecout[kz];
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
	int ix = l.xs;
        if(kz==0){
          output.write("ix,xl,rr,ar,xloc-,cr,xloc+,brinv: {} {} {} {} {} {} {} {}\n",ix,xloc(ix,0).real(),rr(ix,0).real(),ar(l.jy,ix,0).real(),xloc(ix-1,0).real(),cr(l.jy,ix,0).real(),xloc(ix+1,0).real(),brinv(l.jy,ix,0).real());
        }
        xloc(ix, kz) = (rr(ix, kz) - ar(l.jy, ix, kz) * xloc(ix-1, kz)
                       - cr(l.jy, ix, kz) * xloc(ix+1, kz))
                      * brinv(l.jy, ix, kz);
        if(kz==0){
          output.write("xloc after: {}\n",xloc(ix,0).real());
        }
      }
    }
    // Send same data up and down
    MPI_Send(&xloc(l.xs, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1,
             BoutComm::get()); // black never firstX
    if (not l.localmesh->lastX()) {
      MPI_Send(&xloc(l.xs, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, BoutComm::get());
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
        xloc(l.xs-1, kz) = recvecin[kz];
      }
    }
    if (!l.localmesh->lastX()) {
      MPI_Wait(&rreqout, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
          xloc(xe+1, kz) = recvecout[kz];
        }
      }
    }
  }

  // Red processors do work and comms
  if (red and not l.localmesh->lastX()) {
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
	int ix = l.xs;
        xloc(ix, kz) =
            (rr(ix, kz) - ar(l.jy, ix, kz) * xloc(ix-1, kz) - cr(l.jy, ix, kz) * xloc(ix+1, kz))
            * brinv(l.jy, ix, kz);
      }
    }
  }
  if (l.localmesh->lastX()) {
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        // index_start removes branches. On level 0, this is 1, otherwise 0
	int ix = l.xs+1;
        if(kz==0){
          output.write("A ix,xl,rr,ar,xloc-,cr,xloc+,brinv: {} {} {} {} {} {} {} {}\n",ix,xloc(ix,0).real(),rr(ix,0).real(),ar(l.jy,ix,0).real(),xloc(index_start,0).real(),cr(l.jy,ix,0).real(),xloc(ix+1,0).real(),brinv(l.jy,ix,0).real());
        }
        xloc(ix, kz) = (rr(ix, kz) 
		     - ar(l.jy, ix, kz) * xloc(index_start, kz)
                     - cr(l.jy, ix, kz) * xloc(ix+1, kz))
                      * brinv(l.jy, ix, kz);
        if(kz==0){
          output.write("xloc after: {}\n",xloc(ix,0).real());
        }
      }
    }
  }

  if (red or l.localmesh->lastX()) { // red, or last proc when not on level zero
    // Send same data up and down
    if (not l.localmesh->firstX() and not l.localmesh->lastX()) { // excludes last proc
      MPI_Send(&xloc(l.xs, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, BoutComm::get());
    } else if (l.localmesh->lastX() and proc_level != 0) { // last proc on level > 0
      MPI_Send(&xloc(l.xs+1, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, BoutComm::get());
    }
    if (not l.localmesh->lastX()) {
      MPI_Send(&xloc(l.xs, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, BoutComm::get());
    }
  }

//  if (current_level == 0) {
//    // Update boundaries to match interior points
//    // Do this after communication
//    for (int kz = 0; kz < l.nmode; kz++) {
//      if (not l.converged[kz]) {
//        if (l.localmesh->firstX()) {
//          xloc(0, kz) =
//              -l.cvec(l.jy, kz, l.xs - 1) * xloc(1, kz) / l.bvec(l.jy, kz, l.xs - 1);
//        }
//        if (l.localmesh->lastX()) {
//          xloc(3, kz) =
//              -l.avec(l.jy, kz, l.xe + 1) * xloc(2, kz) / l.bvec(l.jy, kz, l.xe + 1);
//        }
//      }
//    }
//  }
}

// Initialization routine for coarser grids. Initialization depends on the grid
// one step finer, lup.
Laplace1DMG::Level::Level(const Laplace1DMG& l, const Level& lup,
                         const std::size_t current_level_in)
    : myproc(lup.myproc), current_level(current_level_in), index_start(l.xs),
      index_end(l.xs+1), included_up(lup.included),
      proc_in_up(lup.proc_in), proc_out_up(lup.proc_out) {

  SCOREP0();

  // 2^current_level
  const auto point_scale = 1 << current_level;

  const auto nguards = l.ncx - (l.localmesh->xend - l.localmesh->xstart + 1);
  const auto nguards_upper = l.ncx - l.localmesh->xend;

  // Number of local x points for this level
  const int nGlobalInternal = (l.localmesh->GlobalNx-nguards);
  ninternal = (l.ncx-nguards) / point_scale;
  if( ninternal < 1 ) ninternal = 1;
  nxloc = nguards + ninternal;
  xe = nxloc - nguards_upper;
  //if(l.localmesh->lastX()) nxloc += 1;
  const int scale = (1 > point_scale/(l.ncx-nguards) ? 1 : point_scale/(l.ncx-nguards));

  // Current level, but offset such that the first level that has 1 point per
  // processor is level zero. This allows us to reuse the logic from the ipt
  // solver for treating the final processor (that has two points).
  proc_level = current_level - log2(l.ncx-nguards);

  output.write("Initialize level {}, nxloc {}\n",current_level,nxloc);
  output.write("l.ncx {}, scale {}, xe {}, ninternal {}, nGlobalInternal {}, point_scale {}, proc_level {}\n",l.ncx,scale,xe,ninternal,nGlobalInternal,point_scale, proc_level);



  // Whether this proc is involved in the multigrid calculation
  included = (ninternal > 1) or (myproc % scale == 0) or l.localmesh->lastX();

  output.write("included {}\n",included);

  if (not included) {
    return;
  }

  // Colouring of processor for Gauss-Seidel
  red = ((myproc / scale) % 2 == 0);
  black = ((myproc / scale) % 2 == 1);

  // The last processor is a special case. It is always included because of
  // the final grid point, which is treated explicitly. Otherwise it should
  // not be included in either the red or black work.
  if (l.localmesh->lastX() and proc_level > 0) {
    red = true;
    black = false;
  }

  // My neighbouring procs
  proc_in = myproc - scale;
  if (l.localmesh->lastX() and proc_level > 0) {
    proc_in += 1;
  }
  const int p = myproc + scale;
  proc_out = (p < l.nproc - 1) ? p : l.nproc - 1;

  output.write("red {}, black {}, proc_in {}, proc_out {}\n",red,black,proc_in,proc_out);

  // Fix indexing on highest proc when only 1 point per proc
  index_start = l.xs;
  index_end = l.xs+1;
  if (l.localmesh->lastX() and proc_level > 0) {
    index_start -= 1;
    index_end += 1;
  }

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
  xloc.reallocate(nxloc, l.nmode);
  residual.reallocate(nxloc, l.nmode);

  // Coefficients for the reduced iterations
  ar.reallocate(l.ny, nxloc, l.nmode);
  br.reallocate(l.ny, nxloc, l.nmode);
  cr.reallocate(l.ny, nxloc, l.nmode);
  rr.reallocate(nxloc, l.nmode);
  brinv.reallocate(l.ny, nxloc, l.nmode);

  auto sendvec = Array<dcomplex>(3 * l.nmode);
  auto recvecin = Array<dcomplex>(3 * l.nmode);
  auto recvecout = Array<dcomplex>(3 * l.nmode);

  for (int kz = 0; kz < l.nmode; kz++) {
    for (int ix = 1; ix < nxloc-1; ix++) {
      // fine grid index
      //int ixf = 2*ix;
      int ixf = 2*ix-2;
      if (l.localmesh->firstX() and ix==2) {
        ar(l.jy, 1, kz) = 0.0;
        br(l.jy, 1, kz) = 0.0;
        cr(l.jy, 1, kz) = 0.0;
        ar(l.jy, 2, kz) = 0.5 * lup.ar(l.jy, 2, kz);
        br(l.jy, 2, kz) = 0.5 * lup.br(l.jy, 2, kz) + 0.25 * lup.cr(l.jy, 2, kz)
                        + 0.25 * lup.ar(l.jy, 3, kz) + 0.125 * lup.br(l.jy, 3, kz);
        cr(l.jy, 2, kz) = 0.25 * lup.cr(l.jy, 2, kz) + 0.125 * lup.br(l.jy, 3, kz)
                        + 0.25 * lup.cr(l.jy, 3, kz);
      } else if (l.localmesh->lastX() and ix==nxloc-2){
	if(proc_level>0) ixf = 2*ix-3; // account for level 0 being special grid
        ar(l.jy, ix, kz) = 0.25 * lup.ar(l.jy, ixf-1, kz) + 0.125 * lup.br(l.jy, ixf-1, kz)
                          + 0.25 * lup.ar(l.jy, ixf, kz);
        br(l.jy, ix, kz) = 0.125 * lup.br(l.jy, ixf-1, kz) + 0.25 * lup.cr(l.jy, ixf-1, kz)
                          + 0.25 * lup.ar(l.jy, ixf, kz) + 0.5 * lup.br(l.jy, ixf, kz);
        cr(l.jy, ix, kz) = 0.5 * lup.cr(l.jy, ixf, kz);
      } else {
	if(l.localmesh->lastX() and proc_level>0) ixf = 2*ix-3; // account for level 0 being special grid
        ar(l.jy, ix, kz) = 0.25 * lup.ar(l.jy, ixf-1, kz) + 0.125 * lup.br(l.jy, ixf-1, kz)
                        + 0.25 * lup.ar(l.jy, ixf, kz);
        br(l.jy, ix, kz) = 0.125 * lup.br(l.jy, ixf-1, kz) + 0.25 * lup.cr(l.jy, ixf-1, kz)
                        + 0.25 * lup.ar(l.jy, ixf, kz) + 0.5 * lup.br(l.jy, ixf, kz)
                        + 0.25 * lup.cr(l.jy, ixf, kz) + 0.25 * lup.ar(l.jy, ixf+1, kz)
                        + 0.125 * lup.br(l.jy, ixf+1, kz);
        cr(l.jy, ix, kz) = 0.25 * lup.cr(l.jy, ixf, kz) + 0.125 * lup.br(l.jy, ixf+1, kz)
                        + 0.25 * lup.cr(l.jy, ixf+1, kz);
      }
      brinv(l.jy, ix, kz) = 1.0 / br(l.jy, ix, kz);
    
///      // Need to communicate my index 1 to this level's neighbours
///      // Index 2 if last proc.
///      if (not l.localmesh->lastX()) {
///        sendvec[kz] = ar(l.jy, l.xs, kz);
///        sendvec[kz + l.nmode] = br(l.jy, l.xs, kz);
///        sendvec[kz + 2 * l.nmode] = cr(l.jy, l.xs, kz);
///      } else {
///        sendvec[kz] = ar(l.jy, xe, kz);
///        sendvec[kz + l.nmode] = br(l.jy, xe, kz);
///        sendvec[kz + 2 * l.nmode] = cr(l.jy, xe, kz);
///      }
    }
  }

///  MPI_Comm comm = BoutComm::get();
///
///  // Communicate in
///  if (not l.localmesh->firstX()) {
///    MPI_Sendrecv(&sendvec[0], 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1, &recvecin[0],
///                 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, comm, MPI_STATUS_IGNORE);
///  }
///
///  // Communicate out
///  if (not l.localmesh->lastX()) {
///    MPI_Sendrecv(&sendvec[0], 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, &recvecout[0],
///                 3 * l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, comm, MPI_STATUS_IGNORE);
///  }
///
///  for (int kz = 0; kz < l.nmode; kz++) {
///    if (not l.localmesh->firstX()) {
///      ar(l.jy, l.xs-1, kz) = recvecin[kz];
///      br(l.jy, l.xs-1, kz) = recvecin[kz + l.nmode];
///      cr(l.jy, l.xs-1, kz) = recvecin[kz + 2 * l.nmode];
///      brinv(l.jy, l.xs-1, kz) = 1.0 / br(l.jy, l.xs-1, kz);
///    }
///    if (not l.localmesh->lastX()) {
///      ar(l.jy, xe+1, kz) = recvecout[kz];
///      br(l.jy, xe+1, kz) = recvecout[kz + l.nmode];
///      cr(l.jy, xe+1, kz) = recvecout[kz + 2 * l.nmode];
///      brinv(l.jy, xe+1, kz) = 1.0 / br(l.jy, xe+1, kz);
///    }
///  }
}

// Init routine for finest level
Laplace1DMG::Level::Level(Laplace1DMG& l)
    : myproc(l.myproc), proc_in(myproc - 1), proc_out(myproc + 1), included(true),
      red(myproc % 2 == 0), black(myproc % 2 == 1), current_level(0), index_start(1),
      index_end(l.xs+1) {

  // Basic definitions for conventional multigrid
  SCOREP0();

  const int ny = l.localmesh->LocalNy;

  const auto nguards = l.ncx - (l.localmesh->xend - l.localmesh->xstart + 1);

  nxloc = l.ncx;
  xe = l.localmesh->xend;
  ninternal = l.ncx-nguards;
  const int nGlobalInternal = (l.localmesh->GlobalNx-nguards);
  proc_level = - log2(l.ncx-nguards);

  std::cout << "Initialize level " <<  current_level << ", nxloc " << nxloc << "\n";
  std::cout << "l.ncx " <<  l.ncx << "\n";
  std::cout << "l.xs " <<  l.xs << "\n";
  std::cout << "xe " <<  xe << "\n";
  
  // Coefficients for the reduced iterations
  ar.reallocate(ny, l.ncx, l.nmode);
  br.reallocate(ny, l.ncx, l.nmode);
  cr.reallocate(ny, l.ncx, l.nmode);
  rr.reallocate(l.ncx, l.nmode);
  brinv.reallocate(ny, l.ncx, l.nmode);

  residual.reallocate(l.ncx, l.nmode);
  residual = 0.0;

  // Define sizes of local coefficients
  xloc.reallocate(l.ncx, l.nmode); // Reduced grid x values


  for (int kz = 0; kz < l.nmode; kz++) {
    for (int ix = 0; ix < l.ncx; ix++) {
      residual(ix, kz) = 0.0;
    }
  }
  // end basic definitions

  for (int kz = 0; kz < l.nmode; kz++) {
    for (int ix = 0; ix < l.ncx; ix++) {
      // Need to offset, as ar always has one guard cell
      ar(l.jy, ix, kz) = l.avec(l.jy, kz, ix);
      br(l.jy, ix, kz) = l.bvec(l.jy, kz, ix);
      cr(l.jy, ix, kz) = l.cvec(l.jy, kz, ix);
      brinv(l.jy, ix, kz) = 1.0/br(l.jy, ix, kz);
    }

    std::cout << "firstX "<<l.localmesh->firstX()<<"\n";
    if( l.localmesh->firstX() ){
      // First grid point is special case: need to eliminate first row so that
      // there are 2^k+1 points in total
      int ix = l.xs; 
      dcomplex b1 = l.bvec(l.jy, kz, ix-1);
      dcomplex b2 = l.bvec(l.jy, kz, ix);
      dcomplex c1 = l.cvec(l.jy, kz, ix-1);
      dcomplex a2 = l.avec(l.jy, kz, ix);
      //std::cout << b1 << " " <<b2<<" "<<a2<<" "<<c1<<"\n";
      // cr unchanged, use xs-1 row to eliminate ar dependence in xs row
      ar(l.jy, ix, kz) = 0.0;
      br(l.jy, ix, kz) = (b2 - c1*a2)/b1;
      brinv(l.jy, ix, kz) = 1.0/br(l.jy, ix, kz);
    }

    if( l.localmesh->firstX() ){
      ar(l.jy, l.xs-2, kz) = 0.0;
      br(l.jy, l.xs-2, kz) = 0.0;
      cr(l.jy, l.xs-2, kz) = 0.0;
      ar(l.jy, l.xs-1, kz) = 0.0;
      br(l.jy, l.xs-1, kz) = 0.0;
      cr(l.jy, l.xs-1, kz) = 0.0;
    }
    if( l.localmesh->lastX() ){
      ar(l.jy, l.xe+2, kz) = 0.0;
      br(l.jy, l.xe+2, kz) = 0.0;
      cr(l.jy, l.xe+2, kz) = 0.0;
    }
  }
}

// Init routine for finest level information that cannot be cached
void Laplace1DMG::Level::init_rhs(Laplace1DMG& l, const Matrix<dcomplex>& bcmplx) {

  SCOREP0();

  //std::cout<<"init rhs\n";
  //std::cout << l.ncx << "\n";
  for (int kz = 0; kz < l.nmode; kz++) {
    //for (int ix = l.localmesh->xstart; ix < l.localmesh->xend; ix++) {
    for (int ix = 0; ix < l.ncx; ix++) {
      //std::cout << ix << " " << kz << "\n";
      //std::cout << bcmplx(kz,ix) << "\n";

      rr(ix, kz) = bcmplx(kz, ix);

    }
  }
//  output.write("rr \n");
//  for(int ix = 0; ix<l.ncx; ix++){
//    output.write("{} ",rr(ix,0).real());
//  }
//  output.write("\n");
}

/*
 * Sum and communicate total residual for the reduced system
 * NB This calculation assumes we are using the finest grid, level 0.
 */
void Laplace1DMG::Level::calculate_total_residual(const Laplace1DMG& l,
                                                 Array<BoutReal>& errornorm,
                                                 Array<bool>& converged) {

  SCOREP0();

  if (current_level != 0) {
    throw BoutException(
        "Laplace1DMG error: calculate_total_residual can only be called on level 0");
  }

  // Communication arrays
  auto subtotal = Array<BoutReal>(l.nmode); // local contribution to residual

  for (int kz = 0; kz < l.nmode; kz++) {
    if (!converged[kz]) {
      errornorm[kz] = 0.0;
      subtotal[kz] = 0.0;
      for (int ix = 2; ix < nxloc-2; ix++) {
        BoutReal w = pow( l.rtol*sqrt(pow(xloc(ix, kz).real(), 2) + pow(xloc(ix, kz).imag(), 2)) + l.atol , 2);
        subtotal[kz] += ( pow(residual(ix, kz).real(), 2) + pow(residual(ix, kz).imag(), 2) ) / w;
      }
    }
  }

  // Communication needed to ensure processors break on same iteration
  MPI_Allreduce(subtotal.begin(), errornorm.begin(), l.nmode, MPI_DOUBLE, MPI_SUM,
                BoutComm::get());

  for (int kz = 0; kz < l.nmode; kz++) {
    if (!converged[kz]) {
	    output.write("\nncx {}\n",l.ncx);
	    output.write("\nGlobal nx {}\n",l.localmesh->GlobalNx);
      errornorm[kz] = sqrt(errornorm[kz]/BoutReal(l.localmesh->GlobalNx));
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
void Laplace1DMG::Level::calculate_residual(const Laplace1DMG& l) {

  SCOREP0();
  if (not included) {
    return;
  }

//  output.write("xloc before residual\n");
//  for(int ix = 0; ix<nxloc; ix++){
//    output.write("{} ",xloc(ix,0).real());
//  }
//  output.write("\n");

  if(l.localmesh->lastX() and proc_level>0){
    // on the last proc, data is in elements xs-1 and xs+1 only
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
	int ix = l.xs-1;
        residual(ix, kz) = rr(ix, kz) - ar(l.jy, ix, kz) * xloc(ix-1, kz)
                        - br(l.jy, ix, kz) * xloc(ix, kz)
                        - cr(l.jy, ix, kz) * xloc(ix+2, kz); // skip up two here
	ix = l.xs-1;
        residual(ix, kz) = rr(ix, kz) - ar(l.jy, ix, kz) * xloc(ix-2, kz) // skip down two here
                        - br(l.jy, ix, kz) * xloc(ix, kz)
                        - cr(l.jy, ix, kz) * xloc(ix+1, kz);
      }
    }
  } else {
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        for (int ix = 1; ix < nxloc-2; ix++) {
          residual(ix, kz) = rr(ix, kz) - ar(l.jy, ix, kz) * xloc(ix-1, kz)
                        - br(l.jy, ix, kz) * xloc(ix, kz)
                        - cr(l.jy, ix, kz) * xloc(ix+1, kz);
        }
        if( l.localmesh->lastX() ){
	  int ix = xe + 1;
          residual(ix, kz) = rr(ix, kz) - ar(l.jy, ix, kz) * xloc(ix-1, kz)
                        - br(l.jy, ix, kz) * xloc(ix, kz)
                        - cr(l.jy, ix, kz) * xloc(ix+1, kz);
	}
      }
    }
  }

//  output.write("residual\n");
//  for(int ix = 0; ix<nxloc; ix++){
//    output.write("{} ",residual(ix,0).real());
//  }
//  output.write("\n");
}

/*
 * Coarsen the fine residual
 */
void Laplace1DMG::Level::coarsen(const Laplace1DMG& l,
                                const Matrix<dcomplex>& fine_residual) {

  SCOREP0();
  if (not included) {
    return;
  }

  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.converged[kz]) {
      if(l.localmesh->lastX() and proc_level>0) {
        residual(3, kz) = 0.25 * fine_residual(2, kz) + 0.5 * fine_residual(3, kz)
                          + 0.25 * fine_residual(4, kz);

      } else {
        int ixend = xe + 1;
        if(l.localmesh->lastX()) ixend = xe + 2;

        for (int ix = l.xs; ix < ixend ; ix++) {
//        int ixf = 2*ix-1;
//	//if(current_level==1 and ix==nxloc-3) ixf = 2*ix-3;
//        residual(ix, kz) = 0.25 * fine_residual(ixf-1, kz) + 0.5 * fine_residual(ixf, kz)
//                          + 0.25 * fine_residual(ixf+1, kz);
//        // Set initial guess for coarse grid levels to zero
//        xloc(ix, kz) = 0.0;
//        // Set RHS equal to residual
//        rr(ix, kz) = residual(ix, kz);

          int ixf = 2*ix-2;
          residual(ix, kz) = 0.25 * fine_residual(ixf-1, kz) + 0.5 * fine_residual(ixf, kz)
                          + 0.25 * fine_residual(ixf+1, kz);
	  output.write("ix {} ixf {}\n",ix,ixf);
	}
      }
    }
  }

  output.write("after coarsening");
  output.write("\nresidual before sync\n");
  for(int ix = 0; ix<nxloc; ix++){
    output.write("{} {}\n",residual(ix,0).real(),ix);
  }
  output.write("\n");

  synchronize_reduced_field(l,residual);

  output.write("\nresidual after sync\n");
  for(int ix = 0; ix<nxloc; ix++){
    output.write("{} {}\n",residual(ix,0).real(),ix);
  }
  output.write("\n");

  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.converged[kz]) {
      for (int ix = 0; ix < nxloc ; ix++) {
        // Set initial guess for coarse grid levels to zero
        xloc(ix, kz) = 0.0;
        // Set RHS equal to residual
        rr(ix, kz) = residual(ix, kz);
      }
    }
  }

  output.write("\nrr\n");
  for(int ix = 0; ix<nxloc; ix++){
    output.write("{} {}\n",residual(ix,0).real(),ix);
  }
  output.write("\n");
      
}

/*
 * Update the solution on the refined grid by adding the error calculated on the coarser
 * grid.
 * Note that this does not update guard cells, so we must synchronize xloc before using
 * it.
 */
void Laplace1DMG::Level::update_solution(const Laplace1DMG& l) {

  SCOREP0();
  if (not included) {
    return;
  }

//  output.write("\nl.fine_error\n");
//	  for(int ix = 0; ix<nxloc; ix++){
//		  output.write("{} ",l.fine_error(ix,0).real());
//	  }
//          output.write("\n");
//	  std::cout<<"\nxloc\n";
//	  for(int ix = 0; ix<nxloc; ix++){
//		  std::cout<<xloc(ix,0)<<" ";
//	  }
//          std::cout<<"\n";

  for (int kz = 0; kz < l.nmode; kz++) {
    if (not l.converged[kz]) {
      for (int ix = 0; ix < nxloc; ix++) {
        xloc(ix, kz) += l.fine_error(ix, kz);
      }
    }
  }
//  output.write("xloc after update\n");
//  for(int ix = 0; ix<nxloc; ix++){
//    output.write("{} ",xloc(ix,0).real());
//  }
//  output.write("\n");
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
void Laplace1DMG::Level::refine(const Laplace1DMG& l, Matrix<dcomplex>& fine_error) {

  SCOREP0();
  Array<dcomplex> sendvec(l.nmode), recvecin(l.nmode), recvecout(l.nmode);
  MPI_Request rreqin, rreqout;

  // Included processors send their contribution to procs that are included on
  // the level above.
  // Special case: last proc sends if on level > 1, but NOT on level 1
  if (included and (not l.localmesh->lastX() or proc_level > 1)) {
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        fine_error(l.xs, kz) = xloc(l.xs, kz);
        sendvec[kz] = xloc(l.xs, kz);
        if (l.localmesh->lastX()) {
          fine_error(l.xs+1, kz) = xloc(l.xs+1, kz);
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
  if ((included_up and not included) or (l.localmesh->lastX() and proc_level == 1)) {
    if (not l.localmesh->firstX()) {
      MPI_Irecv(&recvecin[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_in_up, 0, BoutComm::get(),
                &rreqin);
    }
    if (not l.localmesh->lastX()) {
      MPI_Irecv(&recvecout[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_out_up, 1,
                BoutComm::get(), &rreqout);
    }

    for (int kz = 0; kz < l.nmode; kz++) {
      fine_error(l.xs, kz) = 0.0;
    }

    if (not l.localmesh->firstX()) {
      MPI_Wait(&rreqin, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
          fine_error(l.xs, kz) += 0.5 * recvecin[kz];
        }
      }
    }
    if (not l.localmesh->lastX()) {
      MPI_Wait(&rreqout, MPI_STATUS_IGNORE);
      for (int kz = 0; kz < l.nmode; kz++) {
        if (not l.converged[kz]) {
          fine_error(l.xs, kz) += 0.5 * recvecout[kz];
        }
      }
    }
  }
  // Special case where we need to fill (1,kz) on final proc
  if (l.localmesh->lastX() and proc_level == 1) {
    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        fine_error(l.xs-1, kz) = xloc(l.xs-1, kz);
        fine_error(l.xs, kz) = 0.5 * (xloc(l.xs-1, kz)+xloc(l.xs+1, kz));
        fine_error(l.xs+1, kz) = xloc(l.xs+1, kz);
      }
    }
  }
}

void Laplace1DMG::Level::refine_local(const Laplace1DMG& l, Matrix<dcomplex>& fine_error) {

  SCOREP0();

    for (int kz = 0; kz < l.nmode; kz++) {
      if (not l.converged[kz]) {
        int ixf;
        for (int ix = l.xs; ix < nxloc-1; ix++) {
	  ixf = 2*ix-2;
          fine_error(ixf, kz) = xloc(ix, kz);
          fine_error(ixf+1, kz) = 0.5*(xloc(ix+1,kz) + xloc(ix,kz));
        }
//        for (int ix = 1; ix < nxloc-1; ix++) {
//	  ixf = 2*ix-1;
//          fine_error(ixf-1, kz) = 0.5*(xloc(ix-1,kz) + xloc(ix,kz));
//          fine_error(ixf, kz) = xloc(ix, kz);
//        }
	fine_error(0, kz) = 0.0;
//	// Last point is special case
//	// NB: on level 1, this overrides the fine_error(ixf+1) from the line above
//	int ix = nxloc-3;
//	ixf = 2*ix-2;
//	if(current_level==1) ixf = 2*ix-3;
//	fine_error(ixf, kz) = xloc(ix, kz);
      }
    }

//  Array<dcomplex> sendvec(l.nmode), recvecin(l.nmode), recvecout(l.nmode);
//  MPI_Request rreqin, rreqout;
//
//  // Included processors send their contribution to procs that are included on
//  // the level above.
//  // Special case: last proc sends if on level > 1, but NOT on level 1
//  if (included and (not l.localmesh->lastX() or current_level > 1)) {
//    for (int kz = 0; kz < l.nmode; kz++) {
//      if (not l.converged[kz]) {
//        fine_error(1, kz) = xloc(1, kz);
//        sendvec[kz] = xloc(1, kz);
//        if (l.localmesh->lastX()) {
//          fine_error(2, kz) = xloc(2, kz);
//        }
//      }
//    }
//
//    if (not l.localmesh->lastX()) {
//      MPI_Send(&sendvec[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_out_up, 0, BoutComm::get());
//    }
//    if (not l.localmesh->firstX()) {
//      MPI_Send(&sendvec[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_in_up, 1, BoutComm::get());
//    }
//  }
//
//  // Receive if proc is included on the level above, but not this level.
//  // Special case: last proc receives if on level 1
//  if ((included_up and not included) or (l.localmesh->lastX() and current_level == 1)) {
//    if (not l.localmesh->firstX()) {
//      MPI_Irecv(&recvecin[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_in_up, 0, BoutComm::get(),
//                &rreqin);
//    }
//    if (not l.localmesh->lastX()) {
//      MPI_Irecv(&recvecout[0], l.nmode, MPI_DOUBLE_COMPLEX, proc_out_up, 1,
//                BoutComm::get(), &rreqout);
//    }
//
//    for (int kz = 0; kz < l.nmode; kz++) {
//      fine_error(1, kz) = 0.0;
//    }
//
//    if (not l.localmesh->firstX()) {
//      MPI_Wait(&rreqin, MPI_STATUS_IGNORE);
//      for (int kz = 0; kz < l.nmode; kz++) {
//        if (not l.converged[kz]) {
//          fine_error(1, kz) += 0.5 * recvecin[kz];
//        }
//      }
//    }
//    if (not l.localmesh->lastX()) {
//      MPI_Wait(&rreqout, MPI_STATUS_IGNORE);
//      for (int kz = 0; kz < l.nmode; kz++) {
//        if (not l.converged[kz]) {
//          fine_error(1, kz) += 0.5 * recvecout[kz];
//        }
//      }
//    }
//  }
//  // Special case where we need to fill (1,kz) on final proc
//  if (l.localmesh->lastX() and current_level == 1) {
//    for (int kz = 0; kz < l.nmode; kz++) {
//      if (not l.converged[kz]) {
//        fine_error(1, kz) += 0.5 * xloc(2, kz);
//      }
//    }
//  }
}

/*
 * Synchronize the values of a reduced field(4,nmode) between processors that
 * are neighbours on level l. This assumes each processor's value of
 * field(1,:) is correct, and puts the in-neighbour's value into field(0,:)
 * and out-neighbour's value into field(3,:).
 */
void Laplace1DMG::Level::synchronize_reduced_field(const Laplace1DMG& l,
                                                  Matrix<dcomplex>& field) {

  SCOREP0();
  if (not included) {
    return;
  }

  output.write("\nproc in {} proc out {}\n",proc_in,proc_out);

  MPI_Comm comm = BoutComm::get();
  // Send index 1 to the proc below, unless last proc and not level zero, then send 2
  const int send_in_index = (l.localmesh->lastX() and proc_level > 0) ? l.xs+1 : l.xs;

  // Communicate in
  if (not l.localmesh->firstX()) {
    MPI_Sendrecv(&field(send_in_index, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 1,
                 &field(l.xs-1, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_in, 0, comm,
                 MPI_STATUS_IGNORE);
  }

  // Communicate out
  if (not l.localmesh->lastX()) {
    MPI_Sendrecv(&field(xe, 0), l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 0, &field(xe+1, 0),
                 l.nmode, MPI_DOUBLE_COMPLEX, proc_out, 1, comm, MPI_STATUS_IGNORE);
  }
//  output.write("residual after sync\n");
//  for(int ix = 0; ix<nxloc; ix++){
//    output.write("{} ",field(ix,0).real());
//  }
//  output.write("\n");
}

/*
 * Returns the transpose of a matrix
 */
void Laplace1DMG::transpose(Matrix<dcomplex>& m_t, const Matrix<dcomplex>& m) {
  SCOREP0();
  const auto n1 = std::get<1>(m.shape());
  const auto n2 = std::get<0>(m.shape());
  for (int i1 = 0; i1 < n1; i1++) {
    for (int i2 = 0; i2 < n2; i2++) {
      m_t(i1, i2) = m(i2, i1);
    }
  }
}
