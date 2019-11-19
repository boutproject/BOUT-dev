/**************************************************************************
 * Interface to SUNDIALS CVODE
 *
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

#include "cvode.hxx"

#ifdef BOUT_HAS_CVODE

#include "boutcomm.hxx"
#include "boutexception.hxx"
#include "field3d.hxx"
#include "msg_stack.hxx"
#include "options.hxx"
#include "output.hxx"
#include "unused.hxx"
#include "bout/mesh.hxx"
#include "utils.hxx"

#include <cvode/cvode.h>

#if SUNDIALS_VERSION_MAJOR >= 3
#include <cvode/cvode_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#else
#include <cvode/cvode_spgmr.h>
#endif

#include <cvode/cvode_bbdpre.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>

#include <algorithm>
#include <numeric>
#include <string>

class Field2D;

#define ZERO RCONST(0.)
#define ONE RCONST(1.0)

#ifndef CVODEINT
#if SUNDIALS_VERSION_MAJOR < 3
using CVODEINT = bout::utils::function_traits<CVLocalFn>::arg_t<0>;
#else
using CVODEINT = sunindextype;
#endif
#endif

static int cvode_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data);
static int cvode_bbd_rhs(CVODEINT Nlocal, BoutReal t, N_Vector u, N_Vector du,
                         void* user_data);

static int cvode_pre(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rvec, N_Vector zvec,
                     BoutReal gamma, BoutReal delta, int lr, void* user_data);

#if SUNDIALS_VERSION_MAJOR < 3
// Shim for earlier versions
inline static int cvode_pre_shim(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rvec,
                                 N_Vector zvec, BoutReal gamma, BoutReal delta, int lr,
                                 void* user_data, N_Vector UNUSED(tmp)) {
  return cvode_pre(t, yy, yp, rvec, zvec, gamma, delta, lr, user_data);
}
#else
// Alias for newer versions
constexpr auto& cvode_pre_shim = cvode_pre;
#endif

static int cvode_jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy,
                     void* user_data, N_Vector tmp);

#if SUNDIALS_VERSION_MAJOR < 3
// Shim for earlier versions
inline int CVSpilsSetJacTimes(void* arkode_mem, std::nullptr_t,
                              CVSpilsJacTimesVecFn jtimes) {
  return CVSpilsSetJacTimesVecFn(arkode_mem, jtimes);
}
#endif

#if SUNDIALS_VERSION_MAJOR >= 4
// Shim for newer versions
inline void* CVodeCreate(int lmm, int UNUSED(iter)) { return CVodeCreate(lmm); }
constexpr auto CV_FUNCTIONAL = 0;
constexpr auto CV_NEWTON = 0;
#elif SUNDIALS_VERSION_MAJOR == 3
namespace {
constexpr auto& SUNLinSol_SPGMR = SUNSPGMR;
}
#endif

CvodeSolver::CvodeSolver(Options* opts) : Solver(opts) {
  has_constraints = false; // This solver doesn't have constraints
  canReset = true;
}

CvodeSolver::~CvodeSolver() {
  if (initialised) {
    N_VDestroy_Parallel(uvec);
    CVodeFree(&cvode_mem);
#if SUNDIALS_VERSION_MAJOR >= 3
    SUNLinSolFree(sun_solver);
#endif
#if SUNDIALS_VERSION_MAJOR >= 4
    SUNNonlinSolFree(nonlinear_solver);
#endif
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int CvodeSolver::init(int nout, BoutReal tstep) {
  TRACE("Initialising CVODE solver");

  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;

  // Save nout and tstep for use in run
  NOUT = nout;
  TIMESTEP = tstep;

  output_progress.write("Initialising SUNDIALS' CVODE solver\n");

  // Calculate number of variables (in generic_solver)
  const int local_N = getLocalN();

  // Get total problem size
  int neq;
  if (MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("Allreduce localN -> GlobalN failed!\n");
  }

  output_info.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n", n3Dvars(),
                    n2Dvars(), neq, local_N);

  // Allocate memory
  if ((uvec = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == nullptr)
    throw BoutException("SUNDIALS memory allocation failed\n");

  // Put the variables into uvec
  save_vars(NV_DATA_P(uvec));

  diagnose = (*options)["diagnose"].doc("Print solver diagnostic information?").withDefault(false);
  const auto adams_moulton = (*options)["adams_moulton"]
          .doc("Use Adams Moulton implicit multistep. Otherwise BDF method.")
          .withDefault(false);

  if (adams_moulton) {
    // By default use functional iteration for Adams-Moulton
    output_info.write("\tUsing Adams-Moulton implicit multistep method\n");
  } else {
    // Use Newton iteration for BDF
    output_info.write("\tUsing BDF method\n");
  }

  const auto lmm = adams_moulton ? CV_ADAMS : CV_BDF;
  const auto func_iter = (*options)["func_iter"].withDefault(adams_moulton);
  const auto iter = func_iter ? CV_FUNCTIONAL : CV_NEWTON;

  if ((cvode_mem = CVodeCreate(lmm, iter)) == nullptr)
    throw BoutException("CVodeCreate failed\n");

  // For callbacks, need pointer to solver object
  if (CVodeSetUserData(cvode_mem, this) < 0)
    throw BoutException("CVodeSetUserData failed\n");

  if (CVodeInit(cvode_mem, cvode_rhs, simtime, uvec) < 0)
    throw BoutException("CVodeInit failed\n");

  const auto max_order = (*options)["cvode_max_order"].doc("Maximum order of method to use. < 0 means no limit.").withDefault(-1);
  if (max_order > 0) {
    if (CVodeSetMaxOrd(cvode_mem, max_order) < 0)
      throw BoutException("CVodeSetMaxOrder failed\n");
  }

  const auto stablimdet =
      (*options)["cvode_stability_limit_detection"].withDefault(false);
  if (stablimdet) {
    if (CVodeSetStabLimDet(cvode_mem, stablimdet) < 0)
      throw BoutException("CVodeSetStabLimDet failed\n");
  }

  const auto abstol = (*options)["ATOL"].doc("Absolute tolerance").withDefault(1.0e-12);
  const auto reltol = (*options)["RTOL"].doc("Relative tolerance").withDefault(1.0e-5);
  const auto use_vector_abstol = (*options)["use_vector_abstol"].withDefault(false);
  if (use_vector_abstol) {
    std::vector<BoutReal> f2dtols;
    f2dtols.reserve(f2d.size());
    std::transform(begin(f2d), end(f2d), std::back_inserter(f2dtols),
                   [abstol](const VarStr<Field2D>& f2) {
                     auto f2_options = Options::root()[f2.name];
                     const auto wrong_name = f2_options.isSet("abstol");
                     if (wrong_name) {
                       output_warn << "WARNING: Option 'abstol' for field " << f2.name
                                   << " is deprecated. Please use 'atol' instead\n";
                     }
                     const std::string atol_name = wrong_name ? "abstol" : "atol";
                     return f2_options[atol_name].withDefault(abstol);
                   });

    std::vector<BoutReal> f3dtols;
    f3dtols.reserve(f3d.size());
    std::transform(begin(f3d), end(f3d), std::back_inserter(f3dtols),
                   [abstol](const VarStr<Field3D>& f3) {
                     return Options::root()[f3.name]["atol"].withDefault(abstol);
                   });

    N_Vector abstolvec = N_VNew_Parallel(BoutComm::get(), local_N, neq);
    if (abstolvec == nullptr)
      throw BoutException("SUNDIALS memory allocation (abstol vector) failed\n");

    set_abstol_values(NV_DATA_P(abstolvec), f2dtols, f3dtols);

    if (CVodeSVtolerances(cvode_mem, reltol, abstolvec) < 0)
      throw BoutException("CVodeSVtolerances failed\n");

    N_VDestroy_Parallel(abstolvec);
  } else {
    if (CVodeSStolerances(cvode_mem, reltol, abstol) < 0)
      throw BoutException("CVodeSStolerances failed\n");
  }

  const auto mxsteps = (*options)["mxstep"].doc("Maximum number of internal steps between outputs.").withDefault(500);
  CVodeSetMaxNumSteps(cvode_mem, mxsteps);

  const auto max_timestep = (*options)["max_timestep"].doc("Maximum time step size").withDefault(-1.0);
  if (max_timestep > 0.0) {
    CVodeSetMaxStep(cvode_mem, max_timestep);
  }

  const auto min_timestep = (*options)["min_timestep"].doc("Minimum time step size").withDefault(-1.0);
  if (min_timestep > 0.0) {
    CVodeSetMinStep(cvode_mem, min_timestep);
  }

  const auto start_timestep = (*options)["start_timestep"].doc("Starting time step. < 0 then chosen by CVODE.").withDefault(-1.0);
  if (start_timestep > 0.0) {
    CVodeSetInitStep(cvode_mem, start_timestep);
  }

  const auto mxorder = (*options)["mxorder"].withDefault(-1);
  if (mxorder > 0) {
    CVodeSetMaxOrd(cvode_mem, mxorder);
  }

  /// Newton method can include Preconditioners and Jacobian function
  if (!func_iter) {
    output_info.write("\tUsing Newton iteration\n");
    TRACE("Setting preconditioner");
    const auto maxl = (*options)["maxl"].doc("Maximum number of linear iterations").withDefault(5);
    const auto use_precon = (*options)["use_precon"].doc("Use preconditioner?").withDefault(false);

    if (use_precon) {

      const auto rightprec = (*options)["rightprec"]
                                 .doc("Use right preconditioner? Otherwise use left.")
                                 .withDefault(false);
      const int prectype = rightprec ? PREC_RIGHT : PREC_LEFT;

#if SUNDIALS_VERSION_MAJOR >= 3
      if ((sun_solver = SUNLinSol_SPGMR(uvec, prectype, maxl)) == nullptr)
        throw BoutException("Creating SUNDIALS linear solver failed\n");
      if (CVSpilsSetLinearSolver(cvode_mem, sun_solver) != CV_SUCCESS)
        throw BoutException("CVSpilsSetLinearSolver failed\n");
#else
      if (CVSpgmr(cvode_mem, prectype, maxl) != CVSPILS_SUCCESS)
        throw BoutException("CVSpgmr failed\n");
#endif

      if (!have_user_precon()) {
        output_info.write("\tUsing BBD preconditioner\n");

        /// Get options
        // Compute band_width_default from actually added fields, to allow for multiple
        // Mesh objects
        //
        // Previous implementation was equivalent to:
        //   int MXSUB = mesh->xend - mesh->xstart + 1;
        //   int band_width_default = n3Dvars()*(MXSUB+2);
        const int band_width_default = std::accumulate(
            begin(f3d), end(f3d), 0, [](int a, const VarStr<Field3D>& fvar) {
              Mesh* localmesh = fvar.var->getMesh();
              return a + localmesh->xend - localmesh->xstart + 3;
            });

        const auto mudq = (*options)["mudq"].withDefault(band_width_default);
        const auto mldq = (*options)["mldq"].withDefault(band_width_default);
        const auto mukeep = (*options)["mukeep"].withDefault(n3Dvars() + n2Dvars());
        const auto mlkeep = (*options)["mlkeep"].withDefault(n3Dvars() + n2Dvars());

        if (CVBBDPrecInit(cvode_mem, local_N, mudq, mldq, mukeep, mlkeep, ZERO,
                          cvode_bbd_rhs, nullptr))
          throw BoutException("CVBBDPrecInit failed\n");

      } else {
        output_info.write("\tUsing user-supplied preconditioner\n");

        if (CVSpilsSetPreconditioner(cvode_mem, nullptr, cvode_pre_shim))
          throw BoutException("CVSpilsSetPreconditioner failed\n");
      }
    } else {
      output_info.write("\tNo preconditioning\n");

#if SUNDIALS_VERSION_MAJOR >= 3
      if ((sun_solver = SUNLinSol_SPGMR(uvec, PREC_NONE, maxl)) == nullptr)
        throw BoutException("Creating SUNDIALS linear solver failed\n");
      if (CVSpilsSetLinearSolver(cvode_mem, sun_solver) != CV_SUCCESS)
        throw BoutException("CVSpilsSetLinearSolver failed\n");
#else
      if (CVSpgmr(cvode_mem, PREC_NONE, maxl) != CVSPILS_SUCCESS)
        throw BoutException("CVSpgmr failed\n");
#endif
    }

    /// Set Jacobian-vector multiplication function
    const auto use_jacobian = (*options)["use_jacobian"].withDefault(false);
    if ((use_jacobian) && (jacfunc != nullptr)) {
      output_info.write("\tUsing user-supplied Jacobian function\n");

      if (CVSpilsSetJacTimes(cvode_mem, nullptr, cvode_jac) != CV_SUCCESS)
        throw BoutException("CVSpilsSetJacTimesVecFn failed\n");
    } else
      output_info.write("\tUsing difference quotient approximation for Jacobian\n");
  } else {
    output_info.write("\tUsing Functional iteration\n");
#if SUNDIALS_VERSION_MAJOR >= 4
    if ((nonlinear_solver = SUNNonlinSol_FixedPoint(uvec, 0)) == nullptr)
      throw BoutException("SUNNonlinSol_FixedPoint failed\n");

    if (CVodeSetNonlinearSolver(cvode_mem, nonlinear_solver))
      throw BoutException("CVodeSetNonlinearSolver failed\n");
#endif
  }

  return 0;
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int CvodeSolver::run() {
  TRACE("CvodeSolver::run()");

  if (!initialised)
    throw BoutException("CvodeSolver not initialised\n");

  for (int i = 0; i < NOUT; i++) {

    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP);
    iteration++;

    /// Check if the run succeeded
    if (simtime < 0.0) {
      // Step failed
      throw BoutException("SUNDIALS CVODE timestep failed\n");
    }

    if (diagnose) {
      // Print additional diagnostics
      long int nsteps, nfevals, nniters, npevals, nliters;

      CVodeGetNumSteps(cvode_mem, &nsteps);
      CVodeGetNumRhsEvals(cvode_mem, &nfevals);
      CVodeGetNumNonlinSolvIters(cvode_mem, &nniters);
      CVSpilsGetNumPrecSolves(cvode_mem, &npevals);
      CVSpilsGetNumLinIters(cvode_mem, &nliters);

      output.write(
          "\nCVODE: nsteps %ld, nfevals %ld, nniters %ld, npevals %ld, nliters %ld\n",
          nsteps, nfevals, nniters, npevals, nliters);

      output.write("    -> Newton iterations per step: %e\n",
                   static_cast<BoutReal>(nniters) / static_cast<BoutReal>(nsteps));
      output.write("    -> Linear iterations per Newton iteration: %e\n",
                   static_cast<BoutReal>(nliters) / static_cast<BoutReal>(nniters));
      output.write("    -> Preconditioner evaluations per Newton: %e\n",
                   static_cast<BoutReal>(npevals) / static_cast<BoutReal>(nniters));

      // Last step size
      BoutReal last_step;
      CVodeGetLastStep(cvode_mem, &last_step);

      // Order used in last step
      int last_order;
      CVodeGetLastOrder(cvode_mem, &last_order);

      output.write("    -> Last step size: %e, order: %d\n", last_step, last_order);

      // Local error test failures
      long int num_fails;
      CVodeGetNumErrTestFails(cvode_mem, &num_fails);

      // Number of nonlinear convergence failures
      long int nonlin_fails;
      CVodeGetNumNonlinSolvConvFails(cvode_mem, &nonlin_fails);

      output.write("    -> Local error fails: %ld, nonlinear convergence fails: %ld\n",
                   num_fails, nonlin_fails);

      // Stability limit order reductions
      long int stab_lims;
      CVodeGetNumStabLimOrderReds(cvode_mem, &stab_lims);

      output.write("    -> Stability limit order reductions: %ld\n", stab_lims);
    }

    /// Call the monitor function

    if (call_monitors(simtime, i, NOUT)) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

BoutReal CvodeSolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run({})", tout);

  MPI_Barrier(BoutComm::get());

  pre_Wtime = 0.0;
  pre_ncalls = 0;

  int flag;
  if (!monitor_timestep) {
    // Run in normal mode
    flag = CVode(cvode_mem, tout, uvec, &simtime, CV_NORMAL);
  } else {
    // Run in single step mode, to call timestep monitors
    BoutReal internal_time;
    CVodeGetCurrentTime(cvode_mem, &internal_time);
    while (internal_time < tout) {
      // Run another step
      BoutReal last_time = internal_time;
      flag = CVode(cvode_mem, tout, uvec, &internal_time, CV_ONE_STEP);

      if (flag < 0) {
        throw BoutException("ERROR CVODE solve failed at t = %e, flag = %d\n",
                            internal_time, flag);
      }

      // Call timestep monitor
      call_timestep_monitors(internal_time, internal_time - last_time);
    }
    // Get output at the desired time
    flag = CVodeGetDky(cvode_mem, tout, 0, uvec);
    simtime = tout;
  }

  // Copy variables
  load_vars(NV_DATA_P(uvec));

  // Call rhs function to get extra variables at this time
  run_rhs(simtime);

  if (flag < 0) {
    throw BoutException("ERROR CVODE solve failed at t = %e, flag = %d\n", simtime, flag);
  }

  return simtime;
}

/**************************************************************************
 * RHS function du = F(t, u)
 **************************************************************************/

void CvodeSolver::rhs(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: CvodeSolver::res({})", t);

  // Load state from udata
  load_vars(udata);

  // Get the current timestep
  // Note: CVodeGetCurrentStep updated too late in older versions
  CVodeGetLastStep(cvode_mem, &hcur);

  // Call RHS function
  run_rhs(t);

  // Save derivatives to dudata
  save_derivs(dudata);
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void CvodeSolver::pre(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal* udata,
                      BoutReal* rvec, BoutReal* zvec) {
  TRACE("Running preconditioner: CvodeSolver::pre({})", t);

  BoutReal tstart = MPI_Wtime();

  int N = NV_LOCLENGTH_P(uvec);

  if (!have_user_precon()) {
    // Identity (but should never happen)
    for (int i = 0; i < N; i++)
      zvec[i] = rvec[i];
    return;
  }

  // Load state from udata (as with res function)
  load_vars(udata);

  // Load vector to be inverted into F_vars
  load_derivs(rvec);

  run_precon(t, gamma, delta);

  // Save the solution from F_vars
  save_derivs(zvec);

  pre_Wtime += MPI_Wtime() - tstart;
  pre_ncalls++;
}

/**************************************************************************
 * Jacobian-vector multiplication function
 **************************************************************************/

void CvodeSolver::jac(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata) {
  TRACE("Running Jacobian: CvodeSolver::jac({})", t);

  if (jacfunc == nullptr)
    throw BoutException("No jacobian function supplied!\n");

  // Load state from ydate
  load_vars(ydata);

  // Load vector to be multiplied into F_vars
  load_derivs(vdata);

  // Call function
  (*jacfunc)(t);

  // Save Jv from vars
  save_derivs(Jvdata);
}

/**************************************************************************
 * CVODE RHS functions
 **************************************************************************/

static int cvode_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = NV_DATA_P(u);
  BoutReal* dudata = NV_DATA_P(du);

  auto* s = static_cast<CvodeSolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

/// RHS function for BBD preconditioner
static int cvode_bbd_rhs(CVODEINT UNUSED(Nlocal), BoutReal t, N_Vector u, N_Vector du,
                         void* user_data) {
  return cvode_rhs(t, u, du, user_data);
}

/// Preconditioner function
static int cvode_pre(BoutReal t, N_Vector yy, N_Vector UNUSED(yp), N_Vector rvec,
                     N_Vector zvec, BoutReal gamma, BoutReal delta, int UNUSED(lr),
                     void* user_data) {
  BoutReal* udata = NV_DATA_P(yy);
  BoutReal* rdata = NV_DATA_P(rvec);
  BoutReal* zdata = NV_DATA_P(zvec);

  auto* s = static_cast<CvodeSolver*>(user_data);

  // Calculate residuals
  s->pre(t, gamma, delta, udata, rdata, zdata);

  return 0;
}

/// Jacobian-vector multiplication function
static int cvode_jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector UNUSED(fy),
                     void* user_data, N_Vector UNUSED(tmp)) {
  BoutReal* ydata = NV_DATA_P(y);   ///< System state
  BoutReal* vdata = NV_DATA_P(v);   ///< Input vector
  BoutReal* Jvdata = NV_DATA_P(Jv); ///< Jacobian*vector output

  auto* s = static_cast<CvodeSolver*>(user_data);

  s->jac(t, ydata, vdata, Jvdata);

  return 0;
}

/**************************************************************************
 * vector abstol functions
 **************************************************************************/

void CvodeSolver::set_abstol_values(BoutReal* abstolvec_data,
                                    std::vector<BoutReal>& f2dtols,
                                    std::vector<BoutReal>& f3dtols) {
  int p = 0; // Counter for location in abstolvec_data array

  // All boundaries
  for (const auto& i2d : bout::globals::mesh->getRegion2D("RGN_BNDRY")) {
    loop_abstol_values_op(i2d, abstolvec_data, p, f2dtols, f3dtols, true);
  }
  // Bulk of points
  for (const auto& i2d : bout::globals::mesh->getRegion2D("RGN_NOBNDRY")) {
    loop_abstol_values_op(i2d, abstolvec_data, p, f2dtols, f3dtols, false);
  }
}

void CvodeSolver::loop_abstol_values_op(Ind2D UNUSED(i2d), BoutReal* abstolvec_data,
                                        int& p, std::vector<BoutReal>& f2dtols,
                                        std::vector<BoutReal>& f3dtols, bool bndry) {
  // Loop over 2D variables
  for (std::vector<BoutReal>::size_type i = 0; i < f2dtols.size(); i++) {
    if (bndry && !f2d[i].evolve_bndry) {
      continue;
    }
    abstolvec_data[p] = f2dtols[i];
    p++;
  }

  for (int jz = 0; jz < bout::globals::mesh->LocalNz; jz++) {
    // Loop over 3D variables
    for (std::vector<BoutReal>::size_type i = 0; i < f3dtols.size(); i++) {
      if (bndry && !f3d[i].evolve_bndry) {
        continue;
      }
      abstolvec_data[p] = f3dtols[i];
      p++;
    }
  }
}

void CvodeSolver::resetInternalFields() {
  TRACE("CvodeSolver::resetInternalFields");
  save_vars(NV_DATA_P(uvec));

  if (CVodeReInit(cvode_mem, simtime, uvec) < 0)
    throw BoutException("CVodeReInit failed\n");
}

#endif
