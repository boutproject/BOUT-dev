/**************************************************************************
 * Experimental interface to SUNDIALS ARKode IMEX solver
 *
 * NOTE: ARKode is still in beta testing so use with cautious optimism
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Nick Walkden, nick.walkden@ccfe.ac.uk
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

#include "arkode.hxx"

#ifdef BOUT_HAS_ARKODE

#include "boutcomm.hxx"
#include "boutexception.hxx"
#include "field3d.hxx"
#include "msg_stack.hxx"
#include "options.hxx"
#include "output.hxx"
#include "unused.hxx"
#include "bout/mesh.hxx"
#include "utils.hxx"

#if SUNDIALS_VERSION_MAJOR >= 4
#include <arkode/arkode_arkstep.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#else
#include <arkode/arkode.h>
#if SUNDIALS_VERSION_MAJOR >= 3
#include <arkode/arkode_spils.h>
#else
#include <arkode/arkode_spgmr.h>
#endif
#endif

#include <arkode/arkode_bbdpre.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include <algorithm>
#include <numeric>

class Field2D;

#define ZERO RCONST(0.)
#define ONE RCONST(1.0)

#ifndef ARKODEINT
#if SUNDIALS_VERSION_MAJOR < 3
using ARKODEINT = bout::utils::function_traits<ARKLocalFn>::arg_t<0>;
#else
using ARKODEINT = sunindextype;
#endif
#endif

static int arkode_rhs_explicit(BoutReal t, N_Vector u, N_Vector du, void* user_data);
static int arkode_rhs_implicit(BoutReal t, N_Vector u, N_Vector du, void* user_data);
static int arkode_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data);

static int arkode_bbd_rhs(ARKODEINT Nlocal, BoutReal t, N_Vector u, N_Vector du,
                          void* user_data);
static int arkode_pre(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rvec, N_Vector zvec,
                      BoutReal gamma, BoutReal delta, int lr, void* user_data);
#if SUNDIALS_VERSION_MAJOR < 3
// Shim for earlier versions
inline static int arkode_pre_shim(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rvec,
                                  N_Vector zvec, BoutReal gamma, BoutReal delta, int lr,
                                  void* user_data, N_Vector UNUSED(tmp)) {
  return arkode_pre(t, yy, yp, rvec, zvec, gamma, delta, lr, user_data);
}
#else
// Alias for newer versions
constexpr auto& arkode_pre_shim = arkode_pre;
#endif

static int arkode_jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy,
                      void* user_data, N_Vector tmp);
#if SUNDIALS_VERSION_MAJOR < 4
// Shim for earlier versions
inline int ARKStepSetJacTimes(void* arkode_mem, std::nullptr_t,
                              ARKSpilsJacTimesVecFn jtimes) {
#if SUNDIALS_VERSION_MAJOR < 3
  return ARKSpilsSetJacTimesVecFn(arkode_mem, jtimes);
#else
  return ARKSpilsSetJacTimes(arkode_mem, nullptr, jtimes);
#endif
}
#endif

#if SUNDIALS_VERSION_MAJOR < 4
void* ARKStepCreate(ARKRhsFn fe, ARKRhsFn fi, BoutReal t0, N_Vector y0) {
  auto arkode_mem = ARKodeCreate();

  if (arkode_mem == nullptr)
    throw BoutException("ARKodeCreate failed\n");
  if (ARKodeInit(arkode_mem, fe, fi, t0, y0) != ARK_SUCCESS)
    throw BoutException("ARKodeInit failed\n");
  return arkode_mem;
}

#if SUNDIALS_VERSION_MAJOR == 3
int ARKStepSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, std::nullptr_t) {
  return ARKSpilsSetLinearSolver(arkode_mem, LS);
}

namespace {
constexpr auto& SUNLinSol_SPGMR = SUNSPGMR;
}
#endif

// Aliases for older versions
// In SUNDIALS 4, ARKode has become ARKStep, hence all the renames
constexpr auto& ARKStepEvolve = ARKode;
constexpr auto& ARKStepFree = ARKodeFree;
constexpr auto& ARKStepGetCurrentTime = ARKodeGetCurrentTime;
constexpr auto& ARKStepGetDky = ARKodeGetDky;
constexpr auto& ARKStepGetLastStep = ARKodeGetLastStep;
constexpr auto& ARKStepGetNumLinIters = ARKSpilsGetNumLinIters;
constexpr auto& ARKStepGetNumNonlinSolvIters = ARKodeGetNumNonlinSolvIters;
constexpr auto& ARKStepGetNumPrecEvals = ARKSpilsGetNumPrecEvals;
constexpr auto& ARKStepGetNumRhsEvals = ARKodeGetNumRhsEvals;
constexpr auto& ARKStepGetNumSteps = ARKodeGetNumSteps;
constexpr auto& ARKStepReInit = ARKodeReInit;
constexpr auto& ARKStepSStolerances = ARKodeSStolerances;
constexpr auto& ARKStepSVtolerances = ARKodeSVtolerances;
constexpr auto& ARKStepSetAdaptivityMethod = ARKodeSetAdaptivityMethod;
constexpr auto& ARKStepSetCFLFraction = ARKodeSetCFLFraction;
constexpr auto& ARKStepSetEpsLin = ARKSpilsSetEpsLin;
constexpr auto& ARKStepSetExplicit = ARKodeSetExplicit;
constexpr auto& ARKStepSetFixedPoint = ARKodeSetFixedPoint;
constexpr auto& ARKStepSetFixedStep = ARKodeSetFixedStep;
constexpr auto& ARKStepSetImEx = ARKodeSetImEx;
constexpr auto& ARKStepSetImplicit = ARKodeSetImplicit;
constexpr auto& ARKStepSetInitStep = ARKodeSetInitStep;
constexpr auto& ARKStepSetLinear = ARKodeSetLinear;
constexpr auto& ARKStepSetMaxNumSteps = ARKodeSetMaxNumSteps;
constexpr auto& ARKStepSetMaxStep = ARKodeSetMaxStep;
constexpr auto& ARKStepSetMinStep = ARKodeSetMinStep;
constexpr auto& ARKStepSetOptimalParams = ARKodeSetOptimalParams;
constexpr auto& ARKStepSetOrder = ARKodeSetOrder;
constexpr auto& ARKStepSetPreconditioner = ARKSpilsSetPreconditioner;
constexpr auto& ARKStepSetUserData = ARKodeSetUserData;
#endif

ArkodeSolver::ArkodeSolver(Options* opts) : Solver(opts) {
  has_constraints = false; // This solver doesn't have constraints

  // Add diagnostics to output
  add_int_diagnostic(nsteps, "arkode_nsteps", "Cumulative number of internal steps");
  add_int_diagnostic(nfe_evals, "arkode_nfe_evals",
                     "No. of calls to fe (explicit portion of the right-hand-side "
                     "function) function");
  add_int_diagnostic(nfi_evals, "arkode_nfi_evals",
                     "No. of calls to fi (implicit portion of the right-hand-side "
                     "function) function");
  add_int_diagnostic(nniters, "arkode_nniters", "No. of nonlinear solver iterations");
  add_int_diagnostic(npevals, "arkode_npevals", "No. of preconditioner evaluations");
  add_int_diagnostic(nliters, "arkode_nliters", "No. of linear iterations");
}

ArkodeSolver::~ArkodeSolver() {
  if (initialised) {
    N_VDestroy_Parallel(uvec);
    ARKStepFree(&arkode_mem);
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

int ArkodeSolver::init(int nout, BoutReal tstep) {
  TRACE("Initialising ARKODE solver");

  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;

  // Save nout and tstep for use in run
  NOUT = nout;
  TIMESTEP = tstep;

  output.write("Initialising SUNDIALS' ARKODE solver\n");

  // Calculate number of variables (in generic_solver)
  const int local_N = getLocalN();

  // Get total problem size
  int neq;
  if (MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("Allreduce localN -> GlobalN failed!\n");
  }

  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n", n3Dvars(),
               n2Dvars(), neq, local_N);

  // Allocate memory
  if ((uvec = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == nullptr)
    throw BoutException("SUNDIALS memory allocation failed\n");

  // Put the variables into uvec
  save_vars(NV_DATA_P(uvec));

  diagnose = (*options)["diagnose"].withDefault(false);
  // Maximum number of steps to take between outputs
  const auto mxsteps = (*options)["mxstep"].withDefault(500);
  // Use ImEx capability
  const auto imex = (*options)["imex"].withDefault(true);
  // Solve only explicit part
  const auto solve_explicit = (*options)["explicit"].withDefault(true);
  // Solve only implicit part
  const auto solve_implicit = (*options)["implicit"].withDefault(true);

  ASSERT1(solve_explicit or solve_implicit);

  const auto& explicit_rhs = [imex, solve_explicit]() {
    if (imex) {
      return arkode_rhs_explicit;
    } else {
      return solve_explicit ? arkode_rhs : nullptr;
    }
  }();
  const auto& implicit_rhs = [imex, solve_implicit]() {
    if (imex) {
      return arkode_rhs_implicit;
    } else {
      return solve_implicit ? arkode_rhs : nullptr;
    }
  }();

  if ((arkode_mem = ARKStepCreate(explicit_rhs, implicit_rhs, simtime, uvec)) == nullptr)
    throw BoutException("ARKStepCreate failed\n");

  if (imex and solve_explicit and solve_implicit) {
    output_info.write("\tUsing ARKode ImEx solver \n");
    if (ARKStepSetImEx(arkode_mem) != ARK_SUCCESS)
      throw BoutException("ARKStepSetImEx failed\n");
  } else if (solve_explicit) {
    output_info.write("\tUsing ARKStep Explicit solver \n");
    if (ARKStepSetExplicit(arkode_mem) != ARK_SUCCESS)
      throw BoutException("ARKStepSetExplicit failed\n");
  } else {
    output_info.write("\tUsing ARKStep Implicit solver \n");
    if (ARKStepSetImplicit(arkode_mem) != ARK_SUCCESS)
      throw BoutException("ARKStepSetImplicit failed\n");
  }

  // For callbacks, need pointer to solver object
  if (ARKStepSetUserData(arkode_mem, this) != ARK_SUCCESS)
    throw BoutException("ARKStepSetUserData failed\n");

  // Use linear implicit solver (only evaluates jacobian inversion once
  const auto set_linear = (*options)["set_linear"].withDefault(false);
  if (set_linear) {
    output.write("\tSetting ARKStep implicit solver to Linear\n");
    if (ARKStepSetLinear(arkode_mem, 1) != ARK_SUCCESS)
      throw BoutException("ARKStepSetLinear failed\n");
  }

  // Solve explicit portion in fixed timestep mode
  // NOTE: This is not recommended except for code comparison
  const auto fixed_step = (*options)["fixed_step"].withDefault(false);
  if (fixed_step) {
    // If not given, default to adaptive timestepping
    const auto fixed_timestep = (*options)["timestep"].withDefault(0.0);
    if (ARKStepSetFixedStep(arkode_mem, fixed_timestep) != ARK_SUCCESS)
      throw BoutException("ARKStepSetFixedStep failed\n");
  }

  const auto order = (*options)["order"].withDefault(4);
  if (ARKStepSetOrder(arkode_mem, order) != ARK_SUCCESS)
    throw BoutException("ARKStepSetOrder failed\n");

  const auto cfl_frac = (*options)["cfl_frac"].withDefault(-1.0);
  if (ARKStepSetCFLFraction(arkode_mem, cfl_frac) != ARK_SUCCESS)
    throw BoutException("ARKStepSetCFLFraction failed\n");

  // Set timestep adaptivity function
  const auto adap_method = (*options)["adap_method"].withDefault(0);
  // 0 -> PID adaptivity (default)
  // 1 -> PI
  // 2 -> I
  // 3 -> explicit Gustafsson
  // 4 -> implicit Gustafsson
  // 5 -> ImEx Gustafsson

  if (ARKStepSetAdaptivityMethod(arkode_mem, adap_method, 1, 1, nullptr) != ARK_SUCCESS)
    throw BoutException("ARKStepSetAdaptivityMethod failed\n");

  const auto abstol = (*options)["ATOL"].withDefault(1.0e-12);
  const auto reltol = (*options)["RTOL"].withDefault(1.0e-5);
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

    if (ARKStepSVtolerances(arkode_mem, reltol, abstolvec) != ARK_SUCCESS)
      throw BoutException("ARKStepSVtolerances failed\n");

    N_VDestroy_Parallel(abstolvec);
  } else {
    if (ARKStepSStolerances(arkode_mem, reltol, abstol) != ARK_SUCCESS)
      throw BoutException("ARKStepSStolerances failed\n");
  }

  if (ARKStepSetMaxNumSteps(arkode_mem, mxsteps) != ARK_SUCCESS)
    throw BoutException("ARKStepSetMaxNumSteps failed\n");

  const auto max_timestep = (*options)["max_timestep"].withDefault(-1.);
  if (max_timestep > 0.0) {
    if (ARKStepSetMaxStep(arkode_mem, max_timestep) != ARK_SUCCESS)
      throw BoutException("ARKStepSetMaxStep failed\n");
  }

  const auto min_timestep = (*options)["min_timestep"].withDefault(-1.);
  if (min_timestep > 0.0) {
    if (ARKStepSetMinStep(arkode_mem, min_timestep) != ARK_SUCCESS)
      throw BoutException("ARKStepSetMinStep failed\n");
  }

  const auto start_timestep = (*options)["start_timestep"].withDefault(-1);
  if (start_timestep > 0.0) {
    if (ARKStepSetInitStep(arkode_mem, start_timestep) != ARK_SUCCESS)
      throw BoutException("ARKStepSetInitStep failed");
  }

  // ARKStepSetPredictorMethod(arkode_mem,4);

  /// Newton method can include Preconditioners and Jacobian function
  const auto fixed_point = (*options)["fixed_point"].withDefault(false);

#if SUNDIALS_VERSION_MAJOR < 4
  if (fixed_point) { // Use accellerated fixed point
    output.write("\tUsing accellerated fixed point solver\n");
    if (ARKodeSetFixedPoint(arkode_mem, 3.0))
      throw BoutException("ARKodeSetFixedPoint failed\n");
  } else {
    output.write("\tUsing Newton iteration\n");
    if (ARKodeSetNewton(arkode_mem))
      throw BoutException("ARKodeSetNewton failed\n");
  }
#else
  if (fixed_point) {
    output.write("\tUsing accelerated fixed point solver\n");
    if ((nonlinear_solver = SUNNonlinSol_FixedPoint(uvec, 3)) == nullptr)
      throw BoutException("Creating SUNDIALS fixed point nonlinear solver failed\n");
  } else {
    output.write("\tUsing Newton iteration\n");
    if ((nonlinear_solver = SUNNonlinSol_Newton(uvec)) == nullptr)
      throw BoutException("Creating SUNDIALS Newton nonlinear solver failed\n");
  }
  if (ARKStepSetNonlinearSolver(arkode_mem, nonlinear_solver) != ARK_SUCCESS)
    throw BoutException("ARKStepSetNonlinearSolver failed\n");
#endif

  const auto use_precon = (*options)["use_precon"].withDefault(false);
  const auto maxl = (*options)["maxl"].withDefault(0);

  /// Set Preconditioner
  if (use_precon) {
    const auto rightprec = (*options)["rightprec"].withDefault(false);
    const int prectype = rightprec ? PREC_RIGHT : PREC_LEFT;

#if SUNDIALS_VERSION_MAJOR >= 3
    if ((sun_solver = SUNLinSol_SPGMR(uvec, prectype, maxl)) == nullptr)
      throw BoutException("Creating SUNDIALS linear solver failed\n");
    if (ARKStepSetLinearSolver(arkode_mem, sun_solver, nullptr) != ARK_SUCCESS)
      throw BoutException("ARKStepSetLinearSolver failed\n");
#else
    if (ARKSpgmr(arkode_mem, prectype, maxl) != ARKSPILS_SUCCESS)
      throw BoutException("ARKSpgmr failed\n");
#endif

    if (!hasPreconditioner()) {
      output.write("\tUsing BBD preconditioner\n");

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

      if (ARKBBDPrecInit(arkode_mem, local_N, mudq, mldq, mukeep, mlkeep, ZERO,
                         arkode_bbd_rhs, nullptr)
          != ARK_SUCCESS)
        throw BoutException("ARKBBDPrecInit failed\n");

    } else {
      output.write("\tUsing user-supplied preconditioner\n");

      if (ARKStepSetPreconditioner(arkode_mem, nullptr, arkode_pre_shim) != ARK_SUCCESS)
        throw BoutException("ARKStepSetPreconditioner failed\n");
    }
  } else {
    // Not using preconditioning

    output.write("\tNo preconditioning\n");

#if SUNDIALS_VERSION_MAJOR >= 3
    if ((sun_solver = SUNLinSol_SPGMR(uvec, PREC_NONE, maxl)) == nullptr)
      throw BoutException("Creating SUNDIALS linear solver failed\n");
    if (ARKStepSetLinearSolver(arkode_mem, sun_solver, nullptr) != ARK_SUCCESS)
      throw BoutException("ARKStepSetLinearSolver failed\n");
#else
    if (ARKSpgmr(arkode_mem, PREC_NONE, maxl) != ARKSPILS_SUCCESS)
      throw BoutException("ARKSpgmr failed\n");
#endif
  }

  /// Set Jacobian-vector multiplication function

  const auto use_jacobian = (*options)["use_jacobian"].withDefault(false);
  if (use_jacobian and hasJacobian()) {
    output.write("\tUsing user-supplied Jacobian function\n");

    if (ARKStepSetJacTimes(arkode_mem, nullptr, arkode_jac) != ARK_SUCCESS)
      throw BoutException("ARKStepSetJacTimesVecFn failed\n");
  } else
    output.write("\tUsing difference quotient approximation for Jacobian\n");

  // Use ARKode optimal parameters
  const auto optimize = (*options)["optimize"].withDefault(false);
  if (optimize) {
    output.write("\tUsing ARKode inbuilt optimization\n");
    if (ARKStepSetOptimalParams(arkode_mem) != ARK_SUCCESS)
      throw BoutException("ARKStepSetOptimalParams failed");
  }
  return 0;
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int ArkodeSolver::run() {
  TRACE("ArkodeSolver::run()");

  if (!initialised)
    throw BoutException("ArkodeSolver not initialised\n");

  for (int i = 0; i < NOUT; i++) {

    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP);
    iteration++;

    /// Check if the run succeeded
    if (simtime < 0.0) {
      // Step failed
      output.write("Timestep failed. Aborting\n");

      throw BoutException("ARKode timestep failed\n");
    }

    // Get additional diagnostics
    long int temp_long_int, temp_long_int2;
    ARKStepGetNumSteps(arkode_mem, &temp_long_int);
    nsteps = int(temp_long_int);
    ARKStepGetNumRhsEvals(arkode_mem, &temp_long_int, &temp_long_int2);
    nfe_evals = int(temp_long_int);
    nfi_evals = int(temp_long_int2);
    ARKStepGetNumNonlinSolvIters(arkode_mem, &temp_long_int);
    nniters = int(temp_long_int);
    ARKStepGetNumPrecEvals(arkode_mem, &temp_long_int);
    npevals = int(temp_long_int);
    ARKStepGetNumLinIters(arkode_mem, &temp_long_int);
    nliters = int(temp_long_int);

    if (diagnose) {
      // Print additional diagnostics

      output.write("\nARKODE: nsteps %d, nfe_evals %d, nfi_evals %d, nniters %d, "
                   "npevals %d, nliters %d\n",
                   nsteps, nfe_evals, nfi_evals, nniters, npevals, nliters);

      output.write("    -> Newton iterations per step: %e\n",
                   static_cast<BoutReal>(nniters) / static_cast<BoutReal>(nsteps));
      output.write("    -> Linear iterations per Newton iteration: %e\n",
                   static_cast<BoutReal>(nliters) / static_cast<BoutReal>(nniters));
      output.write("    -> Preconditioner evaluations per Newton: %e\n",
                   static_cast<BoutReal>(npevals) / static_cast<BoutReal>(nniters));
    }

    /// Call the monitor function

    if (call_monitors(simtime, i, NOUT)) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

BoutReal ArkodeSolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run(%e)", tout);

  MPI_Barrier(BoutComm::get());

  pre_Wtime = 0.0;
  pre_ncalls = 0;

  int flag;
  if (!monitor_timestep) {
    // Run in normal mode
    flag = ARKStepEvolve(arkode_mem, tout, uvec, &simtime, ARK_NORMAL);
  } else {
    // Run in single step mode, to call timestep monitors
    BoutReal internal_time;
    ARKStepGetCurrentTime(arkode_mem, &internal_time);
    while (internal_time < tout) {
      // Run another step
      const BoutReal last_time = internal_time;
      flag = ARKStepEvolve(arkode_mem, tout, uvec, &internal_time, ARK_ONE_STEP);

      if (flag != ARK_SUCCESS) {
        output_error.write("ERROR ARKODE solve failed at t = %e, flag = %d\n",
                           internal_time, flag);
        return -1.0;
      }

      // Call timestep monitor
      call_timestep_monitors(internal_time, internal_time - last_time);
    }
    // Get output at the desired time
    flag = ARKStepGetDky(arkode_mem, tout, 0, uvec);
    simtime = tout;
  }

  // Copy variables
  load_vars(NV_DATA_P(uvec));
  // Call rhs function to get extra variables at this time
  run_rhs(simtime);
  // run_diffusive(simtime);
  if (flag != ARK_SUCCESS) {
    output_error.write("ERROR ARKODE solve failed at t = %e, flag = %d\n", simtime, flag);
    return -1.0;
  }

  return simtime;
}

/**************************************************************************
 * Explicit RHS function du = F_E(t, u)
 **************************************************************************/

void ArkodeSolver::rhs_e(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeSolver::rhs_e(%e)", t);

  // Load state from udata
  load_vars(udata);

  // Get the current timestep
  // Note: ARKodeGetCurrentStep updated too late in older versions
  ARKStepGetLastStep(arkode_mem, &hcur);

  // Call RHS function
  run_convective(t);

  // Save derivatives to dudata
  save_derivs(dudata);
}

/**************************************************************************
 *   Implicit RHS function du = F_I(t, u)
 **************************************************************************/

void ArkodeSolver::rhs_i(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeSolver::rhs_i(%e)", t);

  load_vars(udata);
  ARKStepGetLastStep(arkode_mem, &hcur);
  // Call Implicit RHS function
  run_diffusive(t);
  save_derivs(dudata);
}

/**************************************************************************
 *   Full  RHS function du = F(t, u)
 **************************************************************************/
void ArkodeSolver::rhs(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeSolver::rhs(%e)", t);

  load_vars(udata);
  ARKStepGetLastStep(arkode_mem, &hcur);
  // Call Implicit RHS function
  run_rhs(t);
  save_derivs(dudata);
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void ArkodeSolver::pre(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal* udata,
                       BoutReal* rvec, BoutReal* zvec) {
  TRACE("Running preconditioner: ArkodeSolver::pre(%e)", t);

  const BoutReal tstart = MPI_Wtime();

  if (!hasPreconditioner()) {
    // Identity (but should never happen)
    const int N = NV_LOCLENGTH_P(uvec);
    std::copy(rvec, rvec + N, zvec);
    return;
  }

  // Load state from udata (as with res function)
  load_vars(udata);

  // Load vector to be inverted into F_vars
  load_derivs(rvec);

  runPreconditioner(t, gamma, delta);

  // Save the solution from F_vars
  save_derivs(zvec);

  pre_Wtime += MPI_Wtime() - tstart;
  pre_ncalls++;
}

/**************************************************************************
 * Jacobian-vector multiplication function
 **************************************************************************/

void ArkodeSolver::jac(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata) {
  TRACE("Running Jacobian: ArkodeSolver::jac(%e)", t);

  if (not hasJacobian())
    throw BoutException("No jacobian function supplied!\n");

  // Load state from ydate
  load_vars(ydata);

  // Load vector to be multiplied into F_vars
  load_derivs(vdata);

  // Call function
  runJacobian(t);

  // Save Jv from vars
  save_derivs(Jvdata);
}

/**************************************************************************
 * ARKODE explicit RHS functions
 **************************************************************************/

static int arkode_rhs_explicit(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = NV_DATA_P(u);
  BoutReal* dudata = NV_DATA_P(du);

  auto* s = static_cast<ArkodeSolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs_e(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

static int arkode_rhs_implicit(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = NV_DATA_P(u);
  BoutReal* dudata = NV_DATA_P(du);

  auto* s = static_cast<ArkodeSolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs_i(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

static int arkode_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = NV_DATA_P(u);
  BoutReal* dudata = NV_DATA_P(du);

  auto* s = static_cast<ArkodeSolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

/// RHS function for BBD preconditioner
static int arkode_bbd_rhs(ARKODEINT UNUSED(Nlocal), BoutReal t, N_Vector u, N_Vector du,
                          void* user_data) {
  return arkode_rhs_implicit(t, u, du, user_data);
}

/// Preconditioner function
static int arkode_pre(BoutReal t, N_Vector yy, N_Vector UNUSED(yp), N_Vector rvec,
                      N_Vector zvec, BoutReal gamma, BoutReal delta, int UNUSED(lr),
                      void* user_data) {
  BoutReal* udata = NV_DATA_P(yy);
  BoutReal* rdata = NV_DATA_P(rvec);
  BoutReal* zdata = NV_DATA_P(zvec);

  auto* s = static_cast<ArkodeSolver*>(user_data);

  // Calculate residuals
  s->pre(t, gamma, delta, udata, rdata, zdata);

  return 0;
}

/// Jacobian-vector multiplication function
static int arkode_jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
                      N_Vector UNUSED(fy), void* user_data, N_Vector UNUSED(tmp)) {
  BoutReal* ydata = NV_DATA_P(y);   ///< System state
  BoutReal* vdata = NV_DATA_P(v);   ///< Input vector
  BoutReal* Jvdata = NV_DATA_P(Jv); ///< Jacobian*vector output

  auto* s = static_cast<ArkodeSolver*>(user_data);

  s->jac(t, ydata, vdata, Jvdata);

  return 0;
}

/**************************************************************************
 * vector abstol functions
 **************************************************************************/

void ArkodeSolver::set_abstol_values(BoutReal* abstolvec_data,
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

void ArkodeSolver::loop_abstol_values_op(Ind2D UNUSED(i2d), BoutReal* abstolvec_data,
                                         int& p, std::vector<BoutReal>& f2dtols,
                                         std::vector<BoutReal>& f3dtols, bool bndry) {
  // Loop over 2D variables
  for (std::vector<BoutReal>::size_type i = 0; i < f2dtols.size(); i++) {
    if (bndry && !f2d[i].evolve_bndry)
      continue;
    abstolvec_data[p] = f2dtols[i];
    p++;
  }

  for (int jz = 0; jz < bout::globals::mesh->LocalNz; jz++) {
    // Loop over 3D variables
    for (std::vector<BoutReal>::size_type i = 0; i < f3dtols.size(); i++) {
      if (bndry && !f3d[i].evolve_bndry)
        continue;
      abstolvec_data[p] = f3dtols[i];
      p++;
    }
  }
}

#endif
