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

#include "bout/build_config.hxx"

#include "cvode.hxx"

#if BOUT_HAS_CVODE

#include "bout/bout_enum_class.hxx"
#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/msg_stack.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"
#include "bout/unused.hxx"
#include "bout/utils.hxx"

#include "fmt/core.h"

#include <cvode/cvode.h>

#if SUNDIALS_VERSION_MAJOR >= 3
#include <cvode/cvode_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#else
#include <cvode/cvode_spgmr.h>
#endif

#include <cvode/cvode_bbdpre.h>
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

BOUT_ENUM_CLASS(positivity_constraint, none, positive, non_negative, negative,
                non_positive);

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
constexpr auto CV_FUNCTIONAL = 0;
constexpr auto CV_NEWTON = 0;
#endif

#if SUNDIALS_VERSION_MAJOR >= 3
void* CVodeCreate(int lmm, [[maybe_unused]] int iter,
                  [[maybe_unused]] SUNContext context) {
#if SUNDIALS_VERSION_MAJOR == 3
  return CVodeCreate(lmm, iter);
#elif SUNDIALS_VERSION_MAJOR == 4 || SUNDIALS_VERSION_MAJOR == 5
  return CVodeCreate(lmm);
#else
  return CVodeCreate(lmm, context);
#endif
}
#endif

CvodeSolver::CvodeSolver(Options* opts)
    : Solver(opts), diagnose((*options)["diagnose"]
                                 .doc("Print solver diagnostic information?")
                                 .withDefault(false)),
      adams_moulton(
          (*options)["adams_moulton"]
              .doc("Use Adams Moulton implicit multistep. Otherwise BDF method.")
              .withDefault(false)),
      func_iter((*options)["func_iter"]
                    .doc("Use functional iteration instead of Newton")
                    .withDefault(adams_moulton)),
      max_order((*options)["cvode_max_order"]
                    .doc("Maximum order of method to use. < 0 means no limit.")
                    .withDefault(-1)),
      stablimdet((*options)["cvode_stability_limit_detection"].withDefault(false)),
      abstol((*options)["atol"].doc("Absolute tolerance").withDefault(1.0e-12)),
      reltol((*options)["rtol"].doc("Relative tolerance").withDefault(1.0e-5)),
      use_vector_abstol((*options)["use_vector_abstol"]
                            .doc("Use separate absolute tolerance for each field")
                            .withDefault(false)),
      mxsteps((*options)["mxstep"]
                  .doc("Maximum number of internal steps between outputs.")
                  .withDefault(500)),
      max_timestep(
          (*options)["max_timestep"].doc("Maximum time step size").withDefault(-1.0)),
      min_timestep(
          (*options)["min_timestep"].doc("Minimum time step size").withDefault(-1.0)),
      start_timestep((*options)["start_timestep"]
                         .doc("Starting time step. < 0 then chosen by CVODE.")
                         .withDefault(-1.0)),
      mxorder((*options)["mxorder"].doc("Maximum order").withDefault(-1)),
      max_nonlinear_iterations(
          (*options)["max_nonlinear_iterations"]
              .doc("Maximum number of nonlinear iterations allowed by CVODE before "
                   "reducing "
                   "timestep. CVODE default (used if this option is negative) is 3.")
              .withDefault(-1)),
      apply_positivity_constraints(
          (*options)["apply_positivity_constraints"]
              .doc("Use CVODE function CVodeSetConstraints to constrain variables - the "
                   "constraint "
                   "to be applied is set by the positivity_constraint option in the "
                   "subsection for "
                   "each variable")
              .withDefault(false)),
      maxl((*options)["maxl"].doc("Maximum number of linear iterations").withDefault(5)),
      use_precon((*options)["use_precon"].doc("Use preconditioner?").withDefault(false)),
      rightprec((*options)["rightprec"]
                    .doc("Use right preconditioner? Otherwise use left.")
                    .withDefault(false)),
      use_jacobian((*options)["use_jacobian"].withDefault(false)),
      cvode_nonlinear_convergence_coef(
          (*options)["cvode_nonlinear_convergence_coef"]
              .doc("Safety factor used in the nonlinear convergence test")
              .withDefault(0.1)),
      cvode_linear_convergence_coef(
          (*options)["cvode_linear_convergence_coef"]
              .doc("Factor by which the Krylov linear solver’s convergence test constant "
                   "is reduced from the nonlinear solver test constant.")
              .withDefault(0.05)),
      suncontext(static_cast<void*>(&BoutComm::get())) {
  has_constraints = false; // This solver doesn't have constraints
  canReset = true;

  // Add diagnostics to output
  // Needs to be in constructor not init() because init() is called after
  // Solver::outputVars()
  add_int_diagnostic(nsteps, "cvode_nsteps", "Cumulative number of internal steps");
  add_int_diagnostic(nfevals, "cvode_nfevals", "No. of calls to r.h.s.  function");
  add_int_diagnostic(nniters, "cvode_nniters", "No. of nonlinear solver iterations");
  add_int_diagnostic(npevals, "cvode_npevals", "No. of preconditioner solves");
  add_int_diagnostic(nliters, "cvode_nliters", "No. of linear iterations");
  add_BoutReal_diagnostic(last_step, "cvode_last_step",
                          "Step size used for the last step before each output");
  add_int_diagnostic(last_order, "cvode_last_order",
                     "Order used during the last step before each output");
  add_int_diagnostic(num_fails, "cvode_num_fails",
                     "No. of local error test failures that have occurred");
  add_int_diagnostic(nonlin_fails, "cvode_nonlin_fails",
                     "No. of nonlinear convergence failures");
  add_int_diagnostic(stab_lims, "cvode_stab_lims",
                     "No. of order reductions due to stability limit detection");
}

CvodeSolver::~CvodeSolver() {
  if (cvode_initialised) {
    N_VDestroy_Parallel(uvec);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(sun_solver);
    SUNNonlinSolFree(nonlinear_solver);
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int CvodeSolver::init() {
  TRACE("Initialising CVODE solver");

  Solver::init();

  output_progress.write("Initialising SUNDIALS' CVODE solver\n");

  // Calculate number of variables (in generic_solver)
  const int local_N = getLocalN();

  // Get total problem size
  int neq;
  if (bout::globals::mpi->MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    throw BoutException("Allreduce localN -> GlobalN failed!\n");
  }

  output_info.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n",
                    n3Dvars(), n2Dvars(), neq, local_N);

  // Allocate memory
  if ((uvec = N_VNew_Parallel(BoutComm::get(), local_N, neq, suncontext)) == nullptr) {
    throw BoutException("SUNDIALS memory allocation failed\n");
  }

  // Put the variables into uvec
  save_vars(NV_DATA_P(uvec));

  if (adams_moulton) {
    // By default use functional iteration for Adams-Moulton
    output_info.write("\tUsing Adams-Moulton implicit multistep method\n");
  } else {
    // Use Newton iteration for BDF
    output_info.write("\tUsing BDF method\n");
  }

  const auto lmm = adams_moulton ? CV_ADAMS : CV_BDF;
  const auto iter = func_iter ? CV_FUNCTIONAL : CV_NEWTON;

  if ((cvode_mem = CVodeCreate(lmm, iter, suncontext)) == nullptr) {
    throw BoutException("CVodeCreate failed\n");
  }

  // For callbacks, need pointer to solver object
  if (CVodeSetUserData(cvode_mem, this) < 0) {
    throw BoutException("CVodeSetUserData failed\n");
  }

  if (CVodeInit(cvode_mem, cvode_rhs, simtime, uvec) < 0) {
    throw BoutException("CVodeInit failed\n");
  }

  if (max_order > 0) {
    if (CVodeSetMaxOrd(cvode_mem, max_order) < 0) {
      throw BoutException("CVodeSetMaxOrder failed\n");
    }
  }

  if (stablimdet) {
    if (CVodeSetStabLimDet(cvode_mem, stablimdet) < 0) {
      throw BoutException("CVodeSetStabLimDet failed\n");
    }
  }

  if (use_vector_abstol) {
    std::vector<BoutReal> f2dtols;
    f2dtols.reserve(f2d.size());
    std::transform(begin(f2d), end(f2d), std::back_inserter(f2dtols),
                   [this](const VarStr<Field2D>& f2) {
                     auto& f2_options = Options::root()[f2.name];
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
                   [this](const VarStr<Field3D>& f3) {
                     return Options::root()[f3.name]["atol"].withDefault(abstol);
                   });

    N_Vector abstolvec = N_VNew_Parallel(BoutComm::get(), local_N, neq, suncontext);
    if (abstolvec == nullptr) {
      throw BoutException("SUNDIALS memory allocation (abstol vector) failed\n");
    }

    set_vector_option_values(NV_DATA_P(abstolvec), f2dtols, f3dtols);

    if (CVodeSVtolerances(cvode_mem, reltol, abstolvec) < 0) {
      throw BoutException("CVodeSVtolerances failed\n");
    }

    N_VDestroy_Parallel(abstolvec);
  } else {
    if (CVodeSStolerances(cvode_mem, reltol, abstol) < 0) {
      throw BoutException("CVodeSStolerances failed\n");
    }
  }

  CVodeSetMaxNumSteps(cvode_mem, mxsteps);

  if (max_timestep > 0.0) {
    CVodeSetMaxStep(cvode_mem, max_timestep);
  }

  if (min_timestep > 0.0) {
    CVodeSetMinStep(cvode_mem, min_timestep);
  }

  if (start_timestep > 0.0) {
    CVodeSetInitStep(cvode_mem, start_timestep);
  }

  if (mxorder > 0) {
    CVodeSetMaxOrd(cvode_mem, mxorder);
  }

  if (max_nonlinear_iterations > 0) {
    CVodeSetMaxNonlinIters(cvode_mem, max_nonlinear_iterations);
  }

#if not(SUNDIALS_VERSION_MAJOR >= 3 and SUNDIALS_VERSION_MINOR >= 2)
  if (apply_positivity_constraints) {
    throw BoutException("The apply_positivity_constraints option is only available with "
                        "SUNDIALS>=3.2.0");
  }
#else
  if (apply_positivity_constraints) {
    auto f2d_constraints = create_constraints(f2d);
    auto f3d_constraints = create_constraints(f3d);

    N_Vector constraints_vec = N_VNew_Parallel(BoutComm::get(), local_N, neq, suncontext);
    if (constraints_vec == nullptr) {
      throw BoutException("SUNDIALS memory allocation (positivity constraints vector) "
                          "failed\n");
    }

    set_vector_option_values(NV_DATA_P(constraints_vec), f2d_constraints,
                             f3d_constraints);

    if (CVodeSetConstraints(cvode_mem, constraints_vec) < 0) {
      throw BoutException("CVodeSetConstraints failed\n");
    }

    N_VDestroy_Parallel(constraints_vec);
  }
#endif

  /// Newton method can include Preconditioners and Jacobian function
  if (!func_iter) {
    output_info.write("\tUsing Newton iteration\n");
    TRACE("Setting preconditioner");
    if (use_precon) {
      const int prectype = rightprec ? SUN_PREC_RIGHT : SUN_PREC_LEFT;

#if SUNDIALS_VERSION_MAJOR >= 3
      if ((sun_solver = SUNLinSol_SPGMR(uvec, prectype, maxl, suncontext)) == nullptr) {
        throw BoutException("Creating SUNDIALS linear solver failed\n");
      }
      if (CVSpilsSetLinearSolver(cvode_mem, sun_solver) != CV_SUCCESS) {
        throw BoutException("CVSpilsSetLinearSolver failed\n");
      }
#else
      if (CVSpgmr(cvode_mem, prectype, maxl) != CVSPILS_SUCCESS) {
        throw BoutException("CVSpgmr failed\n");
      }
#endif

      if (!hasPreconditioner()) {
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
                          cvode_bbd_rhs, nullptr)) {
          throw BoutException("CVBBDPrecInit failed\n");
        }

      } else {
        output_info.write("\tUsing user-supplied preconditioner\n");

        if (CVSpilsSetPreconditioner(cvode_mem, nullptr, cvode_pre_shim)) {
          throw BoutException("CVSpilsSetPreconditioner failed\n");
        }
      }
    } else {
      output_info.write("\tNo preconditioning\n");

#if SUNDIALS_VERSION_MAJOR >= 3
      if ((sun_solver = SUNLinSol_SPGMR(uvec, SUN_PREC_NONE, maxl, suncontext))
          == nullptr) {
        throw BoutException("Creating SUNDIALS linear solver failed\n");
      }
      if (CVSpilsSetLinearSolver(cvode_mem, sun_solver) != CV_SUCCESS) {
        throw BoutException("CVSpilsSetLinearSolver failed\n");
      }
#else
      if (CVSpgmr(cvode_mem, SUN_PREC_NONE, maxl) != CVSPILS_SUCCESS) {
        throw BoutException("CVSpgmr failed\n");
      }
#endif
    }

    /// Set Jacobian-vector multiplication function
    if (use_jacobian and hasJacobian()) {
      output_info.write("\tUsing user-supplied Jacobian function\n");

      if (CVSpilsSetJacTimes(cvode_mem, nullptr, cvode_jac) != CV_SUCCESS) {
        throw BoutException("CVSpilsSetJacTimesVecFn failed\n");
      }
    } else {
      output_info.write("\tUsing difference quotient approximation for Jacobian\n");
    }
  } else {
    output_info.write("\tUsing Functional iteration\n");
#if SUNDIALS_VERSION_MAJOR >= 4
    if ((nonlinear_solver = SUNNonlinSol_FixedPoint(uvec, 0, suncontext)) == nullptr) {
      throw BoutException("SUNNonlinSol_FixedPoint failed\n");
    }

    if (CVodeSetNonlinearSolver(cvode_mem, nonlinear_solver)) {
      throw BoutException("CVodeSetNonlinearSolver failed\n");
    }
#endif
  }

  // Set internal tolerance factors
  CVodeSetNonlinConvCoef(cvode_mem, cvode_nonlinear_convergence_coef);
  CVodeSetEpsLin(cvode_mem, cvode_linear_convergence_coef);

  cvode_initialised = true;

  return 0;
}

template <class FieldType>
std::vector<BoutReal>
CvodeSolver::create_constraints(const std::vector<VarStr<FieldType>>& fields) {

  std::vector<BoutReal> constraints;
  constraints.reserve(fields.size());
  std::transform(begin(fields), end(fields), std::back_inserter(constraints),
                 [](const VarStr<FieldType>& f) {
                   auto& f_options = Options::root()[f.name];
                   const auto value =
                       f_options["positivity_constraint"]
                           .doc(fmt::format("Constraint to apply to {} if "
                                            "solver:apply_positivity_constraint=true. "
                                            "Possible values are: none (default), "
                                            "positive, non_negative, negative, or "
                                            "non_positive.",
                                            f.name))
                           .withDefault(positivity_constraint::none);
                   switch (value) {
                   case positivity_constraint::none:
                     return 0.0;
                   case positivity_constraint::positive:
                     return 2.0;
                   case positivity_constraint::non_negative:
                     return 1.0;
                   case positivity_constraint::negative:
                     return -2.0;
                   case positivity_constraint::non_positive:
                     return -1.0;
                   default:
                     throw BoutException("Incorrect value for "
                                         "positivity_constraint");
                   }
                 });
  return constraints;
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int CvodeSolver::run() {
  TRACE("CvodeSolver::run()");

  if (!cvode_initialised) {
    throw BoutException("CvodeSolver not initialised\n");
  }

  for (int i = 0; i < getNumberOutputSteps(); i++) {

    /// Run the solver for one output timestep
    simtime = run(simtime + getOutputTimestep());

    /// Check if the run succeeded
    if (simtime < 0.0) {
      // Step failed
      throw BoutException("SUNDIALS CVODE timestep failed\n");
    }

    // Get additional diagnostics
    long int temp_long_int;
    CVodeGetNumSteps(cvode_mem, &temp_long_int);
    nsteps = int(temp_long_int);
    CVodeGetNumRhsEvals(cvode_mem, &temp_long_int);
    nfevals = int(temp_long_int);
    CVodeGetNumNonlinSolvIters(cvode_mem, &temp_long_int);
    nniters = int(temp_long_int);
    CVSpilsGetNumPrecSolves(cvode_mem, &temp_long_int);
    npevals = int(temp_long_int);
    CVSpilsGetNumLinIters(cvode_mem, &temp_long_int);
    nliters = int(temp_long_int);

    // Last step size
    CVodeGetLastStep(cvode_mem, &last_step);

    // Order used in last step
    CVodeGetLastOrder(cvode_mem, &last_order);

    // Local error test failures
    CVodeGetNumErrTestFails(cvode_mem, &temp_long_int);
    num_fails = int(temp_long_int);

    // Number of nonlinear convergence failures
    CVodeGetNumNonlinSolvConvFails(cvode_mem, &temp_long_int);
    nonlin_fails = int(temp_long_int);

    // Stability limit order reductions
    CVodeGetNumStabLimOrderReds(cvode_mem, &temp_long_int);
    stab_lims = int(temp_long_int);

    if (diagnose) {
      // Print additional diagnostics
      output.write("\nCVODE: nsteps {:d}, nfevals {:d}, nniters {:d}, npevals {:d}, "
                   "nliters {:d}\n",
                   nsteps, nfevals, nniters, npevals, nliters);

      output.write("    -> Newton iterations per step: {:e}\n",
                   static_cast<BoutReal>(nniters) / static_cast<BoutReal>(nsteps));
      output.write("    -> Linear iterations per Newton iteration: {:e}\n",
                   static_cast<BoutReal>(nliters) / static_cast<BoutReal>(nniters));
      output.write("    -> Preconditioner evaluations per Newton: {:e}\n",
                   static_cast<BoutReal>(npevals) / static_cast<BoutReal>(nniters));

      output.write("    -> Last step size: {:e}, order: {:d}\n", last_step, last_order);

      output.write("    -> Local error fails: {:d}, nonlinear convergence fails: {:d}\n",
                   num_fails, nonlin_fails);

      output.write("    -> Stability limit order reductions: {:d}\n", stab_lims);
    }

    /// Call the monitor function

    if (call_monitors(simtime, i, getNumberOutputSteps())) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

BoutReal CvodeSolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run({})", tout);

  bout::globals::mpi->MPI_Barrier(BoutComm::get());

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
        throw BoutException("ERROR CVODE solve failed at t = {:e}, flag = {:d}\n",
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
    throw BoutException("ERROR CVODE solve failed at t = {:e}, flag = {:d}\n", simtime,
                        flag);
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

  BoutReal tstart = bout::globals::mpi->MPI_Wtime();

  int N = NV_LOCLENGTH_P(uvec);

  if (!hasPreconditioner()) {
    // Identity (but should never happen)
    for (int i = 0; i < N; i++) {
      zvec[i] = rvec[i];
    }
    return;
  }

  // Load state from udata (as with res function)
  load_vars(udata);

  // Load vector to be inverted into F_vars
  load_derivs(rvec);

  runPreconditioner(t, gamma, delta);

  // Save the solution from F_vars
  save_derivs(zvec);

  pre_Wtime += bout::globals::mpi->MPI_Wtime() - tstart;
  pre_ncalls++;
}

/**************************************************************************
 * Jacobian-vector multiplication function
 **************************************************************************/

void CvodeSolver::jac(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata) {
  TRACE("Running Jacobian: CvodeSolver::jac({})", t);

  if (not hasJacobian()) {
    throw BoutException("No jacobian function supplied!\n");
  }

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
 * CVODE vector option functions
 **************************************************************************/

void CvodeSolver::set_vector_option_values(BoutReal* option_data,
                                           std::vector<BoutReal>& f2dtols,
                                           std::vector<BoutReal>& f3dtols) {
  int p = 0; // Counter for location in option_data array

  // All boundaries
  for (const auto& i2d : bout::globals::mesh->getRegion2D("RGN_BNDRY")) {
    loop_vector_option_values_op(i2d, option_data, p, f2dtols, f3dtols, true);
  }
  // Bulk of points
  for (const auto& i2d : bout::globals::mesh->getRegion2D("RGN_NOBNDRY")) {
    loop_vector_option_values_op(i2d, option_data, p, f2dtols, f3dtols, false);
  }
}

void CvodeSolver::loop_vector_option_values_op(Ind2D UNUSED(i2d), BoutReal* option_data,
                                               int& p, std::vector<BoutReal>& f2dtols,
                                               std::vector<BoutReal>& f3dtols,
                                               bool bndry) {
  // Loop over 2D variables
  for (std::vector<BoutReal>::size_type i = 0; i < f2dtols.size(); i++) {
    if (bndry && !f2d[i].evolve_bndry) {
      continue;
    }
    option_data[p] = f2dtols[i];
    p++;
  }

  for (int jz = 0; jz < bout::globals::mesh->LocalNz; jz++) {
    // Loop over 3D variables
    for (std::vector<BoutReal>::size_type i = 0; i < f3dtols.size(); i++) {
      if (bndry && !f3d[i].evolve_bndry) {
        continue;
      }
      option_data[p] = f3dtols[i];
      p++;
    }
  }
}

void CvodeSolver::resetInternalFields() {
  TRACE("CvodeSolver::resetInternalFields");
  save_vars(NV_DATA_P(uvec));

  if (CVodeReInit(cvode_mem, simtime, uvec) < 0) {
    throw BoutException("CVodeReInit failed\n");
  }
}

#endif
