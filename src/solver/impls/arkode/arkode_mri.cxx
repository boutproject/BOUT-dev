/**************************************************************************
 * Experimental interface to SUNDIALS ARKode MRI solver
 *
 * NOTE: ARKode is still in beta testing so use with cautious optimism
 *
 **************************************************************************
 * Copyright 2010-2024 BOUT++ contributors
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

#include "arkode_mri.hxx"

#if BOUT_HAS_ARKODE

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

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_bbdpre.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include <algorithm>
#include <numeric>

class Field2D;

// NOLINTBEGIN(readability-identifier-length)
namespace {
int arkode_rhs_s_explicit(BoutReal t, N_Vector u, N_Vector du, void* user_data);
int arkode_rhs_s_implicit(BoutReal t, N_Vector u, N_Vector du, void* user_data);
int arkode_rhs_f_explicit(BoutReal t, N_Vector u, N_Vector du, void* user_data);
int arkode_rhs_f_implicit(BoutReal t, N_Vector u, N_Vector du, void* user_data);
int arkode_s_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data);
int arkode_f_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data);
int arkode_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data);

int arkode_s_bbd_rhs(sunindextype Nlocal, BoutReal t, N_Vector u, N_Vector du,
                   void* user_data);
int arkode_f_bbd_rhs(sunindextype Nlocal, BoutReal t, N_Vector u, N_Vector du,
                   void* user_data);                   
int arkode_s_pre(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rvec, N_Vector zvec,
               BoutReal gamma, BoutReal delta, int lr, void* user_data);
int arkode_f_pre(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rvec, N_Vector zvec,
               BoutReal gamma, BoutReal delta, int lr, void* user_data);               

int arkode_s_jac(N_Vector v, N_Vector Jv, BoutReal t, N_Vector y, N_Vector fy,
               void* user_data, N_Vector tmp);
int arkode_f_jac(N_Vector v, N_Vector Jv, BoutReal t, N_Vector y, N_Vector fy,
               void* user_data, N_Vector tmp);               
} // namespace
// NOLINTEND(readability-identifier-length)

ArkodeMRISolver::ArkodeMRISolver(Options* opts)
    : Solver(opts), diagnose((*options)["diagnose"]
                                 .doc("Print some additional diagnostics")
                                 .withDefault(false)),
      mxsteps((*options)["mxstep"]
                  .doc("Maximum number of steps to take between outputs")
                  .withDefault(500)),
      treatment((*options)["treatment"]
                    .doc("Use default capability (imex) or provide a specific treatment: "
                         "implicit or explicit")
                    .withDefault(Treatment::ImEx)),
      set_linear(
          (*options)["set_linear"]
              .doc("Use linear implicit solver (only evaluates jacobian inversion once)")
              .withDefault(false)),
      fixed_step((*options)["fixed_step"]
                     .doc("Solve explicit portion in fixed timestep mode. NOTE: This is "
                          "not recommended except for code comparison")
                     .withDefault(false)),
      order((*options)["order"].doc("Order of internal step").withDefault(4)),
      adap_method(
          (*options)["adap_method"]
              .doc("Set timestep adaptivity function: pid, pi, i, explicit_gustafsson,  "
                   "implicit_gustafsson, imex_gustafsson.")
              .withDefault(AdapMethod::PID)),
      abstol((*options)["atol"].doc("Absolute tolerance").withDefault(1.0e-12)),
      reltol((*options)["rtol"].doc("Relative tolerance").withDefault(1.0e-5)),
      suncontext(createSUNContext(BoutComm::get())) {
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

ArkodeMRISolver::~ArkodeMRISolver() {
  N_VDestroy(uvec);
  ARKodeFree(&arkode_mem);
  SUNLinSolFree(sun_solver);
  SUNNonlinSolFree(nonlinear_solver);

#if SUNDIALS_CONTROLLER_SUPPORT
  SUNAdaptController_Destroy(controller);
#endif
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int ArkodeMRISolver::init() {
  TRACE("Initialising ARKODE solver");

  Solver::init();

  output.write("Initialising SUNDIALS' ARKODE solver\n");

  // Calculate number of variables (in generic_solver)
  const int local_N = getLocalN();

  // Get total problem size
  int neq;
  if (bout::globals::mpi->MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    throw BoutException("Allreduce localN -> GlobalN failed!\n");
  }

  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n", n3Dvars(),
               n2Dvars(), neq, local_N);

  // Allocate memory
  uvec = callWithSUNContext(N_VNew_Parallel, suncontext, BoutComm::get(), local_N, neq);
  if (uvec == nullptr) {
    throw BoutException("SUNDIALS memory allocation failed\n");
  }

  // Put the variables into uvec
  save_vars(N_VGetArrayPointer(uvec));

  switch (treatment) {
  case Treatment::ImEx:
    arkode_mem = callWithSUNContext(ARKStepCreate, suncontext, arkode_rhs_s_explicit,
                                    arkode_rhs_s_implicit, simtime, uvec);
    break;
  case Treatment::Explicit:
    arkode_mem =
        callWithSUNContext(ARKStepCreate, suncontext, arkode_s_rhs, nullptr, simtime, uvec);
    break;
  case Treatment::Implicit:
    arkode_mem =
        callWithSUNContext(ARKStepCreate, suncontext, nullptr, arkode_s_rhs, simtime, uvec);
    break;
  default:
    throw BoutException("Invalid treatment: {}\n", toString(treatment));
  }
  if (arkode_mem == nullptr) {
    throw BoutException("ARKStepCreate failed\n");
  }

  switch (treatment) {
  case Treatment::ImEx:
    output_info.write("\tUsing ARKode ImEx solver \n");
    if (ARKStepSetImEx(arkode_mem) != ARK_SUCCESS) {
      throw BoutException("ARKodeSetImEx failed\n");
    }
    break;
  case Treatment::Explicit:
    output_info.write("\tUsing ARKode Explicit solver \n");
    if (ARKStepSetExplicit(arkode_mem) != ARK_SUCCESS) {
      throw BoutException("ARKodeSetExplicit failed\n");
    }
    break;
  case Treatment::Implicit:
    output_info.write("\tUsing ARKode Implicit solver \n");
    if (ARKStepSetImplicit(arkode_mem) != ARK_SUCCESS) {
      throw BoutException("ARKodeSetImplicit failed\n");
    }
    break;
  default:
    throw BoutException("Invalid treatment: {}\n", toString(treatment));
  }

  switch (inner_treatment) {
  case Treatment::ImEx:
    arkode_mem = callWithSUNContext(ARKStepCreate, suncontext, arkode_rhs_f_explicit,
                                    arkode_rhs_f_implicit, simtime, uvec);
    break;
  case Treatment::Explicit:
    arkode_mem =
        callWithSUNContext(ARKStepCreate, suncontext, arkode_f_rhs, nullptr, simtime, uvec);
    break;
  case Treatment::Implicit:
    arkode_mem =
        callWithSUNContext(ARKStepCreate, suncontext, nullptr, arkode_f_rhs, simtime, uvec);
    break;
  default:
    throw BoutException("Invalid treatment: {}\n", toString(treatment));
  }
  if (arkode_mem == nullptr) {
    throw BoutException("ARKStepCreate failed\n");
  }

  switch (inner_treatment) {
  case Treatment::ImEx:
    output_info.write("\tUsing ARKode ImEx solver \n");
    if (ARKStepSetImEx(arkode_mem) != ARK_SUCCESS) {
      throw BoutException("ARKodeSetImEx failed\n");
    }
    break;
  case Treatment::Explicit:
    output_info.write("\tUsing ARKode Explicit solver \n");
    if (ARKStepSetExplicit(arkode_mem) != ARK_SUCCESS) {
      throw BoutException("ARKodeSetExplicit failed\n");
    }
    break;
  case Treatment::Implicit:
    output_info.write("\tUsing ARKode Implicit solver \n");
    if (ARKStepSetImplicit(arkode_mem) != ARK_SUCCESS) {
      throw BoutException("ARKodeSetImplicit failed\n");
    }
    break;
  default:
    throw BoutException("Invalid treatment: {}\n", toString(treatment));
  }

  // For callbacks, need pointer to solver object
  if (ARKodeSetUserData(arkode_mem, this) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetUserData failed\n");
  }

  if (ARKodeSetLinear(arkode_mem, set_linear) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetLinear failed\n");
  }

  if (fixed_step) {
    // If not given, default to adaptive timestepping
    const auto fixed_timestep = (*options)["timestep"].withDefault(0.0);
    if (ARKodeSetFixedStep(arkode_mem, fixed_timestep) != ARK_SUCCESS) {
      throw BoutException("ARKodeSetFixedStep failed\n");
    }
  }

  if (ARKodeSetOrder(arkode_mem, order) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetOrder failed\n");
  }

#if SUNDIALS_CONTROLLER_SUPPORT
  switch (adap_method) {
  case AdapMethod::PID:
    controller = SUNAdaptController_PID(suncontext);
    break;
  case AdapMethod::PI:
    controller = SUNAdaptController_PI(suncontext);
    break;
  case AdapMethod::I:
    controller = SUNAdaptController_I(suncontext);
    break;
  case AdapMethod::Explicit_Gustafsson:
    controller = SUNAdaptController_ExpGus(suncontext);
    break;
  case AdapMethod::Implicit_Gustafsson:
    controller = SUNAdaptController_ImpGus(suncontext);
    break;
  case AdapMethod::ImEx_Gustafsson:
    controller = SUNAdaptController_ImExGus(suncontext);
    break;
  default:
    throw BoutException("Invalid adap_method\n");
  }

  if (ARKodeSetAdaptController(arkode_mem, controller) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetAdaptController failed\n");
  }

  if (ARKodeSetAdaptivityAdjustment(arkode_mem, 0) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetAdaptivityAdjustment failed\n");
  }
#else
  int adap_method_int;
  // Could cast to underlying integer, but this is more explicit
  switch (adap_method) {
  case AdapMethod::PID:
    adap_method_int = 0;
    break;
  case AdapMethod::PI:
    adap_method_int = 1;
    break;
  case AdapMethod::I:
    adap_method_int = 2;
    break;
  case AdapMethod::Explicit_Gustafsson:
    adap_method_int = 3;
    break;
  case AdapMethod::Implicit_Gustafsson:
    adap_method_int = 4;
    break;
  case AdapMethod::ImEx_Gustafsson:
    adap_method_int = 5;
    break;
  default:
    throw BoutException("Invalid adap_method\n");
  }

  if (ARKodeSetAdaptivityMethod(arkode_mem, adap_method_int, 1, 1, nullptr)
      != ARK_SUCCESS) {
    throw BoutException("ARKodeSetAdaptivityMethod failed\n");
  }
#endif

  if (use_vector_abstol) {
    std::vector<BoutReal> f2dtols;
    f2dtols.reserve(f2d.size());
    std::transform(begin(f2d), end(f2d), std::back_inserter(f2dtols),
                   [abstol = abstol](const VarStr<Field2D>& f2) {
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
                   [abstol = abstol](const VarStr<Field3D>& f3) {
                     return Options::root()[f3.name]["atol"].withDefault(abstol);
                   });

    N_Vector abstolvec = N_VClone(uvec);
    if (abstolvec == nullptr) {
      throw BoutException("SUNDIALS memory allocation (abstol vector) failed\n");
    }

    set_abstol_values(N_VGetArrayPointer(abstolvec), f2dtols, f3dtols);

    if (ARKodeSVtolerances(arkode_mem, reltol, abstolvec) != ARK_SUCCESS) {
      throw BoutException("ARKodeSVtolerances failed\n");
    }

    N_VDestroy(abstolvec);
  } else {
    if (ARKodeSStolerances(arkode_mem, reltol, abstol) != ARK_SUCCESS) {
      throw BoutException("ARKodeSStolerances failed\n");
    }
  }

  if (ARKodeSetMaxNumSteps(arkode_mem, mxsteps) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetMaxNumSteps failed\n");
  }

  if (treatment == Treatment::ImEx or treatment == Treatment::Implicit) {
    {
      output.write("\tUsing Newton iteration\n");

      const auto prectype =
          use_precon ? (rightprec ? SUN_PREC_RIGHT : SUN_PREC_LEFT) : SUN_PREC_NONE;
      sun_solver = callWithSUNContext(SUNLinSol_SPGMR, suncontext, uvec, prectype, maxl);
      if (sun_solver == nullptr) {
        throw BoutException("Creating SUNDIALS linear solver failed\n");
      }
      if (ARKodeSetLinearSolver(arkode_mem, sun_solver, nullptr) != ARKLS_SUCCESS) {
        throw BoutException("ARKodeSetLinearSolver failed\n");
      }

      /// Set Preconditioner
      if (use_precon) {
        if (hasPreconditioner()) {
          output.write("\tUsing user-supplied preconditioner\n");

          if (ARKodeSetPreconditioner(arkode_mem, nullptr, arkode_s_pre)
              != ARKLS_SUCCESS) {
            throw BoutException("ARKodeSetPreconditioner failed\n");
          }
        } else {
          output.write("\tUsing BBD preconditioner\n");

          /// Get options
          // Compute band_width_default from actually added fields, to allow for multiple
          // Mesh objects
          //
          // Previous implementation was equivalent to:
          //   int MXSUB = mesh->xend - mesh->xstart + 1;
          //   int band_width_default = n3Dvars()*(MXSUB+2);
          const int band_width_default = std::accumulate(
              begin(f3d), end(f3d), 0, [](int acc, const VarStr<Field3D>& fvar) {
                Mesh* localmesh = fvar.var->getMesh();
                return acc + localmesh->xend - localmesh->xstart + 3;
              });

          const auto mudq = (*options)["mudq"]
                                .doc("Upper half-bandwidth to be used in the difference "
                                     "quotient Jacobian approximation")
                                .withDefault(band_width_default);
          const auto mldq = (*options)["mldq"]
                                .doc("Lower half-bandwidth to be used in the difference "
                                     "quotient Jacobian approximation")
                                .withDefault(band_width_default);
          const auto mukeep = (*options)["mukeep"]
                                  .doc("Upper half-bandwidth of the retained banded "
                                       "approximate Jacobian block")
                                  .withDefault(n3Dvars() + n2Dvars());
          const auto mlkeep = (*options)["mlkeep"]
                                  .doc("Lower half-bandwidth of the retained banded "
                                       "approximate Jacobian block")
                                  .withDefault(n3Dvars() + n2Dvars());

          if (ARKBBDPrecInit(arkode_mem, local_N, mudq, mldq, mukeep, mlkeep, 0,
                             arkode_s_bbd_rhs, nullptr)
              != ARKLS_SUCCESS) {
            throw BoutException("ARKBBDPrecInit failed\n");
          }
        }
      } else {
        // Not using preconditioning
        output.write("\tNo preconditioning\n");
      }
    }

    /// Set Jacobian-vector multiplication function
    output.write("\tUsing difference quotient approximation for Jacobian\n");
  }

  return 0;
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int ArkodeMRISolver::run() {
  TRACE("ArkodeMRISolver::run()");

  if (!initialised) {
    throw BoutException("ArkodeMRISolver not initialised\n");
  }

  for (int i = 0; i < getNumberOutputSteps(); i++) {

    /// Run the solver for one output timestep
    simtime = run(simtime + getOutputTimestep());

    /// Check if the run succeeded
    if (simtime < 0.0) {
      // Step failed
      output.write("Timestep failed. Aborting\n");

      throw BoutException("ARKode timestep failed\n");
    }

    // Get additional diagnostics
    long int temp_long_int, temp_long_int2;
    ARKodeGetNumSteps(arkode_mem, &temp_long_int);
    nsteps = int(temp_long_int);
    ARKStepGetNumRhsEvals(arkode_mem, &temp_long_int, &temp_long_int2); //Change after the release
    nfe_evals = int(temp_long_int);
    nfi_evals = int(temp_long_int2);
    if (treatment == Treatment::ImEx or treatment == Treatment::Implicit) {
      ARKodeGetNumNonlinSolvIters(arkode_mem, &temp_long_int);
      nniters = int(temp_long_int);
      ARKodeGetNumPrecEvals(arkode_mem, &temp_long_int);
      npevals = int(temp_long_int);
      ARKodeGetNumLinIters(arkode_mem, &temp_long_int);
      nliters = int(temp_long_int);
    }

    if (diagnose) {
      output.write("\nARKODE: nsteps {:d}, nfe_evals {:d}, nfi_evals {:d}, nniters {:d}, "
                   "npevals {:d}, nliters {:d}\n",
                   nsteps, nfe_evals, nfi_evals, nniters, npevals, nliters);
      if (treatment == Treatment::ImEx or treatment == Treatment::Implicit) {
        output.write("    -> Newton iterations per step: {:e}\n",
                     static_cast<BoutReal>(nniters) / static_cast<BoutReal>(nsteps));
        output.write("    -> Linear iterations per Newton iteration: {:e}\n",
                     static_cast<BoutReal>(nliters) / static_cast<BoutReal>(nniters));
        output.write("    -> Preconditioner evaluations per Newton: {:e}\n",
                     static_cast<BoutReal>(npevals) / static_cast<BoutReal>(nniters));
      }
    }

    if (call_monitors(simtime, i, getNumberOutputSteps())) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

BoutReal ArkodeMRISolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run({:e})", tout);

  bout::globals::mpi->MPI_Barrier(BoutComm::get());

  pre_Wtime_s = 0.0;
  pre_ncalls_s = 0;

  int flag;
  if (!monitor_timestep) {
    // Run in normal mode
    flag = ARKodeEvolve(arkode_mem, tout, uvec, &simtime, ARK_NORMAL);
  } else {
    // Run in single step mode, to call timestep monitors
    BoutReal internal_time;
    ARKodeGetCurrentTime(arkode_mem, &internal_time);
    while (internal_time < tout) {
      // Run another step
      const BoutReal last_time = internal_time;
      flag = ARKodeEvolve(arkode_mem, tout, uvec, &internal_time, ARK_ONE_STEP);

      if (flag != ARK_SUCCESS) {
        output_error.write("ERROR ARKODE solve failed at t = {:e}, flag = {:d}\n",
                           internal_time, flag);
        return -1.0;
      }

      // Call timestep monitor
      call_timestep_monitors(internal_time, internal_time - last_time);
    }
    // Get output at the desired time
    flag = ARKodeGetDky(arkode_mem, tout, 0, uvec);
    simtime = tout;
  }

  // Copy variables
  load_vars(N_VGetArrayPointer(uvec));
  // Call rhs function to get extra variables at this time
  run_rhs(simtime);
  // run_diffusive(simtime);
  if (flag != ARK_SUCCESS) {
    output_error.write("ERROR ARKODE solve failed at t = {:e}, flag = {:d}\n", simtime,
                       flag);
    return -1.0;
  }

  return simtime;
}

/**************************************************************************
 * Explicit RHS function du = F^s_E(t, u)
 **************************************************************************/

void ArkodeMRISolver::rhs_se(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeMRISolver::rhs_e({:e})", t);

  // Load state from udata
  load_vars(udata);

  // Get the current timestep
  // Note: ARKodeGetCurrentStep updated too late in older versions
  ARKodeGetLastStep(arkode_mem, &hcur);

  // Call RHS function
  run_convective(t);

  // Save derivatives to dudata
  save_derivs(dudata);
}

/**************************************************************************
 *   Implicit RHS function du = F^s_I(t, u)
 **************************************************************************/

void ArkodeMRISolver::rhs_si(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeMRISolver::rhs_si({:e})", t);

  load_vars(udata);
  ARKodeGetLastStep(arkode_mem, &hcur);
  // Call Implicit RHS function
  run_diffusive(t);
  save_derivs(dudata);
}

/**************************************************************************
 * Explicit RHS function du = F^f_E(t, u)
 **************************************************************************/

void ArkodeMRISolver::rhs_fe(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeMRISolver::rhs_e({:e})", t);

  // Load state from udata
  load_vars(udata);

  // Get the current timestep
  // Note: ARKodeGetCurrentStep updated too late in older versions
  ARKodeGetLastStep(arkode_mem, &hcur);

  // Call RHS function
  run_convective(t);

  // Save derivatives to dudata
  save_derivs(dudata);
}

/**************************************************************************
 *   Implicit RHS function du = F^f_I(t, u)
 **************************************************************************/

void ArkodeMRISolver::rhs_fi(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeMRISolver::rhs_si({:e})", t);

  load_vars(udata);
  ARKodeGetLastStep(arkode_mem, &hcur);
  // Call Implicit RHS function
  run_diffusive(t);
  save_derivs(dudata);
}

/**************************************************************************
 *   Full  RHS function du = F(t, u)
 **************************************************************************/
void ArkodeMRISolver::rhs(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeMRISolver::rhs({:e})", t);

  load_vars(udata);
  ARKodeGetLastStep(arkode_mem, &hcur);
  // Call Implicit RHS function
  run_rhs(t);
  save_derivs(dudata);
}

/**************************************************************************
 * Preconditioner functions
 **************************************************************************/

void ArkodeMRISolver::pre_s(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal* udata,
                       BoutReal* rvec, BoutReal* zvec) {
  TRACE("Running preconditioner: ArkodeMRISolver::pre({:e})", t);

  const BoutReal tstart = bout::globals::mpi->MPI_Wtime();

  if (!hasPreconditioner()) {
    // Identity (but should never happen)
    const auto length = N_VGetLocalLength_Parallel(uvec);
    std::copy(rvec, rvec + length, zvec);
    return;
  }

  // Load state from udata (as with res function)
  load_vars(udata);

  // Load vector to be inverted into F_vars
  load_derivs(rvec);

  runPreconditioner(t, gamma, delta);

  // Save the solution from F_vars
  save_derivs(zvec);

  pre_Wtime_s += bout::globals::mpi->MPI_Wtime() - tstart;
  pre_ncalls_s++;
}

void ArkodeMRISolver::pre_f(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal* udata,
                       BoutReal* rvec, BoutReal* zvec) {
  TRACE("Running preconditioner: ArkodeMRISolver::pre({:e})", t);

  const BoutReal tstart = bout::globals::mpi->MPI_Wtime();

  if (!hasPreconditioner()) {
    // Identity (but should never happen)
    const auto length = N_VGetLocalLength_Parallel(uvec);
    std::copy(rvec, rvec + length, zvec);
    return;
  }

  // Load state from udata (as with res function)
  load_vars(udata);

  // Load vector to be inverted into F_vars
  load_derivs(rvec);

  runPreconditioner(t, gamma, delta);

  // Save the solution from F_vars
  save_derivs(zvec);

  pre_Wtime_s += bout::globals::mpi->MPI_Wtime() - tstart;
  pre_ncalls_s++;
}

/**************************************************************************
 * Jacobian-vector multiplication functions
 **************************************************************************/

void ArkodeMRISolver::jac_s(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata) {
  TRACE("Running Jacobian: ArkodeMRISolver::jac({:e})", t);

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

void ArkodeMRISolver::jac_f(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata) {
  TRACE("Running Jacobian: ArkodeMRISolver::jac({:e})", t);

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
 * ARKODE explicit RHS functions
 **************************************************************************/

// NOLINTBEGIN(readability-identifier-length)
namespace {
int arkode_rhs_s_explicit(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = N_VGetArrayPointer(u);
  BoutReal* dudata = N_VGetArrayPointer(du);

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs_se(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

int arkode_rhs_s_implicit(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = N_VGetArrayPointer(u);
  BoutReal* dudata = N_VGetArrayPointer(du);

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs_si(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

int arkode_rhs_f_explicit(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = N_VGetArrayPointer(u);
  BoutReal* dudata = N_VGetArrayPointer(du);

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs_fe(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

int arkode_rhs_f_implicit(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = N_VGetArrayPointer(u);
  BoutReal* dudata = N_VGetArrayPointer(du);

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs_fi(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

int arkode_s_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = N_VGetArrayPointer(u);
  BoutReal* dudata = N_VGetArrayPointer(du);

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs_s(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

int arkode_f_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = N_VGetArrayPointer(u);
  BoutReal* dudata = N_VGetArrayPointer(du);

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs_f(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

int arkode_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data) {

  BoutReal* udata = N_VGetArrayPointer(u);
  BoutReal* dudata = N_VGetArrayPointer(du);

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs(t, udata, dudata);
  } catch (BoutRhsFail& error) {
    return 1;
  }
  return 0;
}

/// RHS function for BBD preconditioner
int arkode_s_bbd_rhs(sunindextype UNUSED(Nlocal), BoutReal t, N_Vector u, N_Vector du,
                   void* user_data) {
  return arkode_rhs_s_implicit(t, u, du, user_data);
}

int arkode_f_bbd_rhs(sunindextype UNUSED(Nlocal), BoutReal t, N_Vector u, N_Vector du,
                   void* user_data) {
  return arkode_rhs_f_implicit(t, u, du, user_data);
}

/// Preconditioner function
int arkode_s_pre(BoutReal t, N_Vector yy, N_Vector UNUSED(yp), N_Vector rvec, N_Vector zvec,
               BoutReal gamma, BoutReal delta, int UNUSED(lr), void* user_data) {
  BoutReal* udata = N_VGetArrayPointer(yy);
  BoutReal* rdata = N_VGetArrayPointer(rvec);
  BoutReal* zdata = N_VGetArrayPointer(zvec);

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  // Calculate residuals
  s->pre_s(t, gamma, delta, udata, rdata, zdata);

  return 0;
}

int arkode_f_pre(BoutReal t, N_Vector yy, N_Vector UNUSED(yp), N_Vector rvec, N_Vector zvec,
               BoutReal gamma, BoutReal delta, int UNUSED(lr), void* user_data) {
  BoutReal* udata = N_VGetArrayPointer(yy);
  BoutReal* rdata = N_VGetArrayPointer(rvec);
  BoutReal* zdata = N_VGetArrayPointer(zvec);

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  // Calculate residuals
  s->pre_f(t, gamma, delta, udata, rdata, zdata);

  return 0;
}

/// Jacobian-vector multiplication function
int arkode_s_jac(N_Vector v, N_Vector Jv, BoutReal t, N_Vector y, N_Vector UNUSED(fy),
               void* user_data, N_Vector UNUSED(tmp)) {
  BoutReal* ydata = N_VGetArrayPointer(y);   ///< System state
  BoutReal* vdata = N_VGetArrayPointer(v);   ///< Input vector
  BoutReal* Jvdata = N_VGetArrayPointer(Jv); ///< Jacobian*vector output

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  s->jac_s(t, ydata, vdata, Jvdata);

  return 0;
}

int arkode_f_jac(N_Vector v, N_Vector Jv, BoutReal t, N_Vector y, N_Vector UNUSED(fy),
               void* user_data, N_Vector UNUSED(tmp)) {
  BoutReal* ydata = N_VGetArrayPointer(y);   ///< System state
  BoutReal* vdata = N_VGetArrayPointer(v);   ///< Input vector
  BoutReal* Jvdata = N_VGetArrayPointer(Jv); ///< Jacobian*vector output

  auto* s = static_cast<ArkodeMRISolver*>(user_data);

  s->jac_f(t, ydata, vdata, Jvdata);

  return 0;
}
} // namespace
// NOLINTEND(readability-identifier-length)

/**************************************************************************
 * vector abstol functions
 **************************************************************************/

void ArkodeMRISolver::set_abstol_values(BoutReal* abstolvec_data,
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

void ArkodeMRISolver::loop_abstol_values_op(Ind2D UNUSED(i2d), BoutReal* abstolvec_data,
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

#endif
