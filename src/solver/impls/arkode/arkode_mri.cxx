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
                    .withDefault(MRI_Treatment::ImEx)),
      inner_treatment((*options)["inner_treatment"]
                    .doc("Use default capability (imex) or provide a specific inner_treatment: "
                         "implicit or explicit")
                    .withDefault(MRI_Treatment::ImEx)),
      set_linear(
          (*options)["set_linear"]
              .doc("Use linear implicit solver (only evaluates jacobian inversion once)")
              .withDefault(false)),
      inner_set_linear(
          (*options)["inner_set_linear"]
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
              .withDefault(MRI_AdapMethod::PID)),
      abstol((*options)["atol"].doc("Absolute tolerance").withDefault(1.0e-12)),
      reltol((*options)["rtol"].doc("Relative tolerance").withDefault(1.0e-5)),
      use_vector_abstol((*options)["use_vector_abstol"]
                            .doc("Use separate absolute tolerance for each field")
                            .withDefault(false)),
      use_precon((*options)["use_precon"]
                     .doc("Use user-supplied preconditioner function")
                     .withDefault(false)),
      inner_use_precon((*options)["inner_use_precon"]
                     .doc("Use user-supplied preconditioner function")
                     .withDefault(false)),
      maxl(
          (*options)["maxl"].doc("Number of Krylov basis vectors to use").withDefault(0)),
      inner_maxl(
          (*options)["inner_maxl"].doc("Number of Krylov basis vectors to use").withDefault(0)),
      rightprec((*options)["rightprec"]
                    .doc("Use right preconditioning instead of left preconditioning")
                    .withDefault(false)),
      suncontext(createSUNContext(BoutComm::get())) {
  has_constraints = false; // This solver doesn't have constraints

  // Add diagnostics to output
  add_int_diagnostic(nsteps, "arkode_nsteps", "Cumulative number of internal steps");
  add_int_diagnostic(inner_nsteps, "arkode_inner_nsteps", "Cumulative number of inner internal steps");
  add_int_diagnostic(nfe_evals, "arkode_nfe_evals",
                     "No. of calls to fe (explicit portion of the right-hand-side "
                     "function) function");
  add_int_diagnostic(inner_nfe_evals, "arkode_inner_nfe_evals",
                     "No. of calls to fe (explicit portion of the inner right-hand-side "
                     "function) function");
  add_int_diagnostic(nfi_evals, "arkode_nfi_evals",
                     "No. of calls to fi (implicit portion of the right-hand-side "
                     "function) function");
  add_int_diagnostic(inner_nfi_evals, "arkode_inner_nfi_evals",
                     "No. of calls to fi (implicit portion of inner the right-hand-side "
                     "function) function");
  add_int_diagnostic(nniters, "arkode_nniters", "No. of nonlinear solver iterations");
  add_int_diagnostic(inner_nniters, "arkode_inner_nniters", "No. of inner nonlinear solver iterations");
  add_int_diagnostic(npevals, "arkode_npevals", "No. of preconditioner evaluations");
  add_int_diagnostic(inner_npevals, "arkode_inner_npevals", "No. of inner preconditioner evaluations");
  add_int_diagnostic(nliters, "arkode_nliters", "No. of linear iterations");
  add_int_diagnostic(inner_nliters, "arkode_inner_nliters", "No. of inner linear iterations");
}

ArkodeMRISolver::~ArkodeMRISolver() {
  N_VDestroy(uvec);
  ARKodeFree(&arkode_mem);
  SUNLinSolFree(sun_solver);
  SUNNonlinSolFree(nonlinear_solver);

  SUNAdaptController_Destroy(controller);
  SUNAdaptController_Destroy(inner_controller);
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int ArkodeMRISolver::init() {
  TRACE("Initialising ARKODE MRI solver");

  Solver::init();

  output.write("Initialising SUNDIALS' ARKODE MRI solver\n");

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

  switch (inner_treatment) {
  case MRI_Treatment::ImEx:
    inner_arkode_mem = callWithSUNContext(ARKStepCreate, suncontext, arkode_rhs_f_explicit,
                                    arkode_rhs_f_implicit, simtime, uvec);
    output_info.write("\tUsing ARKode ImEx inner solver \n");
    break;
  case MRI_Treatment::Explicit:
    inner_arkode_mem =
        callWithSUNContext(ARKStepCreate, suncontext, arkode_f_rhs, nullptr, simtime, uvec);
    output_info.write("\tUsing ARKode Explicit inner solver \n");
    break;
  case MRI_Treatment::Implicit:
    inner_arkode_mem =
        callWithSUNContext(ARKStepCreate, suncontext, nullptr, arkode_f_rhs, simtime, uvec);
    output_info.write("\tUsing ARKode Implicit inner solver \n");
    break;
  default:
    throw BoutException("Invalid inner_treatment: {}\n", toString(inner_treatment));
  }
  if (inner_arkode_mem == nullptr) {
    throw BoutException("ARKStepCreate failed\n");
  }

  // For callbacks, need pointer to solver object
  if (ARKodeSetUserData(inner_arkode_mem, this) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetUserData failed\n");
  }

  if(inner_treatment != MRI_Treatment::Explicit)
    if (ARKodeSetLinear(inner_arkode_mem, set_linear) != ARK_SUCCESS) {
      throw BoutException("ARKodeSetLinear failed\n");
    }

  if (fixed_step) {
    // If not given, default to adaptive timestepping
    const auto inner_fixed_timestep = (*options)["inner_timestep"].withDefault(0.0);
    if (ARKodeSetFixedStep(inner_arkode_mem, inner_fixed_timestep) != ARK_SUCCESS) {
      throw BoutException("ARKodeSetFixedStep failed\n");
    }
  }

  if (ARKodeSetOrder(inner_arkode_mem, order) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetOrder failed\n");
  }

  if (ARKStepCreateMRIStepInnerStepper(inner_arkode_mem, &inner_stepper) != ARK_SUCCESS) {
    throw BoutException("ARKStepCreateMRIStepInnerStepper failed\n");
  }

  // Initialize the slow integrator. Specify the explicit slow right-hand side
  // function in y'=fe(t,y)+fi(t,y)+ff(t,y), the inital time T0, the
  // initial dependent variable vector y, and the fast integrator.

  switch (treatment) {
  case MRI_Treatment::ImEx:
    arkode_mem = callWithSUNContext(MRIStepCreate, suncontext, arkode_rhs_s_explicit, arkode_rhs_s_implicit, 
                                    simtime, uvec, inner_stepper);
    output_info.write("\tUsing ARKode ImEx solver \n");
    break;
  case MRI_Treatment::Explicit:
    arkode_mem = callWithSUNContext(MRIStepCreate, suncontext, arkode_s_rhs, nullptr, 
                                    simtime, uvec, inner_stepper);
    output_info.write("\tUsing ARKode Explicit solver \n");
    break;
  case MRI_Treatment::Implicit:
    arkode_mem = callWithSUNContext(MRIStepCreate, suncontext, nullptr, arkode_s_rhs,
                                    simtime, uvec, inner_stepper);
    output_info.write("\tUsing ARKode Implicit solver \n");
    break;
  default:
    throw BoutException("Invalid treatment: {}\n", toString(treatment));
  }
  if (arkode_mem == nullptr) {
    throw BoutException("MRIStepCreate failed\n");
  }

  // For callbacks, need pointer to solver object
  if (ARKodeSetUserData(arkode_mem, this) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetUserData failed\n");
  }

  if(treatment != MRI_Treatment::Explicit)
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

  switch (adap_method) {
  case MRI_AdapMethod::PID:
    controller = SUNAdaptController_PID(suncontext);
    inner_controller = SUNAdaptController_PID(suncontext);
    break;
  case MRI_AdapMethod::PI:
    controller = SUNAdaptController_PI(suncontext);
    inner_controller = SUNAdaptController_PI(suncontext);
    break;
  case MRI_AdapMethod::I:
    controller = SUNAdaptController_I(suncontext);
    inner_controller = SUNAdaptController_I(suncontext);
    break;
  case MRI_AdapMethod::Explicit_Gustafsson:
    controller = SUNAdaptController_ExpGus(suncontext);
    inner_controller = SUNAdaptController_ExpGus(suncontext);
    break;
  case MRI_AdapMethod::Implicit_Gustafsson:
    controller = SUNAdaptController_ImpGus(suncontext);
    inner_controller = SUNAdaptController_ImpGus(suncontext);
    break;
  case MRI_AdapMethod::ImEx_Gustafsson:
    controller = SUNAdaptController_ImExGus(suncontext);
    inner_controller = SUNAdaptController_ImExGus(suncontext);
    break;
  default:
    throw BoutException("Invalid adap_method\n");
  }

  // if (ARKodeSetAdaptController(arkode_mem, controller) != ARK_SUCCESS) {
  //   throw BoutException("ARKodeSetAdaptController failed\n");
  // }

  // if (ARKodeSetAdaptivityAdjustment(arkode_mem, 0) != ARK_SUCCESS) {
  //   throw BoutException("ARKodeSetAdaptivityAdjustment failed\n");
  // }

  if (ARKodeSetFixedStep(arkode_mem, 0.0001) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetAdaptController failed\n");
  }

  if (ARKodeSetAdaptController(inner_arkode_mem, controller) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetAdaptController failed\n");
  }

  if (ARKodeSetAdaptivityAdjustment(inner_arkode_mem, 0) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetAdaptivityAdjustment failed\n");
  }

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
    if (ARKodeSVtolerances(inner_arkode_mem, reltol, abstolvec) != ARK_SUCCESS) {
      throw BoutException("ARKodeSVtolerances failed\n");
    }

    N_VDestroy(abstolvec);
  } else {
    if (ARKodeSStolerances(arkode_mem, reltol, abstol) != ARK_SUCCESS) {
      throw BoutException("ARKodeSStolerances failed\n");
    }
    if (ARKodeSStolerances(inner_arkode_mem, reltol, abstol) != ARK_SUCCESS) {
      throw BoutException("ARKodeSStolerances failed\n");
    }
  }

  if (ARKodeSetMaxNumSteps(arkode_mem, mxsteps) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetMaxNumSteps failed\n");
  }
  if (ARKodeSetMaxNumSteps(inner_arkode_mem, mxsteps) != ARK_SUCCESS) {
    throw BoutException("ARKodeSetMaxNumSteps failed\n");
  }

  if (inner_treatment == MRI_Treatment::ImEx or inner_treatment == MRI_Treatment::Implicit) {
    {
      output.write("\tUsing Newton iteration\n");

      const auto prectype =
          inner_use_precon ? (rightprec ? SUN_PREC_RIGHT : SUN_PREC_LEFT) : SUN_PREC_NONE;
      inner_sun_solver = callWithSUNContext(SUNLinSol_SPGMR, suncontext, uvec, prectype, inner_maxl);
      if (inner_sun_solver == nullptr) {
        throw BoutException("Creating SUNDIALS linear solver failed\n");
      }
      if (ARKodeSetLinearSolver(inner_arkode_mem, inner_sun_solver, nullptr) != ARKLS_SUCCESS) {
        throw BoutException("ARKodeSetLinearSolver failed\n");
      }

      /// Set Preconditioner
      if (inner_use_precon) {
        if (hasPreconditioner()) {  // change to inner_hasPreconditioner when it is available
          output.write("\tUsing user-supplied preconditioner\n");

          if (ARKodeSetPreconditioner(inner_arkode_mem, nullptr, arkode_f_pre)
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

          if (ARKBBDPrecInit(inner_arkode_mem, local_N, mudq, mldq, mukeep, mlkeep, 0,
                             arkode_f_bbd_rhs, nullptr)
              != ARKLS_SUCCESS) {
            throw BoutException("ARKBBDPrecInit failed\n");
          }
        }
      } else {
        // Not using preconditioning
        output.write("\tNo inner preconditioning\n");
      }
    }

    /// Set Jacobian-vector multiplication function
    output.write("\tUsing difference quotient approximation for Jacobian\n");
  }


  if (treatment == MRI_Treatment::ImEx or treatment == MRI_Treatment::Implicit) {
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
    MRIStepGetNumRhsEvals(arkode_mem, &temp_long_int, &temp_long_int2); //Change after the release
    nfe_evals = int(temp_long_int);
    nfi_evals = int(temp_long_int2);
    if (treatment == MRI_Treatment::ImEx or treatment == MRI_Treatment::Implicit) {
      ARKodeGetNumNonlinSolvIters(arkode_mem, &temp_long_int);
      nniters = int(temp_long_int);
      ARKodeGetNumPrecEvals(arkode_mem, &temp_long_int);
      npevals = int(temp_long_int);
      ARKodeGetNumLinIters(arkode_mem, &temp_long_int);
      nliters = int(temp_long_int);
    }

    ARKodeGetNumSteps(inner_arkode_mem, &temp_long_int);
    inner_nsteps = int(temp_long_int);
    ARKStepGetNumRhsEvals(inner_arkode_mem, &temp_long_int, &temp_long_int2); //Change after the release
    inner_nfe_evals = int(temp_long_int);
    inner_nfi_evals = int(temp_long_int2);
    if (inner_treatment == MRI_Treatment::ImEx or inner_treatment == MRI_Treatment::Implicit) {
      ARKodeGetNumNonlinSolvIters(inner_arkode_mem, &temp_long_int);
      inner_nniters = int(temp_long_int);
      ARKodeGetNumPrecEvals(inner_arkode_mem, &temp_long_int);
      inner_npevals = int(temp_long_int);
      ARKodeGetNumLinIters(inner_arkode_mem, &temp_long_int);
      inner_nliters = int(temp_long_int);
    }

    if (diagnose) {
      output.write("\nARKODE: nsteps {:d}, nfe_evals {:d}, nfi_evals {:d}, nniters {:d}, "
                   "npevals {:d}, nliters {:d}\n",
                   nsteps, nfe_evals, nfi_evals, nniters, npevals, nliters);
      if (treatment == MRI_Treatment::ImEx or treatment == MRI_Treatment::Implicit) {
        output.write("    -> Newton iterations per step: {:e}\n",
                     static_cast<BoutReal>(nniters) / static_cast<BoutReal>(nsteps));
        output.write("    -> Linear iterations per Newton iteration: {:e}\n",
                     static_cast<BoutReal>(nliters) / static_cast<BoutReal>(nniters));
        output.write("    -> Preconditioner evaluations per Newton: {:e}\n",
                     static_cast<BoutReal>(npevals) / static_cast<BoutReal>(nniters));
      }

      output.write("\nARKODE Inner: inner_nsteps {:d}, inner_nfe_evals {:d}, inner_nfi_evals {:d}, inner_nniters {:d}, "
                   "inner_npevals {:d}, inner_nliters {:d}\n",
                   inner_nsteps, inner_nfe_evals, inner_nfi_evals, inner_nniters, inner_npevals, inner_nliters);
      if (inner_treatment == MRI_Treatment::ImEx or inner_treatment == MRI_Treatment::Implicit) {
        output.write("    -> Inner Newton iterations per step: {:e}\n",
                     static_cast<BoutReal>(inner_nniters) / static_cast<BoutReal>(inner_nsteps));
        output.write("    -> Inner Linear iterations per Newton iteration: {:e}\n",
                     static_cast<BoutReal>(inner_nliters) / static_cast<BoutReal>(inner_nniters));
        output.write("    -> Inner Preconditioner evaluations per Newton: {:e}\n",
                     static_cast<BoutReal>(inner_npevals) / static_cast<BoutReal>(inner_nniters));
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
  run_rhs_se(t);

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
  run_rhs_si(t);
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
  run_rhs_fe(t);

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
  run_rhs_fi(t);
  save_derivs(dudata);
}

/**************************************************************************
 * Slow RHS function du = F^s(t, u)
 **************************************************************************/

void ArkodeMRISolver::rhs_s(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeMRISolver::rhs_e({:e})", t);

  // Load state from udata
  load_vars(udata);

  // Get the current timestep
  // Note: ARKodeGetCurrentStep updated too late in older versions
  ARKodeGetLastStep(arkode_mem, &hcur);

  // Call RHS function
  // run_rhs_s(t);
  run_rhs_s(t);

  // Save derivatives to dudata
  save_derivs(dudata);
}

/**************************************************************************
 * Fast RHS function du = F^f(t, u)
 **************************************************************************/

void ArkodeMRISolver::rhs_f(BoutReal t, BoutReal* udata, BoutReal* dudata) {
  TRACE("Running RHS: ArkodeMRISolver::rhs_e({:e})", t);

  // Load state from udata
  load_vars(udata);

  // Get the current timestep
  // Note: ARKodeGetCurrentStep updated too late in older versions
  ARKodeGetLastStep(arkode_mem, &hcur);

  // Call RHS function
  // run_rhs_f(t);
  run_rhs_f(t);

  // Save derivatives to dudata
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
