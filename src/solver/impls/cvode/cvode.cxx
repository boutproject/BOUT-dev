/**************************************************************************
 * Interface to SUNDIALS CVODE
 *
 *
 **************************************************************************
 * Copyright 2010-2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#include "bout/build_defines.hxx"
#include <cstddef>

#include "cvode.hxx"

#if BOUT_HAS_CVODE

#include "../../sundials_nvector_interface.hxx"

#include "bout/bout_enum_class.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/build_config.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/mpi_wrapper.hxx"
#include "bout/msg_stack.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"
#include "bout/petsc_preconditioner.hxx"
#include "bout/petsclib.hxx"
#include "bout/region.hxx"
#include "bout/solver.hxx"
#include "bout/sundials_backports.hxx"
#include "bout/unused.hxx"

#include "fmt/core.h"
#include "fmt/format.h"

#include <cvode/cvode.h>
#include <cvode/cvode_bbdpre.h>
#include <cvode/cvode_ls.h>

#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>

#include <algorithm>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

BOUT_ENUM_CLASS_NS(bout, positivity_constraint, none, positive, non_negative, negative,
                   non_positive);

BOUT_ENUM_CLASS_NS(bout, linear_solver, gmres, fgmres, tfqmr, bcgs);

// NOLINTBEGIN(readability-identifier-length)
namespace {
int cvode_linear_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data);
int cvode_nonlinear_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data);
int cvode_bbd_rhs(sunindextype Nlocal, BoutReal t, N_Vector u, N_Vector du,
                  void* user_data);

int cvode_pre(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rvec, N_Vector zvec,
              BoutReal gamma, BoutReal delta, int lr, void* user_data);

int cvode_jac(N_Vector v, N_Vector Jv, BoutReal t, N_Vector y, N_Vector fy,
              void* user_data, N_Vector tmp);

#if SUNDIALS_VERSION_MAJOR >= 7
// Used by solver:cvode_failure_counters
/// Wrap a Krylov solver in a SUNLinearSolver that delegates every operation and
/// reports the return flag of each solve to CvodeSolver::countKrylovFlag. Takes
/// ownership of \p inner.
SUNLinearSolver makeCountingLinearSolver(SUNLinearSolver inner, CvodeSolver* owner,
                                         SUNContext ctx);
/// Solve function of the Newton solver: calls SUNNonlinSolSolve_Newton and reports
/// a failed solve to CvodeSolver::countNewtonSolveFailure
int countingNewtonSolve(SUNNonlinearSolver nls, N_Vector y0, N_Vector ycor, N_Vector w,
                        sunrealtype tol, sunbooleantype call_lsetup, void* mem);
#endif
} // namespace
// NOLINTEND(readability-identifier-length)

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
                    .doc("Maximum order of method to use. <= 0 means default limit.")
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
          (*options)["max_timestep"].doc("Maximum time step size").withDefault(0.0)),
      min_timestep(
          (*options)["min_timestep"].doc("Minimum time step size").withDefault(0.0)),
      start_timestep((*options)["start_timestep"]
                         .doc("Starting time step. = 0 then chosen by CVODE.")
                         .withDefault(0.0)),
      mxorder((*options)["mxorder"].doc("Maximum order").withDefault(-1)),
      max_nonlinear_iterations(
          (*options)["max_nonlinear_iterations"]
              .doc("Maximum number of nonlinear iterations allowed by CVODE before "
                   "reducing timestep.")
              .withDefault(3)),
      lsetup_frequency(
          (*options)["cvode_lsetup_frequency"]
              .doc("Linear solver setup frequency (CVodeSetLSetupFrequency). "
                   "0 uses the SUNDIALS default.")
              .withDefault(0)),
      jac_eval_frequency((*options)["cvode_jac_eval_frequency"]
                             .doc("Jacobian/preconditioner evaluation frequency "
                                  "(CVodeSetJacEvalFrequency). "
                                  "0 uses the SUNDIALS default.")
                             .withDefault(0)),
      apply_positivity_constraints(
          (*options)["apply_positivity_constraints"]
              .doc("Use CVODE function CVodeSetConstraints to constrain variables - the "
                   "constraint "
                   "to be applied is set by the positivity_constraint option in the "
                   "subsection for "
                   "each variable")
              .withDefault(false)),
      maxl((*options)["maxl"].doc("Maximum number of linear iterations").withDefault(5)),
      rightprec((*options)["rightprec"]
                    .doc("Use right preconditioner? Otherwise use left.")
                    .withDefault(false)),
      use_jacobian((*options)["use_jacobian"].withDefault(false)),
      precon_method(
          (*options)["cvode_precon_method"]
              .doc("Preconditioner to use with CVODE Newton iteration. "
                   "Choices: none (default), auto, user, petsc, bbd. "
                   "auto prefers user (if supplied), then petsc (if available), "
                   "then bbd.")
              .withDefault(CvodePreconMethod::none)),
      nvector_type((*options)["nvector"]
                       .doc("N_Vector backend to use: sundials or manyvector")
                       .withDefault(NVectorType::Sundials)),
      cvode_nonlinear_convergence_coef(
          (*options)["cvode_nonlinear_convergence_coef"]
              .doc("Safety factor used in the nonlinear convergence test")
              .withDefault(0.1)),
      cvode_linear_convergence_coef(
          (*options)["cvode_linear_convergence_coef"]
              .doc("Factor by which the Krylov linear solver’s convergence test constant "
                   "is reduced from the nonlinear solver test constant.")
              .withDefault(0.05)),
      save_jacobian((*options)["save_jacobian"]
                        .doc("Save PETSc Jacobian diagnostics to datadir")
                        .withDefault(false)),
      jacobian_export_kind((*options)["jacobian_export_kind"]
                               .doc("Which Jacobian to save for CVODE: system or rhs "
                                    "(scaled is not supported)")
                               .withDefault(bout::JacobianExportKind::system)),
      jacobian_export_trigger(
          (*options)["jacobian_export_trigger"]
              .doc("When to save CVODE Jacobians: output or linear_setup")
              .withDefault(bout::CvodeJacobianExportTrigger::linear_setup)),
      failure_counters(
          (*options)["cvode_failure_counters"]
              .doc("Count the causes of CVODE step failures (Newton failures caused by "
                   "Krylov solves, RHS failures or other causes; constraint "
                   "rejections) and print them after each output. Requires SUNDIALS "
                   ">= 7 and Newton iteration. Does not change the solution.")
              .withDefault(false)),
      suncontext(createSUNContext(BoutComm::get())) {
  supports_constraints = false; // This solver doesn't have constraints
  canReset = true;

  if ((*options)["use_precon"].isSet()) {
    throw BoutException("solver:use_precon is deprecated for CVODE and is now "
                        "ignored. To reproduce behaviour of use_precon=true, "
                        "set cvode_precon_method=user. This enables the Hermes-3 "
                        "physical preconditioners.\n");
  }

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
  if (failure_counters) {
    add_int_diagnostic(fc_krylov_solves, "cvode_krylov_solves", "No. of Krylov solves");
    add_int_diagnostic(fc_krylov_accepted, "cvode_krylov_accepted",
                       "No. of Krylov solves that only reduced the residual in the "
                       "first Newton iteration (accepted by CVODE)");
    add_int_diagnostic(fc_krylov_fail_newton, "cvode_krylov_fail_newton",
                       "No. of Newton failures caused by a Krylov solve that did not "
                       "converge or failed recoverably");
    add_int_diagnostic(fc_rhs_fail_newton, "cvode_rhs_fail_newton",
                       "No. of Newton failures caused by a recoverable RHS failure");
    add_int_diagnostic(fc_newton_notconv, "cvode_newton_notconv",
                       "No. of other Newton failures (divergence, iteration limit or "
                       "failed preconditioner setup)");
    add_int_diagnostic(fc_step_solve_fails, "cvode_step_solve_fails",
                       "No. of steps rejected by a failed Newton solve or by the "
                       "inequality constraints");
    add_int_diagnostic(fc_constraint_fails, "cvode_constraint_fails",
                       "No. of steps rejected by the inequality constraints");
    add_int_diagnostic(fc_lin_conv_fails, "cvode_lin_conv_fails",
                       "No. of Krylov solves that did not converge, including the "
                       "accepted ones");
    add_int_diagnostic(fc_jtimes, "cvode_jtimes", "No. of Jacobian-vector products");
    add_BoutReal_diagnostic(fc_pre_wtime, "cvode_precon_wtime",
                            "Wall time in the preconditioner during the last output "
                            "interval (this processor) [s]");
  }
}

CvodeSolver::~CvodeSolver() {
  if (cvode_initialised) {
    N_VDestroy(uvec);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(sun_solver);
    SUNNonlinSolFree(nonlinear_solver);
#if BOUT_HAS_PETSC
    if (petsc_ksp != nullptr) {
      KSPDestroy(&petsc_ksp);
    }
    if (petsc_r != nullptr) {
      VecDestroy(&petsc_r);
    }
    if (petsc_z != nullptr) {
      VecDestroy(&petsc_z);
    }
    if (petsc_x != nullptr) {
      VecDestroy(&petsc_x);
    }
    if (petsc_f != nullptr) {
      VecDestroy(&petsc_f);
    }
#endif
  }
}

#if BOUT_HAS_PETSC
void CvodeSolver::exportMatrixAndMetadata(bout::JacobianExportKind kind,

                                          Mat jacobian) {

  Solver::writeOnceJacobianMetadata("cvode");
  Solver::writeJacobianMatrix(kind, jacobian);
}

void CvodeSolver::saveDiagnosticJacobian(bout::JacobianExportKind kind, Vec x, BoutReal t,
                                         BoutReal gamma) {
  if (kind == bout::JacobianExportKind::scaled) {
    throw BoutException("solver:jacobian_export_kind=scaled is not supported for CVODE");
  }

  petsc_t = t;
  petsc_gamma = gamma;

  PetscPreconditioner diagnostic_preconditioner;
  Field3D index = globalIndex(0);
  BOUT_DO_PETSC(diagnostic_preconditioner.createJacobianPattern(
      index, *options, getLocalN(), n2Dvars(), n3Dvars(), BoutComm::get()));

  if (kind == bout::JacobianExportKind::rhs) {
    BOUT_DO_PETSC(diagnostic_preconditioner.updateColoring(
        CvodeSolver::petscFormRhsFunction, this));
  } else {
    BOUT_DO_PETSC(
        diagnostic_preconditioner.updateColoring(CvodeSolver::petscFormFunction, this));
  }

  Vec diagnostic_f{nullptr};
  BOUT_DO_PETSC(VecDuplicate(x, &diagnostic_f));
  if (kind == bout::JacobianExportKind::rhs) {
    BOUT_DO_PETSC(CvodeSolver::petscFormRhsFunction(nullptr, x, diagnostic_f, this));
  } else {
    BOUT_DO_PETSC(CvodeSolver::petscFormFunction(nullptr, x, diagnostic_f, this));
  }
  BOUT_DO_PETSC(MatFDColoringSetF(diagnostic_preconditioner.coloring(), diagnostic_f));

  Mat diagnostic_jacobian = diagnostic_preconditioner.jacobian();
  BOUT_DO_PETSC(MatZeroEntries(diagnostic_jacobian));
  BOUT_DO_PETSC(MatFDColoringApply(diagnostic_jacobian,
                                   diagnostic_preconditioner.coloring(), x, nullptr));
  BOUT_DO_PETSC(MatAssemblyBegin(diagnostic_jacobian, MAT_FINAL_ASSEMBLY));
  BOUT_DO_PETSC(MatAssemblyEnd(diagnostic_jacobian, MAT_FINAL_ASSEMBLY));
  exportMatrixAndMetadata(kind, diagnostic_jacobian);
  BOUT_DO_PETSC(VecDestroy(&diagnostic_f));
}

void CvodeSolver::maybeExportJacobian(Mat system_jacobian, Vec x, BoutReal t,
                                      BoutReal gamma) {
  if (!save_jacobian
      or jacobian_export_trigger != bout::CvodeJacobianExportTrigger::linear_setup) {
    return;
  }

  if (jacobian_export_kind == bout::JacobianExportKind::system) {
    exportMatrixAndMetadata(jacobian_export_kind, system_jacobian);
    return;
  }

  saveDiagnosticJacobian(jacobian_export_kind, x, t, gamma);
}

void CvodeSolver::maybeExportOutputJacobian(BoutReal t) {
  if (!save_jacobian
      or jacobian_export_trigger != bout::CvodeJacobianExportTrigger::output) {
    return;
  }

  Vec x{nullptr};
  BOUT_DO_PETSC(VecCreate(BoutComm::get(), &x));
  BOUT_DO_PETSC(VecSetSizes(x, getLocalN(), PETSC_DETERMINE));
  BOUT_DO_PETSC(VecSetType(x, VECMPI));
  BOUT_DO_PETSC(VecSetUp(x));

  BoutReal* xdata = nullptr;
  BOUT_DO_PETSC(VecGetArray(x, &xdata));
  save_vars(xdata);
  BOUT_DO_PETSC(VecRestoreArray(x, &xdata));

  saveDiagnosticJacobian(bout::JacobianExportKind::rhs, x, t, 0.0);
  BOUT_DO_PETSC(VecDestroy(&x));
}
#endif

/**************************************************************************
 * Initialise
 **************************************************************************/

int CvodeSolver::init() {
  const auto backend = nvector_backend();

  Solver::init();
  backend.ensure_manyvector_available();

  output_progress.write("Initialising SUNDIALS' CVODE solver\n");

  // Calculate number of variables (in generic_solver)
  const int local_N = getLocalN();

  // Get total problem size
  int neq{0};
  if (bout::globals::mpi->MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())
      != 0) {
    throw BoutException("Allreduce localN -> GlobalN failed!\n");
  }

  output_info.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n",
                    n3Dvars(), n2Dvars(), neq, local_N);
  output_info.write("\tUsing {} N_Vector backend\n", backend.backend_name());

  // Allocate memory
  uvec = backend.create_state_vector(local_N, neq);

  if (adams_moulton) {
    // By default use functional iteration for Adams-Moulton
    output_info.write("\tUsing Adams-Moulton implicit multistep method\n");
  } else {
    // Use Newton iteration for BDF
    output_info.write("\tUsing BDF method\n");
  }

  const auto lmm = adams_moulton ? CV_ADAMS : CV_BDF;

  cvode_mem = callWithSUNContext(CVodeCreate, suncontext, lmm);
  if (cvode_mem == nullptr) {
    throw BoutException("CVodeCreate failed\n");
  }

  // For callbacks, need pointer to solver object
  if (CVodeSetUserData(cvode_mem, this) != CV_SUCCESS) {
    throw BoutException("CVodeSetUserData failed\n");
  }

#if SUNDIALS_VERSION_MAJOR >= 6
  // Set the default RHS to linear, then pass nonlinear rhs to NL solver
  if (CVodeInit(cvode_mem, cvode_linear_rhs, simtime, uvec) != CV_SUCCESS) {
    throw BoutException("CVodeInit failed\n");
  }
#else
  if (CVodeInit(cvode_mem, cvode_nonlinear_rhs, simtime, uvec) != CV_SUCCESS) {
    throw BoutException("CVodeInit failed\n");
  }
#endif

  if (failure_counters) {
#if SUNDIALS_VERSION_MAJOR < 7
    throw BoutException("solver:cvode_failure_counters requires SUNDIALS >= 7");
#else
    if (func_iter) {
      throw BoutException("solver:cvode_failure_counters requires Newton iteration "
                          "(solver:func_iter = false)");
    }
    // Replace the default Newton solver of CVODE by one created in the same way as
    // in CVodeInit, so that the Krylov wrapper can read the Newton iteration index
    // and failed solves are counted. Attached before CVodeSetMaxNonlinIters,
    // because CVodeSetNonlinearSolver resets the iteration limit.
    nonlinear_solver = callWithSUNContext(SUNNonlinSol_Newton, suncontext, uvec);
    if (nonlinear_solver == nullptr) {
      throw BoutException("SUNNonlinSol_Newton failed\n");
    }
    nonlinear_solver->ops->solve = countingNewtonSolve;
    if (CVodeSetNonlinearSolver(cvode_mem, nonlinear_solver) != CV_SUCCESS) {
      throw BoutException("CVodeSetNonlinearSolver failed\n");
    }
    output_info.write("\tCounting CVODE failures (solver:cvode_failure_counters)\n");
#endif
  }

  if (mxorder > 0) {
    output_warn << "WARNING: Option 'mxorder' is deprecated. Please use "
                   "'cvode_max_order' instead\n";
    if (CVodeSetMaxOrd(cvode_mem, mxorder) != CV_SUCCESS) {
      throw BoutException("CVodeSetMaxOrder failed\n");
    }
  }
  if (max_order > 0) {
    if (CVodeSetMaxOrd(cvode_mem, max_order) != CV_SUCCESS) {
      throw BoutException("CVodeSetMaxOrder failed\n");
    }
  }

  if (CVodeSetStabLimDet(cvode_mem, static_cast<int>(stablimdet)) != CV_SUCCESS) {
    throw BoutException("CVodeSetStabLimDet failed\n");
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

    N_Vector abstolvec =
        backend.clone_vector_like(uvec, "SUNDIALS memory allocation (abstol vector)");

    backend.fill_vector_values(abstolvec, f2dtols, f3dtols);

    if (CVodeSVtolerances(cvode_mem, reltol, abstolvec) != CV_SUCCESS) {
      throw BoutException("CVodeSVtolerances failed\n");
    }

    N_VDestroy(abstolvec);
  } else {
    if (CVodeSStolerances(cvode_mem, reltol, abstol) != CV_SUCCESS) {
      throw BoutException("CVodeSStolerances failed\n");
    }
  }

  if (CVodeSetMaxNumSteps(cvode_mem, mxsteps) != CV_SUCCESS) {
    throw BoutException("CVodeSetMaxNumSteps failed\n");
  }

  if (CVodeSetMaxStep(cvode_mem, max_timestep) != CV_SUCCESS) {
    throw BoutException("CVodeSetMaxStep failed\n");
  }

  if (CVodeSetMinStep(cvode_mem, min_timestep) != CV_SUCCESS) {
    throw BoutException("CVodeSetMinStep failed\n");
  }

  if (CVodeSetInitStep(cvode_mem, start_timestep) != CV_SUCCESS) {
    throw BoutException("CVodeSetInitStep failed\n");
  }

  if (CVodeSetMaxNonlinIters(cvode_mem, max_nonlinear_iterations) != CV_SUCCESS) {
    throw BoutException("CVodeSetMaxNonlinIters failed\n");
  }

#if SUNDIALS_VERSION_MAJOR >= 6
  if (CVodeSetLSetupFrequency(cvode_mem, lsetup_frequency) != CV_SUCCESS) {
    throw BoutException("CVodeSetLSetupFrequency failed\n");
  }
#endif

  if (apply_positivity_constraints) {
    auto f2d_constraints = create_constraints(f2d);
    auto f3d_constraints = create_constraints(f3d);

    N_Vector constraints_vec = backend.clone_vector_like(
        uvec, "SUNDIALS memory allocation (positivity constraints vector)");

    backend.fill_vector_values(constraints_vec, f2d_constraints, f3d_constraints);

    if (CVodeSetConstraints(cvode_mem, constraints_vec) != CV_SUCCESS) {
      throw BoutException("CVodeSetConstraints failed\n");
    }

    N_VDestroy(constraints_vec);
  }

  CvodePreconMethod selected_precon = precon_method;
  if (selected_precon == CvodePreconMethod::Auto) {
    if (hasPreconditioner()) {
      selected_precon = CvodePreconMethod::user;
    } else if (bout::build::has_petsc) {
      selected_precon = CvodePreconMethod::petsc;
    } else {
      selected_precon = CvodePreconMethod::bbd;
    }
  }

  if (save_jacobian) {
#if !BOUT_HAS_PETSC
    throw BoutException("solver:save_jacobian for CVODE requires PETSc support.");
#else
    if (jacobian_export_kind == bout::JacobianExportKind::scaled) {
      throw BoutException(
          "solver:jacobian_export_kind=scaled is not supported for CVODE because "
          "CVODE does not currently apply solver-coordinate scaling.");
    }
    if (jacobian_export_kind == bout::JacobianExportKind::system
        and jacobian_export_trigger != bout::CvodeJacobianExportTrigger::linear_setup) {
      throw BoutException("solver:jacobian_export_kind=system for CVODE requires "
                          "solver:jacobian_export_trigger=linear_setup.");
    }
    if (jacobian_export_trigger == bout::CvodeJacobianExportTrigger::output
        and jacobian_export_kind != bout::JacobianExportKind::rhs) {
      throw BoutException("solver:jacobian_export_trigger=output for CVODE currently "
                          "supports only solver:jacobian_export_kind=rhs.");
    }
    if (jacobian_export_trigger == bout::CvodeJacobianExportTrigger::linear_setup) {
      if (func_iter) {
        throw BoutException("solver:jacobian_export_trigger=linear_setup for CVODE "
                            "requires Newton iteration (set solver:func_iter=false).");
      }
      if (nvector_type == NVectorType::ManyVector) {
        throw BoutException("solver:jacobian_export_trigger=linear_setup for CVODE is "
                            "not supported with solver:nvector=manyvector.");
      }
      if (selected_precon != CvodePreconMethod::petsc) {
        throw BoutException("solver:jacobian_export_trigger=linear_setup for CVODE "
                            "currently requires solver:cvode_precon_method=petsc "
                            "(or auto resolving to petsc).");
      }
    }
    if (!petsc_lib) {
      petsc_lib = std::make_unique<PetscLib>();
    }
#endif
  }

  /// Newton method can include Preconditioners and Jacobian function
  if (func_iter) {
    output_info.write("\tUsing Functional iteration\n");
    nonlinear_solver = callWithSUNContext(SUNNonlinSol_FixedPoint, suncontext, uvec, 0);
    if (nonlinear_solver == nullptr) {
      throw BoutException("SUNNonlinSol_FixedPoint failed\n");
    }

    if (CVodeSetNonlinearSolver(cvode_mem, nonlinear_solver) != 0) {
      throw BoutException("CVodeSetNonlinearSolver failed\n");
    }
  } else {
    output_info.write("\tUsing Newton iteration\n");

    auto prectype = SUN_PREC_NONE;
    if (selected_precon != CvodePreconMethod::none) {
      if (rightprec) {
        prectype = SUN_PREC_RIGHT;
      } else {
        prectype = SUN_PREC_LEFT;
      }
    }

    switch ((*options)["linear_solver"]
                .doc("Set linear solver type. Default is gmres.")
                .withDefault(bout::linear_solver::gmres)) {
    case bout::linear_solver::gmres:
      sun_solver = callWithSUNContext(SUNLinSol_SPGMR, suncontext, uvec, prectype, maxl);
      break;
    case bout::linear_solver::fgmres:
      sun_solver = callWithSUNContext(SUNLinSol_SPFGMR, suncontext, uvec, prectype, maxl);
      break;
    case bout::linear_solver::tfqmr:
      sun_solver =
          callWithSUNContext(SUNLinSol_SPTFQMR, suncontext, uvec, prectype, maxl);
      break;
    case bout::linear_solver::bcgs:
      sun_solver = callWithSUNContext(SUNLinSol_SPBCGS, suncontext, uvec, prectype, maxl);
      break;
    };
    if (sun_solver == nullptr) {
      throw BoutException("Creating SUNDIALS linear solver failed\n");
    }
#if SUNDIALS_VERSION_MAJOR >= 7
    if (failure_counters) {
      sun_solver = makeCountingLinearSolver(sun_solver, this, suncontext);
    }
#endif
    if (CVodeSetLinearSolver(cvode_mem, sun_solver, nullptr) != CVLS_SUCCESS) {
      throw BoutException("CVodeSetLinearSolver failed\n");
    }

#if SUNDIALS_VERSION_MAJOR >= 6
    if (CVodeSetJacEvalFrequency(cvode_mem, jac_eval_frequency) != CVLS_SUCCESS) {
      throw BoutException("CVodeSetJacEvalFrequency failed\n");
    }
#endif

    if (selected_precon == CvodePreconMethod::none) {
      output_info.write("\tNo preconditioning\n");

    } else if (selected_precon == CvodePreconMethod::user) {
      if (!hasPreconditioner()) {
        throw BoutException(
            "solver:cvode_precon_method=user requested, but no user preconditioner "
            "has been supplied.");
      }

      output_info.write("\tUsing user-supplied preconditioner\n");

      if (CVodeSetPreconditioner(cvode_mem, nullptr, cvode_pre) != CVLS_SUCCESS) {
        throw BoutException("CVodeSetPreconditioner failed\n");
      }
    } else if (selected_precon == CvodePreconMethod::petsc) {
#if !BOUT_HAS_PETSC
      throw BoutException(
          "solver:cvode_precon_method=petsc requested, but BOUT++ was not "
          "configured with PETSc.");
#else
      output_info.write("\tUsing PETSc coloring preconditioner\n");

      if (!petsc_lib) {
        petsc_lib = std::make_unique<PetscLib>();
      }

      petsc_global_N = neq;
      petsc_rhs_tmp.resize(local_N);

      if (petsc_ksp == nullptr) {
        PetscCall(KSPCreate(BoutComm::get(), &petsc_ksp));
        PetscCall(KSPSetType(petsc_ksp, KSPPREONLY));
        PetscCall(KSPSetOptionsPrefix(petsc_ksp, "cvode_petscpre_"));
        PetscCall(KSPSetFromOptions(petsc_ksp));
      }

      if (petsc_x == nullptr) {
        PetscCall(VecCreateMPI(BoutComm::get(), local_N, petsc_global_N, &petsc_x));
        PetscCall(VecDuplicate(petsc_x, &petsc_f));
        PetscCall(VecCreateMPI(BoutComm::get(), local_N, petsc_global_N, &petsc_r));
        PetscCall(VecCreateMPI(BoutComm::get(), local_N, petsc_global_N, &petsc_z));
      }

      Field3D index = globalIndex(0);
      PetscCall(petsc_preconditioner.createJacobianPattern(
          index, *options, local_N, n2Dvars(), n3Dvars(), BoutComm::get()));
      PetscCall(
          petsc_preconditioner.updateColoring(CvodeSolver::petscFormFunction, this));
      PetscCall(MatFDColoringSetF(petsc_preconditioner.coloring(), petsc_f));

      if (CVodeSetPreconditioner(cvode_mem, petscPSetup, petscPSolve) != CVLS_SUCCESS) {
        throw BoutException("CVodeSetPreconditioner (PETSc coloring) failed\n");
      }
#endif // BOUT_HAS_PETSC

    } else if (selected_precon == CvodePreconMethod::bbd) {
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
            const Mesh* localmesh = fvar.var->getMesh();
            return a + localmesh->xend - localmesh->xstart + 3;
          });

      const auto mudq = (*options)["mudq"].withDefault(band_width_default);
      const auto mldq = (*options)["mldq"].withDefault(band_width_default);
      const auto mukeep = (*options)["mukeep"].withDefault(n3Dvars() + n2Dvars());
      const auto mlkeep = (*options)["mlkeep"].withDefault(n3Dvars() + n2Dvars());

      if (CVBBDPrecInit(cvode_mem, local_N, mudq, mldq, mukeep, mlkeep, 0.0,
                        cvode_bbd_rhs, nullptr)
          != CVLS_SUCCESS) {
        throw BoutException("CVBBDPrecInit failed\n");
      }
    }

    /// Set Jacobian-vector multiplication function
    if (use_jacobian and hasJacobian()) {
      output_info.write("\tUsing user-supplied Jacobian function\n");

      if (CVodeSetJacTimes(cvode_mem, nullptr, cvode_jac) != CVLS_SUCCESS) {
        throw BoutException("CVodeSetJacTimes failed\n");
      }
    } else {
      output_info.write("\tUsing difference quotient approximation for Jacobian\n");
    }
  }

#if SUNDIALS_VERSION_MAJOR >= 6
  // Set the RHS function to be used in the nonlinear solver
  CVodeSetNlsRhsFn(cvode_mem, cvode_nonlinear_rhs);
#endif

  // Set internal tolerance factors
  if (CVodeSetNonlinConvCoef(cvode_mem, cvode_nonlinear_convergence_coef) != CV_SUCCESS) {
    throw BoutException("CVodeSetNonlinConvCoef failed\n");
  }

  if (CVodeSetEpsLin(cvode_mem, cvode_linear_convergence_coef) != CV_SUCCESS) {
    throw BoutException("CVodeSetEpsLin failed\n");
  }

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
                           .withDefault(bout::positivity_constraint::none);
                   switch (value) {
                   case bout::positivity_constraint::none:
                     return 0.0;
                   case bout::positivity_constraint::positive:
                     return 2.0;
                   case bout::positivity_constraint::non_negative:
                     return 1.0;
                   case bout::positivity_constraint::negative:
                     return -2.0;
                   case bout::positivity_constraint::non_positive:
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

  if (!cvode_initialised) {
    throw BoutException("CvodeSolver not initialised\n");
  }

  for (int i = 1; i <= getNumberOutputSteps(); i++) {

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
    CVodeGetNumPrecSolves(cvode_mem, &temp_long_int);
    npevals = int(temp_long_int);
    CVodeGetNumLinIters(cvode_mem, &temp_long_int);
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

    if (failure_counters) {
      updateFailureCounters();
    }

    /// Call the monitor function

    if (call_monitors(simtime, i, getNumberOutputSteps()) != 0) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

BoutReal CvodeSolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run({})", tout);
  const auto backend = nvector_backend();

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
      const BoutReal last_time = internal_time;
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
  backend.copy_state_from_vector(uvec);

  if (flag < 0) {
    throw BoutException("ERROR CVODE solve failed at t = {:e}, flag = {:d}\n", simtime,
                        flag);
  }

#if BOUT_HAS_PETSC
  maybeExportOutputJacobian(simtime);
#endif

  // Call rhs function to get extra variables at this time
  run_rhs(simtime);

  return simtime;
}

/**************************************************************************
 * RHS function du = F(t, u)
 **************************************************************************/

void CvodeSolver::rhs(BoutReal t, N_Vector u, N_Vector du, bool linear) {
  TRACE("Running RHS: CvodeSolver::res({})", t);
  const auto backend = nvector_backend();

  // Load state from udata
  backend.copy_state_from_vector(u);

  // Get the current timestep
  // Note: CVodeGetCurrentStep updated too late in older versions
  CVodeGetLastStep(cvode_mem, &hcur);

  // Call RHS function
  run_rhs(t, linear);

  // Save derivatives to dudata
  backend.copy_state_to_vector(u);
  backend.copy_deriv_to_vector(du);
}

void CvodeSolver::rhs(BoutReal t, BoutReal* u, BoutReal* du, bool linear) {
  TRACE("Running RHS: CvodeSolver::res({})", t);

  // Load state
  load_vars(u);

  // Get the current timestep
  // Note: CVodeGetCurrentStep updated too late in older versions
  CVodeGetLastStep(cvode_mem, &hcur);

  // Call RHS function
  run_rhs(t, linear);

  // Save derivatives to du
  save_derivs(du);
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void CvodeSolver::pre(BoutReal t, BoutReal gamma, BoutReal delta, N_Vector u,
                      N_Vector rvec, N_Vector zvec) {
  TRACE("Running preconditioner: CvodeSolver::pre({})", t);
  const auto backend = nvector_backend();

  const BoutReal tstart = bout::globals::mpi->MPI_Wtime();

  if (!hasPreconditioner()) {
    // Identity (but should never happen)
    N_VScale(1.0, rvec, zvec);
    return;
  }

  // Load state from udata (as with res function)
  backend.copy_state_from_vector(u);

  // Load vector to be inverted into F_vars
  backend.copy_deriv_from_vector(rvec);

  runPreconditioner(t, gamma, delta);

  // Save the solution from F_vars
  backend.copy_state_to_vector(u);
  backend.copy_deriv_to_vector(zvec);

  pre_Wtime += bout::globals::mpi->MPI_Wtime() - tstart;
  pre_ncalls++;
}

/**************************************************************************
 * Jacobian-vector multiplication function
 **************************************************************************/

void CvodeSolver::jac(BoutReal t, N_Vector y, N_Vector v, N_Vector Jv) {
  TRACE("Running Jacobian: CvodeSolver::jac({})", t);
  const auto backend = nvector_backend();

  if (not hasJacobian()) {
    throw BoutException("No jacobian function supplied!\n");
  }

  // Load state from ydate
  backend.copy_state_from_vector(y);

  // Load vector to be multiplied into F_vars
  backend.copy_deriv_from_vector(v);

  // Call function
  runJacobian(t);

  // Save Jv from vars
  backend.copy_state_to_vector(y);
  backend.copy_deriv_to_vector(Jv);
}

/**************************************************************************
 * CVODE RHS functions
 **************************************************************************/

// NOLINTBEGIN(readability-identifier-length)
namespace {
int cvode_linear_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data) {
  auto* s = static_cast<CvodeSolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs(t, u, du, true);
  } catch (BoutRhsFail& error) {
    s->countRhsFailure(true);
    return 1;
  }
  return 0;
}

int cvode_nonlinear_rhs(BoutReal t, N_Vector u, N_Vector du, void* user_data) {
  auto* s = static_cast<CvodeSolver*>(user_data);

  // Calculate RHS function
  try {
    s->rhs(t, u, du, false);
  } catch (BoutRhsFail& error) {
    s->countRhsFailure(false);
    return 1;
  }
  return 0;
}

/// RHS function for BBD preconditioner
int cvode_bbd_rhs(sunindextype UNUSED(Nlocal), BoutReal t, N_Vector u, N_Vector du,
                  void* user_data) {
  return cvode_linear_rhs(t, u, du, user_data);
}

/// Preconditioner function
int cvode_pre(BoutReal t, N_Vector yy, N_Vector UNUSED(yp), N_Vector rvec, N_Vector zvec,
              BoutReal gamma, BoutReal delta, int UNUSED(lr), void* user_data) {
  auto* s = static_cast<CvodeSolver*>(user_data);

  // Calculate residuals
  s->pre(t, gamma, delta, yy, rvec, zvec);

  return 0;
}

/// Jacobian-vector multiplication function
int cvode_jac(N_Vector v, N_Vector Jv, BoutReal t, N_Vector y, N_Vector UNUSED(fy),
              void* user_data, N_Vector UNUSED(tmp)) {
  auto* s = static_cast<CvodeSolver*>(user_data);

  s->jac(t, y, v, Jv);

  return 0;
}
} // namespace
// NOLINTEND(readability-identifier-length)

#if BOUT_HAS_PETSC
PetscErrorCode CvodeSolver::petscFormFunction(void* UNUSED(dummy), Vec x, Vec f,
                                              void* ctx) {
  auto* s = static_cast<CvodeSolver*>(ctx);

  PetscInt length = 0;
  PetscCall(VecGetLocalSize(x, &length));
  s->petsc_rhs_tmp.resize(static_cast<std::size_t>(length));

  const BoutReal* xdata = nullptr;
  PetscCall(VecGetArrayRead(x, &xdata));

  BoutReal* fdata = nullptr;
  PetscCall(VecGetArray(f, &fdata));

  try {
    s->rhs(s->petsc_t, const_cast<BoutReal*>(xdata), s->petsc_rhs_tmp.data(), true);
  } catch (BoutRhsFail&) {
    PetscCall(VecRestoreArrayRead(x, &xdata));
    PetscCall(VecRestoreArray(f, &fdata));
    return 1;
  }

  for (PetscInt i = 0; i < length; ++i) {
    fdata[i] = xdata[i] - (s->petsc_gamma * s->petsc_rhs_tmp[i]);
  }

  PetscCall(VecRestoreArrayRead(x, &xdata));
  PetscCall(VecRestoreArray(f, &fdata));

  return PETSC_SUCCESS;
}

PetscErrorCode CvodeSolver::petscFormRhsFunction(void* UNUSED(dummy), Vec x, Vec f,
                                                 void* ctx) {
  auto* s = static_cast<CvodeSolver*>(ctx);

  PetscInt length = 0;
  PetscCall(VecGetLocalSize(x, &length));
  s->petsc_rhs_tmp.resize(static_cast<std::size_t>(length));

  const BoutReal* xdata = nullptr;
  PetscCall(VecGetArrayRead(x, &xdata));

  BoutReal* fdata = nullptr;
  PetscCall(VecGetArray(f, &fdata));

  try {
    s->rhs(s->petsc_t, const_cast<BoutReal*>(xdata), s->petsc_rhs_tmp.data(), true);
  } catch (BoutRhsFail&) {
    PetscCall(VecRestoreArrayRead(x, &xdata));
    PetscCall(VecRestoreArray(f, &fdata));
    return 1;
  }

  for (PetscInt i = 0; i < length; ++i) {
    fdata[i] = s->petsc_rhs_tmp[i];
  }

  PetscCall(VecRestoreArrayRead(x, &xdata));
  PetscCall(VecRestoreArray(f, &fdata));

  return PETSC_SUCCESS;
}

int CvodeSolver::petscPSetup(BoutReal t, N_Vector yy, N_Vector UNUSED(yp),
                             CvodeBool UNUSED(jok), CvodeBool* jcurPtr, BoutReal gamma,
                             void* user_data) {
  auto* s = static_cast<CvodeSolver*>(user_data);

  s->petsc_t = t;
  s->petsc_gamma = gamma;

  if (s->nvector_type == NVectorType::ManyVector) {
    throw BoutException("solver:cvode_precon_method=petsc is not supported with "
                        "solver:nvector=manyvector");
  }

  auto* ydata = N_VGetArrayPointer(yy);

  PetscErrorCode ierr = VecPlaceArray(s->petsc_x, ydata);
  if (ierr != 0) {
    return 1;
  }

  ierr = petscFormFunction(nullptr, s->petsc_x, s->petsc_f, s);
  if (ierr != 0) {
    VecResetArray(s->petsc_x);
    return 1;
  }

  ierr = MatFDColoringSetF(s->petsc_preconditioner.coloring(), s->petsc_f);
  if (ierr != 0) {
    VecResetArray(s->petsc_x);
    return 1;
  }

  Mat J = s->petsc_preconditioner.jacobian();
  ierr = MatZeroEntries(J);
  if (ierr != 0) {
    VecResetArray(s->petsc_x);
    return 1;
  }

  ierr = MatFDColoringApply(J, s->petsc_preconditioner.coloring(), s->petsc_x, nullptr);
  if (ierr != 0) {
    VecResetArray(s->petsc_x);
    return 1;
  }

  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  if (ierr == 0) {
    ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  }
  if (ierr != 0) {
    VecResetArray(s->petsc_x);
    return 1;
  }

  s->maybeExportJacobian(J, s->petsc_x, t, gamma);

  ierr = KSPSetOperators(s->petsc_ksp, J, J);
  if (ierr == 0) {
    ierr = KSPSetUp(s->petsc_ksp);
  }

  VecResetArray(s->petsc_x);

  if (ierr != 0) {
    return 1;
  }

  if (jcurPtr != nullptr) {
    *jcurPtr = 1;
  }

  return 0;
}

int CvodeSolver::petscPSolve(BoutReal UNUSED(t), N_Vector UNUSED(yy), N_Vector UNUSED(yp),
                             N_Vector rvec, N_Vector zvec, BoutReal UNUSED(gamma),
                             BoutReal UNUSED(delta), int UNUSED(lr), void* user_data) {
  auto* s = static_cast<CvodeSolver*>(user_data);

  auto* rdata = N_VGetArrayPointer(rvec);
  auto* zdata = N_VGetArrayPointer(zvec);

  PetscErrorCode ierr = VecPlaceArray(s->petsc_r, rdata);
  if (ierr == 0) {
    ierr = VecPlaceArray(s->petsc_z, zdata);
  }
  if (ierr == 0) {
    ierr = KSPSolve(s->petsc_ksp, s->petsc_r, s->petsc_z);
  }

  VecResetArray(s->petsc_r);
  VecResetArray(s->petsc_z);

  return ierr == 0 ? 0 : 1;
}
#endif

/**************************************************************************
 * Failure counters (solver:cvode_failure_counters)
 **************************************************************************/

bout::CvodeKrylovOutcome bout::classifyCvodeKrylovFlag(int flag, int newton_iteration) {
  // Same classification as in cvLsSolve (SUNDIALS cvode_ls.c): a reduced residual
  // is accepted in the first Newton iteration; every other recoverable flag fails
  // the Newton iteration
  switch (flag) {
  case 0: // SUN_SUCCESS
    return CvodeKrylovOutcome::converged;
  case SUNLS_RES_REDUCED:
    return (newton_iteration == 0) ? CvodeKrylovOutcome::accepted
                                   : CvodeKrylovOutcome::newton_failure;
  case SUNLS_CONV_FAIL:
  case SUNLS_ATIMES_FAIL_REC:
  case SUNLS_PSOLVE_FAIL_REC:
  case SUNLS_PACKAGE_FAIL_REC:
  case SUNLS_QRFACT_FAIL:
  case SUNLS_LUFACT_FAIL:
    return CvodeKrylovOutcome::newton_failure;
  default:
    return CvodeKrylovOutcome::unrecoverable;
  }
}

void CvodeSolver::countKrylovFlag(int flag) {
  if (!failure_counters) {
    return;
  }
  ++fc_krylov_solves;
  int curiter = 0;
  if (flag != 0 and nonlinear_solver != nullptr) {
    SUNNonlinSolGetCurIter(nonlinear_solver, &curiter);
  }
  switch (bout::classifyCvodeKrylovFlag(flag, curiter)) {
  case bout::CvodeKrylovOutcome::accepted:
    ++fc_krylov_accepted;
    break;
  case bout::CvodeKrylovOutcome::newton_failure:
    ++fc_krylov_fail_newton;
    if (flag == SUNLS_RES_REDUCED) {
      ++fc_krylov_resred_later;
    } else if (flag == SUNLS_CONV_FAIL) {
      ++fc_krylov_convfail;
    }
    break;
  case bout::CvodeKrylovOutcome::converged:
  case bout::CvodeKrylovOutcome::unrecoverable:
    break;
  }
}

void CvodeSolver::countNewtonSolveFailure() {
  if (failure_counters) {
    ++fc_newton_solve_fails;
  }
}

void CvodeSolver::countRhsFailure(bool linear) {
  if (!failure_counters) {
    return;
  }
  if (linear) {
    ++fc_rhs_fail_linear;
  } else {
    ++fc_rhs_fail_newton;
  }
}

void CvodeSolver::updateFailureCounters() {
  long int temp_long_int = 0;
#if SUNDIALS_VERSION_AT_LEAST(6, 2, 0)
  CVodeGetNumStepSolveFails(cvode_mem, &temp_long_int);
  fc_step_solve_fails = static_cast<int>(temp_long_int);
#endif
  CVodeGetNumLinConvFails(cvode_mem, &temp_long_int);
  fc_lin_conv_fails = static_cast<int>(temp_long_int);
  CVodeGetNumJtimesEvals(cvode_mem, &temp_long_int);
  fc_jtimes = static_cast<int>(temp_long_int);
  CVodeGetNumLinRhsEvals(cvode_mem, &temp_long_int);
  fc_lin_rhs = static_cast<int>(temp_long_int);
  CVodeGetNumPrecEvals(cvode_mem, &temp_long_int);
  fc_prec_setups = static_cast<int>(temp_long_int);

  // A Newton iteration (failures counted by CVODE in nonlin_fails) stops at its
  // first failure: in the residual (RHS), in the linear solve (Krylov), in the
  // convergence test (divergence), at the iteration limit, or in the
  // preconditioner setup. The last three are not visible here and are counted
  // together as the remainder.
  fc_newton_notconv = nonlin_fails - fc_krylov_fail_newton - fc_rhs_fail_newton;
  // A step is rejected by a failed Newton solve or by the constraints. One Newton
  // solve can contain more than one failed Newton iteration (a retry after a new
  // preconditioner setup), so nonlin_fails cannot be used here.
  fc_constraint_fails = fc_step_solve_fails - fc_newton_solve_fails;
  // pre_Wtime is reset at the start of each output interval
  fc_pre_wtime = pre_Wtime;

  // Cumulative values, and the values in this output interval (differences to the
  // previous output; a counter below its previous value was reset by CVodeReInit)
  const std::vector<long> now{nsteps,
                              num_fails,
                              fc_step_solve_fails,
                              fc_newton_solve_fails,
                              fc_constraint_fails,
                              nonlin_fails,
                              fc_krylov_fail_newton,
                              fc_krylov_resred_later,
                              fc_krylov_convfail,
                              fc_rhs_fail_newton,
                              fc_newton_notconv,
                              fc_krylov_solves,
                              fc_lin_conv_fails,
                              fc_krylov_accepted,
                              npevals,
                              fc_jtimes,
                              fc_lin_rhs,
                              fc_rhs_fail_linear,
                              fc_prec_setups};
  if (fc_previous.size() != now.size()) {
    fc_previous.assign(now.size(), 0);
  }
  std::vector<long> interval(now.size());
  for (std::size_t k = 0; k < now.size(); ++k) {
    interval[k] = (now[k] >= fc_previous[k]) ? now[k] - fc_previous[k] : now[k];
  }
  fc_previous = now;
  // "<value in this interval> (<cumulative value>)"
  const auto both = [&](std::size_t k) {
    return fmt::format("{:d} ({:d})", interval[k], now[k]);
  };

  output.write("    -> Failure counters, this output interval (cumulative):\n");
  output.write("       steps {}, err-test-fails {}, step-solve-fails {} = "
               "newton-solve-fails {} + constraint-rejects {}\n",
               both(0), both(1), both(2), both(3), both(4));
  output.write("       newton-failures {} = krylov-not-converged {} [res-reduced {:d}, "
               "conv-fail {:d}, other {:d}] + rhs-recoverable {} + "
               "newton-not-converged {}\n",
               both(5), both(6), interval[7], interval[8],
               interval[6] - interval[7] - interval[8], both(9), both(10));
  output.write("       krylov-solves {}, krylov-failures-all {}, "
               "krylov-accepted-reduced {}\n",
               both(11), both(12), both(13));
  output.write("       psolve {}, psolve-wall {:.3f} s, jtimes {}, rhs-for-jtimes {}, "
               "rhs-linear-recoverable {}, prec-setups {}\n",
               both(14), fc_pre_wtime, both(15), both(16), both(17), both(18));
}

void CvodeSolver::resetInternalFields() {
  const auto backend = nvector_backend();

  backend.copy_state_to_vector(uvec);

  if (CVodeReInit(cvode_mem, simtime, uvec) != CV_SUCCESS) {
    throw BoutException("CVodeReInit failed\n");
  }
  if (failure_counters) {
    // CVodeReInit resets the CVODE counters; reset the own counters as well
    fc_krylov_solves = fc_krylov_accepted = fc_krylov_fail_newton = 0;
    fc_krylov_resred_later = fc_krylov_convfail = 0;
    fc_rhs_fail_newton = fc_rhs_fail_linear = fc_newton_solve_fails = 0;
    fc_previous.clear();
  }
}

#if SUNDIALS_VERSION_MAJOR >= 7
// NOLINTBEGIN(readability-identifier-length)
namespace {
int countingNewtonSolve(SUNNonlinearSolver nls, N_Vector y0, N_Vector ycor, N_Vector w,
                        sunrealtype tol, sunbooleantype call_lsetup, void* mem) {
  const int flag = SUNNonlinSolSolve_Newton(nls, y0, ycor, w, tol, call_lsetup, mem);
  if (flag != SUN_SUCCESS) {
    // CVODE passes its memory block, which holds the solver as user data
    void* user_data = nullptr;
    if (CVodeGetUserData(mem, &user_data) == CV_SUCCESS and user_data != nullptr) {
      static_cast<CvodeSolver*>(user_data)->countNewtonSolveFailure();
    }
  }
  return flag;
}

struct CountingLinearSolverContent {
  SUNLinearSolver inner;
  CvodeSolver* owner;
};

SUNLinearSolver innerSolver(SUNLinearSolver S) {
  return static_cast<CountingLinearSolverContent*>(S->content)->inner;
}

SUNLinearSolver_Type countingGetType(SUNLinearSolver S) {
  return SUNLinSolGetType(innerSolver(S));
}
SUNLinearSolver_ID countingGetID(SUNLinearSolver S) {
  return SUNLinSolGetID(innerSolver(S));
}
SUNErrCode countingSetATimes(SUNLinearSolver S, void* A_data, SUNATimesFn ATimes) {
  return SUNLinSolSetATimes(innerSolver(S), A_data, ATimes);
}
SUNErrCode countingSetPreconditioner(SUNLinearSolver S, void* P_data, SUNPSetupFn Pset,
                                     SUNPSolveFn Psol) {
  return SUNLinSolSetPreconditioner(innerSolver(S), P_data, Pset, Psol);
}
SUNErrCode countingSetScalingVectors(SUNLinearSolver S, N_Vector s1, N_Vector s2) {
  return SUNLinSolSetScalingVectors(innerSolver(S), s1, s2);
}
SUNErrCode countingSetZeroGuess(SUNLinearSolver S, sunbooleantype onoff) {
  return SUNLinSolSetZeroGuess(innerSolver(S), onoff);
}
SUNErrCode countingInitialize(SUNLinearSolver S) {
  return SUNLinSolInitialize(innerSolver(S));
}
int countingSetup(SUNLinearSolver S, SUNMatrix A) {
  return SUNLinSolSetup(innerSolver(S), A);
}
int countingSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b,
                  sunrealtype tol) {
  auto* content = static_cast<CountingLinearSolverContent*>(S->content);
  const int flag = SUNLinSolSolve(content->inner, A, x, b, tol);
  content->owner->countKrylovFlag(flag);
  return flag;
}
int countingNumIters(SUNLinearSolver S) { return SUNLinSolNumIters(innerSolver(S)); }
sunrealtype countingResNorm(SUNLinearSolver S) {
  return SUNLinSolResNorm(innerSolver(S));
}
sunindextype countingLastFlag(SUNLinearSolver S) {
  return SUNLinSolLastFlag(innerSolver(S));
}
SUNErrCode countingSpace(SUNLinearSolver S, long int* lenrwLS, long int* leniwLS) {
  return SUNLinSolSpace(innerSolver(S), lenrwLS, leniwLS);
}
N_Vector countingResid(SUNLinearSolver S) { return SUNLinSolResid(innerSolver(S)); }
SUNErrCode countingFree(SUNLinearSolver S) {
  if (S == nullptr) {
    return SUN_SUCCESS;
  }
  auto* content = static_cast<CountingLinearSolverContent*>(S->content);
  if (content != nullptr) {
    SUNLinSolFree(content->inner);
    delete content;
    S->content = nullptr;
  }
  SUNLinSolFreeEmpty(S);
  return SUN_SUCCESS;
}

SUNLinearSolver makeCountingLinearSolver(SUNLinearSolver inner, CvodeSolver* owner,
                                         SUNContext ctx) {
  SUNLinearSolver S = SUNLinSolNewEmpty(ctx);
  if (S == nullptr) {
    throw BoutException("SUNLinSolNewEmpty failed\n");
  }
  S->content = new CountingLinearSolverContent{inner, owner};
  // Offer exactly the operations of the wrapped solver (all operations of
  // SUNLinearSolver in SUNDIALS 7.2), because CVODE checks which ones exist
  const auto* ops = inner->ops;
  S->ops->gettype = (ops->gettype != nullptr) ? countingGetType : nullptr;
  S->ops->getid = (ops->getid != nullptr) ? countingGetID : nullptr;
  S->ops->setatimes = (ops->setatimes != nullptr) ? countingSetATimes : nullptr;
  S->ops->setpreconditioner =
      (ops->setpreconditioner != nullptr) ? countingSetPreconditioner : nullptr;
  S->ops->setscalingvectors =
      (ops->setscalingvectors != nullptr) ? countingSetScalingVectors : nullptr;
  S->ops->setzeroguess = (ops->setzeroguess != nullptr) ? countingSetZeroGuess : nullptr;
  S->ops->initialize = (ops->initialize != nullptr) ? countingInitialize : nullptr;
  S->ops->setup = (ops->setup != nullptr) ? countingSetup : nullptr;
  S->ops->solve = (ops->solve != nullptr) ? countingSolve : nullptr;
  S->ops->numiters = (ops->numiters != nullptr) ? countingNumIters : nullptr;
  S->ops->resnorm = (ops->resnorm != nullptr) ? countingResNorm : nullptr;
  S->ops->lastflag = (ops->lastflag != nullptr) ? countingLastFlag : nullptr;
  S->ops->space = (ops->space != nullptr) ? countingSpace : nullptr;
  S->ops->resid = (ops->resid != nullptr) ? countingResid : nullptr;
  S->ops->free = countingFree;
  return S;
}
} // namespace
// NOLINTEND(readability-identifier-length)
#endif
#endif
