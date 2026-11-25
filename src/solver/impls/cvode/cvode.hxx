/**************************************************************************
 * Interface to SUNDIALS CVODE
 *
 **************************************************************************
 * Copyright 2010 - 2026 BOUT++ contributors
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

#ifndef BOUT_SUNDIAL_SOLVER_H
#define BOUT_SUNDIAL_SOLVER_H

#include "bout/bout_enum_class.hxx"
#include "bout/build_defines.hxx"
#include "bout/solver.hxx"
#include <memory>

#if not BOUT_HAS_CVODE

namespace {
RegisterUnavailableSolver
    registerunavailablecvode("cvode", "BOUT++ was not configured with CVODE/SUNDIALS");
}

#else

#include "../../sundials_nvector_interface.hxx"
#include "bout/bout_types.hxx"
#include "bout/petsc_preconditioner.hxx"
#include "bout/region.hxx"
#include "bout/sundials_backports.hxx"

#if BOUT_HAS_PETSC
#include "bout/petsclib.hxx"

#include <petscksp.h>
#include <petscvec.h>
#endif

#include <string>
#include <vector>

class CvodeSolver;
class Options;

namespace {
RegisterSolver<CvodeSolver> registersolvercvode("cvode");
}

// Preconditioner selection for CVODE.
// Note: String comparisons are case-insensitive so "Auto" avoids conflict with keyword
BOUT_ENUM_CLASS(CvodePreconMethod, none, Auto, user, petsc, bbd);
BOUT_ENUM_CLASS_NS(bout, CvodeJacobianExportTrigger,
                   output,      ///< Export once per solver output timestep
                   linear_setup ///< Export whenever CVODE rebuilds linear solver data
);

#if SUNDIALS_VERSION_AT_LEAST(6, 0, 0)
using CvodeBool = sunbooleantype;
#else
using CvodeBool = booleantype;
#endif

namespace bout {
/// What CVODE does with the return flag of a Krylov solve (as in cvLsSolve,
/// SUNDIALS cvode_ls.c). Used by solver:cvode_failure_counters
enum class CvodeKrylovOutcome {
  converged,      ///< SUN_SUCCESS
  accepted,       ///< SUNLS_RES_REDUCED in the first Newton iteration
  newton_failure, ///< Other recoverable flags: the Newton iteration fails
  unrecoverable   ///< CVODE stops with an error
};

/// Classify the return flag \p flag of a Krylov solve in Newton iteration
/// \p newton_iteration (0 is the first)
CvodeKrylovOutcome classifyCvodeKrylovFlag(int flag, int newton_iteration);
} // namespace bout

class CvodeSolver : public Solver {
public:
  explicit CvodeSolver(Options* opts = nullptr);
  ~CvodeSolver() override;

  BoutReal getCurrentTimestep() override { return hcur; }

  int init() override;
  int run() override;
  BoutReal run(BoutReal tout);

  void resetInternalFields() override;

  // These functions are used internally (but need to be public)
  void rhs(BoutReal t, N_Vector u, N_Vector du, bool linear);
  /// This version will copy data to/from a linear vector into BOUT++ fields
  void rhs(BoutReal t, BoutReal* u, BoutReal* du, bool linear);
  void pre(BoutReal t, BoutReal gamma, BoutReal delta, N_Vector u, N_Vector rvec,
           N_Vector zvec);
  void jac(BoutReal t, N_Vector y, N_Vector v, N_Vector Jv);

  /// Failure counters (solver:cvode_failure_counters): count the return flag of one
  /// Krylov solve, classified by bout::classifyCvodeKrylovFlag with the iteration
  /// index of the Newton solver
  void countKrylovFlag(int flag);
  /// Failure counters: one call of the Newton solver failed
  void countNewtonSolveFailure();
  /// Failure counters: a recoverable RHS failure (BoutRhsFail) in the linear RHS
  /// (Jacobian-vector products, BBD preconditioner) or in the nonlinear RHS
  /// (Newton residual)
  void countRhsFailure(bool linear);

private:
#if BOUT_HAS_PETSC
  /// PETSc callback for the CVODE system function ``x - gamma * rhs(t, x)``.
  static PetscErrorCode petscFormFunction(void* dummy, Vec x, Vec f, void* ctx);
  /// PETSc callback for the raw CVODE RHS ``rhs(t, x)`` with ``linear=true``.
  static PetscErrorCode petscFormRhsFunction(void* dummy, Vec x, Vec f, void* ctx);
  static int petscPSetup(BoutReal t, N_Vector yy, N_Vector yp, CvodeBool jok,
                         CvodeBool* jcurPtr, BoutReal gamma, void* user_data);
  static int petscPSolve(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rvec,
                         N_Vector zvec, BoutReal gamma, BoutReal delta, int lr,
                         void* user_data);
  /// Write one Jacobian matrix and the shared metadata JSON.
  void exportMatrixAndMetadata(bout::JacobianExportKind kind, Mat jacobian);
  /// Build and save a diagnostic Jacobian of the requested kind.
  void saveDiagnosticJacobian(bout::JacobianExportKind kind, Vec x, BoutReal t,
                              BoutReal gamma);
  /// Save Jacobians during CVODE linear solver setup callbacks.
  void maybeExportJacobian(Mat system_jacobian, Vec x, BoutReal t, BoutReal gamma);
  /// Save Jacobians once per solver output timestep.
  void maybeExportOutputJacobian(BoutReal t);
#endif

  BoutReal hcur; //< Current internal timestep

  bool diagnose{false}; //< Output additional diagnostics

  N_Vector uvec{nullptr};   //< Values
  void* cvode_mem{nullptr}; //< CVODE internal memory block

  BoutReal pre_Wtime{0.0}; //< Time in preconditioner
  int pre_ncalls{0};       //< Number of calls to preconditioner

  /// Use Adams Moulton implicit multistep. Otherwise BDF method
  bool adams_moulton;
  /// Use functional iteration instead of Newton
  bool func_iter;
  /// Maximum order of method to use. < 0 means no limit
  int max_order;
  bool stablimdet;
  /// Absolute tolerance
  BoutReal abstol;
  /// Relative tolerance
  BoutReal reltol;
  /// Use separate absolute tolerance for each field
  bool use_vector_abstol;
  /// Maximum number of internal steps between outputs.
  int mxsteps;
  /// Maximum time step size
  BoutReal max_timestep;
  /// Minimum time step size
  BoutReal min_timestep;
  /// Starting time step. < 0 then chosen by CVODE.
  BoutReal start_timestep;
  /// Maximum order
  int mxorder;
  /// Maximum number of nonlinear iterations allowed by CVODE before
  /// reducing timestep. CVODE default (used if this option is
  /// negative) is 3
  int max_nonlinear_iterations;
  /// Max number of steps between calls to linear solver setup
  long int lsetup_frequency;
  /// Max number of steps between Jacobian/preconditioner evaluations
  long int jac_eval_frequency;
  /// Use CVODE function CVodeSetConstraints to constrain variables -
  /// the constraint to be applied is set by the positivity_constraint
  /// option in the subsection for each variable
  bool apply_positivity_constraints;
  /// Maximum number of linear iterations
  int maxl;
  /// Use right preconditioner? Otherwise use left.
  bool rightprec;
  bool use_jacobian;
  CvodePreconMethod precon_method;
  NVectorType nvector_type;
  BoutReal cvode_nonlinear_convergence_coef;
  BoutReal cvode_linear_convergence_coef;
  bool save_jacobian; ///< Save PETSc Jacobian diagnostics to ``datadir``?
  bout::JacobianExportKind
      jacobian_export_kind; ///< Export ``system`` or ``rhs`` Jacobian
  bout::CvodeJacobianExportTrigger
      jacobian_export_trigger; ///< Export on ``output`` or ``linear_setup``
  /// Count the causes of CVODE step failures (solver:cvode_failure_counters)
  bool failure_counters;

  // Diagnostics from CVODE
  int nsteps{0};
  int nfevals{0};
  int nniters{0};
  int npevals{0};
  int nliters{0};
  BoutReal last_step{0.0};
  int last_order{0};
  int num_fails{0};
  int nonlin_fails{0};
  int stab_lims{0};

  // Failure counters (only with failure_counters), cumulative like the
  // diagnostics above
  int fc_krylov_solves{0};       ///< Krylov solves (SUNLinSolSolve calls)
  int fc_krylov_accepted{0};     ///< SUNLS_RES_REDUCED in Newton iteration 0 (accepted)
  int fc_krylov_fail_newton{0};  ///< Krylov solves that failed a Newton iteration
  int fc_krylov_resred_later{0}; ///< of these: SUNLS_RES_REDUCED in a later iteration
  int fc_krylov_convfail{0};     ///< of these: SUNLS_CONV_FAIL
  int fc_rhs_fail_newton{0};     ///< recoverable failures of the nonlinear RHS
  int fc_rhs_fail_linear{0};     ///< recoverable failures of the linear RHS
  int fc_newton_notconv{0};      ///< other Newton failures (divergence, ...)
  int fc_newton_solve_fails{0};  ///< failed calls of the Newton solver
  int fc_step_solve_fails{0};    ///< CVodeGetNumStepSolveFails
  int fc_constraint_fails{0};    ///< steps rejected by the inequality constraints
  int fc_lin_conv_fails{0};      ///< CVodeGetNumLinConvFails
  int fc_jtimes{0};              ///< CVodeGetNumJtimesEvals
  int fc_lin_rhs{0};             ///< CVodeGetNumLinRhsEvals
  int fc_prec_setups{0};         ///< CVodeGetNumPrecEvals
  BoutReal fc_pre_wtime{0.0};    ///< preconditioner wall time, last output interval
  /// Cumulative values at the previous output, for the values per output interval
  std::vector<long> fc_previous;
  /// Update the failure counters from CVODE and print them (called from run())
  void updateFailureCounters();

  bool cvode_initialised{false};

  template <class FieldType>
  std::vector<BoutReal> create_constraints(const std::vector<VarStr<FieldType>>& fields);

  SundialsNVectorInterface nvector_backend() {
    return SundialsNVectorInterface(*this, suncontext, nvector_type);
  }

  /// SPGMR solver structure
  SUNLinearSolver sun_solver{nullptr};
  /// Solver for functional iterations for Adams-Moulton
  SUNNonlinearSolver nonlinear_solver{nullptr};
  /// Context for SUNDIALS memory allocations
  sundials::Context suncontext;

#if BOUT_HAS_PETSC
  // PETSc-coloring-based preconditioning and Jacobian diagnostics for CVODE
  std::unique_ptr<PetscLib> petsc_lib;
  PetscPreconditioner petsc_preconditioner;
  KSP petsc_ksp{nullptr};
  Vec petsc_r{nullptr};
  Vec petsc_z{nullptr};
  Vec petsc_x{nullptr};
  Vec petsc_f{nullptr};
  PetscInt petsc_global_N{0};
  BoutReal petsc_gamma{0.0};
  BoutReal petsc_t{0.0};
  std::vector<BoutReal> petsc_rhs_tmp;
#endif
};

#endif // BOUT_HAS_CVODE
#endif // BOUT_SUNDIAL_SOLVER_H
