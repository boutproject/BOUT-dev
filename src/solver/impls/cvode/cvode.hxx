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
#include "bout/region.hxx"
#include "bout/sundials_backports.hxx"

#if BOUT_HAS_PETSC
#include "bout/petsc_preconditioner.hxx"
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

#if SUNDIALS_VERSION_AT_LEAST(6, 0, 0)
using CvodeBool = sunbooleantype;
#else
using CvodeBool = booleantype;
#endif

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

private:
#if BOUT_HAS_PETSC
  static PetscErrorCode petscFormFunction(void* dummy, Vec x, Vec f, void* ctx);
  static int petscPSetup(BoutReal t, N_Vector yy, N_Vector yp, CvodeBool jok,
                         CvodeBool* jcurPtr, BoutReal gamma, void* user_data);
  static int petscPSolve(BoutReal t, N_Vector yy, N_Vector yp, N_Vector rvec,
                         N_Vector zvec, BoutReal gamma, BoutReal delta, int lr,
                         void* user_data);
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
  // PETSc-coloring-based preconditioning for CVODE
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
