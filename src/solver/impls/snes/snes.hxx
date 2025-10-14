/**************************************************************************
 *
 * Finds the steady-state solution of a set of equations
 * using PETSc for the SNES interface
 *
 **************************************************************************
 * Copyright 2015-2025 BOUT++ contributors
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

#ifndef BOUT_SNES_SOLVER_H
#define BOUT_SNES_SOLVER_H

#include <bout/build_defines.hxx>
#include <bout/solver.hxx>

#if BOUT_HAS_PETSC

class SNESSolver;

#include "mpi.h"

#include <bout/bout_enum_class.hxx>
#include <bout/bout_types.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/petsclib.hxx>

#include <petsc.h>
#include <petscsnes.h>

namespace {
RegisterSolver<SNESSolver> registersolversnes("snes");
RegisterSolver<SNESSolver> registersolverbeuler("beuler");
} // namespace

BOUT_ENUM_CLASS(BoutSnesEquationForm, pseudo_transient, rearranged_backward_euler,
                backward_euler, direct_newton);

BOUT_ENUM_CLASS(BoutPTCStrategy,
                inverse_residual, ///< dt = pseudo_alpha / residual
                history_based,    ///< Grow/shrink dt based on residual decrease/increase
                hybrid); ///< Combine inverse_residual and history_based strategies

/// Uses PETSc's SNES interface to find a steady state solution to a
/// nonlinear ODE by integrating in time with Backward Euler
class SNESSolver : public Solver {
public:
  explicit SNESSolver(Options* opts = nullptr);
  ~SNESSolver() override = default;

  int init() override;
  int run() override;

  /// Nonlinear function. This is called by PETSc SNES object
  /// via a static C-style function. For implicit
  /// time integration this function calculates:
  ///
  ///     f = (x - gamma*G(x)) - rhs
  ///
  /// The form depends on equation_form
  ///
  ///
  /// @param[in] x  The state vector
  /// @param[out] f  The vector for the result f(x)
  /// @param[in] linear  Specifies that the SNES solver is in a linear (KSP) inner loop,
  ///                    so the operator should be linearised if possible
  PetscErrorCode snes_function(Vec x, Vec f, bool linear); ///< Nonlinear function

  /// Preconditioner. Called by PCapply
  /// via a C-style static function.
  ///
  /// @param[in] x  The vector to be operated on
  /// @param[out] f  The result of the operation
  PetscErrorCode precon(Vec x, Vec f);

  /// Scale an approximate Jacobian,
  /// and update the internal RHS scaling factors
  /// This is called by SNESComputeJacobianScaledColor with the
  /// finite difference approximated Jacobian.
  PetscErrorCode scaleJacobian(Mat B);

  /// Save diagnostics to output
  void outputVars(Options& output_options, bool save_repeat = true) override;

private:
  PetscErrorCode FDJinitialise();         ///< Finite Difference Jacobian initialise
  PetscErrorCode FDJpruneJacobian();      ///< Remove small elements from the Jacobian
  PetscErrorCode FDJrestoreFromPruning(); ///< Restore Jacobian to original pattern

  /// Call the physics model RHS function
  ///
  /// @param[in] x       The state vector. Will be scaled if scale_vars=true
  /// @param[out] f      The vector for the result f(x)
  /// @param[in] linear  Specifies that the SNES solver is in a linear (KSP) inner loop
  PetscErrorCode rhs_function(Vec x, Vec f, bool linear);

  BoutReal timestep;     ///< Internal timestep
  BoutReal dt;           ///< Current timestep used in snes_function.
  BoutReal dt_min_reset; ///< If dt falls below this, reset solve
  BoutReal max_timestep; ///< Maximum timestep

  /// Form of the equation to solve
  BoutSnesEquationForm equation_form;

  std::string snes_type;
  BoutReal atol; ///< Absolute tolerance
  BoutReal rtol; ///< Relative tolerance
  BoutReal stol; ///< Convergence tolerance
  int maxf; ///< Maximum number of function evaluations allowed in the solver (default: 10000)

  int maxits;               ///< Maximum nonlinear iterations
  int lower_its, upper_its; ///< Limits on iterations for timestep adjustment

  BoutReal timestep_factor_on_failure;
  BoutReal timestep_factor_on_upper_its;
  BoutReal timestep_factor_on_lower_its;

  // Pseudo-Transient Continuation (PTC) variables
  // These are used if equation_form = pseudo_transient
  BoutPTCStrategy pseudo_strategy;  ///< Strategy to use when setting timesteps
  BoutReal pseudo_alpha;            ///< dt = alpha / residual
  BoutReal pseudo_growth_factor;    ///< Timestep increase 1.1 - 1.2
  BoutReal pseudo_reduction_factor; ///< Timestep decrease 0.5
  BoutReal pseudo_max_ratio;        ///< Maximum timestep ratio between neighboring cells
  Vec dt_vec;                       ///< Each quantity can have its own timestep

  /// Initialize the Pseudo-Transient Continuation method
  PetscErrorCode initPseudoTimestepping();
  /// Update dt_vec based on new solution x
  PetscErrorCode updatePseudoTimestepping(Vec x);
  /// Decide the next pseudo-timestep. Called by updatePseudoTimestepping
  BoutReal updatePseudoTimestep(BoutReal previous_timestep, BoutReal previous_residual,
                                BoutReal current_residual);
  BoutReal updatePseudoTimestep_inverse_residual(BoutReal previous_timestep,
                                                 BoutReal current_residual);
  BoutReal updatePseudoTimestep_history_based(BoutReal previous_timestep,
                                              BoutReal previous_residual,
                                              BoutReal current_residual);
  Field3D pseudo_residual; ///< Diagnostic output
  Field2D pseudo_residual_2d;
  Field3D pseudo_timestep;

  ///< PID controller parameters
  bool pid_controller; ///< Use PID controller?
  int target_its;      ///< Target number of nonlinear iterations for the PID controller.
  ///< Use with caution! Not tested values.
  BoutReal kP; ///< (0.6 - 0.8) Proportional parameter (main response to current step)
  BoutReal kI; ///< (0.2 - 0.4) Integral parameter (smooths history of changes)
  BoutReal kD; ///< (0.1 - 0.3) Derivative (dampens oscillation - optional)

  int nl_its_prev;
  int nl_its_prev2;

  BoutReal pid(BoutReal timestep, int nl_its, BoutReal max_dt); ///< Updates the timestep

  bool diagnose;          ///< Output additional diagnostics
  bool diagnose_failures; ///< Print diagnostics on SNES failures

  int nlocal; ///< Number of variables on local processor
  int neq;    ///< Number of variables in total

  PetscLib lib; ///< Handles initialising, finalising PETSc
  Vec snes_f;   ///< Used by SNES to store function
  Vec snes_x;   ///< Result of SNES
  Vec x0;       ///< Solution at start of current timestep
  Vec delta_x;  ///< Change in solution
  Vec output_x; ///< Solution to output. Used if interpolating.

  bool predictor;       ///< Use linear predictor?
  Vec x1;               ///< Previous solution
  BoutReal time1{-1.0}; ///< Time of previous solution

  SNES snes;                         ///< SNES context
  Mat Jmf;                           ///< Matrix Free Jacobian
  Mat Jfd;                           ///< Finite Difference Jacobian
  MatFDColoring fdcoloring{nullptr}; ///< Matrix coloring context
                                     ///< Jacobian evaluation

  bool use_precon;                ///< Use preconditioner
  std::string ksp_type;           ///< Linear solver type
  bool kspsetinitialguessnonzero; ///< Set initial guess to non-zero
  int maxl;                       ///< Maximum linear iterations
  std::string pc_type;            ///< Preconditioner type
  std::string pc_hypre_type;      ///< Hypre preconditioner type
  std::string line_search_type;   ///< Line search type

  bool matrix_free;          ///< Use matrix free Jacobian
  bool matrix_free_operator; ///< Use matrix free Jacobian in the operator?
  int lag_jacobian;          ///< Re-use Jacobian
  bool use_coloring;         ///< Use matrix coloring

  bool jacobian_recalculated; ///< Flag set when Jacobian is recalculated
  bool prune_jacobian;        ///< Remove small elements in the Jacobian?
  BoutReal prune_abstol;      ///< Prune values with absolute values smaller than this
  BoutReal prune_fraction;    ///< Prune if fraction of small elements is larger than this
  bool jacobian_pruned{false}; ///< Has the Jacobian been pruned?
  Mat Jfd_original;            ///< Used to reset the Jacobian if over-pruned
  void updateColoring();       ///< Updates the coloring using Jfd

  bool scale_rhs;          ///< Scale time derivatives?
  Vec rhs_scaling_factors; ///< Factors to multiply RHS function
  Vec jac_row_inv_norms;   ///< 1 / Norm of the rows of the Jacobian

  bool scale_vars;         ///< Scale individual variables?
  Vec var_scaling_factors; ///< Factors to multiply variables when passing to user
  Vec scaled_x;            ///< The values passed to the user RHS
};

#else

namespace {
RegisterUnavailableSolver registerunavailablesnes("snes",
                                                  "BOUT++ was not configured with PETSc");
RegisterUnavailableSolver
    registerunavailablebeuler("beuler", "BOUT++ was not configured with PETSc");
} // namespace

#endif // BOUT_HAS_PETSC

#endif // BOUT_SNES_SOLVER_H
