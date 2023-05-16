/**************************************************************************
 *
 * Finds the steady-state solution of a set of equations
 * using PETSc for the SNES interface
 *
 **************************************************************************
 * Copyright 2015, 2021 B.D.Dudson
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

#ifndef __SNES_SOLVER_H__
#define __SNES_SOLVER_H__

#include <bout/build_config.hxx>
#include <bout/solver.hxx>

#if BOUT_HAS_PETSC

class SNESSolver;

#include "mpi.h"

#include <bout/bout_enum_class.hxx>
#include <bout/bout_types.hxx>
#include <bout/petsclib.hxx>

#include <petsc.h>
#include <petscsnes.h>

namespace {
RegisterSolver<SNESSolver> registersolversnes("snes");
RegisterSolver<SNESSolver> registersolverbeuler("beuler");
} // namespace

BOUT_ENUM_CLASS(BoutSnesEquationForm, pseudo_transient, rearranged_backward_euler,
                backward_euler, direct_newton);

/// Uses PETSc's SNES interface to find a steady state solution to a
/// nonlinear ODE by integrating in time with Backward Euler
class SNESSolver : public Solver {
public:
  explicit SNESSolver(Options* opts = nullptr);
  ~SNESSolver() = default;

  int init() override;
  int run() override;

  /// Nonlinear function. This is called by PETSc SNES object
  /// via a static C-style function. For implicit
  /// time integration this function calculates:
  ///
  ///     f = (x - gamma*G(x)) - rhs
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

private:
  BoutReal timestep;     ///< Internal timestep
  BoutReal dt;           ///< Current timestep used in snes_function
  BoutReal dt_min_reset; ///< If dt falls below this, reset solve
  BoutReal max_timestep; ///< Maximum timestep

  std::string snes_type;
  BoutReal atol; ///< Absolute tolerance
  BoutReal rtol; ///< Relative tolerance
  BoutReal stol; ///< Convergence tolerance

  int maxits;               ///< Maximum nonlinear iterations
  int lower_its, upper_its; ///< Limits on iterations for timestep adjustment

  bool diagnose;          ///< Output additional diagnostics
  bool diagnose_failures; ///< Print diagnostics on SNES failures

  int nlocal; ///< Number of variables on local processor
  int neq;    ///< Number of variables in total

  /// Form of the equation to solve
  BoutSnesEquationForm equation_form;

  PetscLib lib; ///< Handles initialising, finalising PETSc
  Vec snes_f;   ///< Used by SNES to store function
  Vec snes_x;   ///< Result of SNES
  Vec x0;       ///< Solution at start of current timestep
  Vec delta_x;  ///< Change in solution

  bool predictor;       ///< Use linear predictor?
  Vec x1;               ///< Previous solution
  BoutReal time1{-1.0}; ///< Time of previous solution

  SNES snes;                ///< SNES context
  Mat Jmf;                  ///< Matrix-free Jacobian
  MatFDColoring fdcoloring; ///< Matrix coloring context, used for finite difference
                            ///< Jacobian evaluation

  bool use_precon;                ///< Use preconditioner
  std::string ksp_type;           ///< Linear solver type
  bool kspsetinitialguessnonzero; ///< Set initial guess to non-zero
  int maxl;                       ///< Maximum linear iterations
  std::string pc_type;            ///< Preconditioner type
  std::string line_search_type;   ///< Line search type
  bool matrix_free;               ///< Use matrix free Jacobian
  int lag_jacobian;               ///< Re-use Jacobian
  bool use_coloring;              ///< Use matrix coloring
};

#else

namespace {
RegisterUnavailableSolver registerunavailablesnes("snes",
                                                  "BOUT++ was not configured with PETSc");
RegisterUnavailableSolver
    registerunavailablebeuler("beuler", "BOUT++ was not configured with PETSc");
} // namespace

#endif // BOUT_HAS_PETSC

#endif // __SNES_SOLVER_H__
