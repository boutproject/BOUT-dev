/**************************************************************************
 * Interface to PETSc solver
 *
 **************************************************************************
 * Copyright 2010 - 2025 BOUT++ contributors
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

#ifndef BOUT_PETSC_SOLVER_H
#define BOUT_PETSC_SOLVER_H

#include "bout/build_defines.hxx"
#include "bout/solver.hxx"

#if not BOUT_HAS_PETSC

namespace {
RegisterUnavailableSolver
    registerunavailablepetsc("petsc", "BOUT++ was not configured with PETSc");
}

#else

class PetscSolver;

#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/petsclib.hxx>
#include <bout/vector2d.hxx>
#include <bout/vector3d.hxx>

#include <petsc.h>
#include <petscts.h>

#include <vector>

namespace {
RegisterSolver<PetscSolver> registersolverpetsc("petsc");
}

class PetscSolver : public Solver {
public:
  PetscSolver(Options* opts = nullptr);
  ~PetscSolver();

  int init() override;
  int run() override;

  // These functions used internally (but need to be public)

  /// Wrapper for the RHS function
  PetscErrorCode rhs(PetscReal t, Vec globalin, Vec globalout, bool linear);
  /// Wrapper for the preconditioner
  PetscErrorCode pre(Vec x, Vec y);

  // Call back functions that need to access internal state
  friend PetscErrorCode PetscMonitor(TS, PetscInt, PetscReal, Vec, void* ctx);

  friend PetscErrorCode solver_ijacobian(TS, BoutReal, Vec, Vec, PetscReal shift, Mat J,
                                         Mat Jpre, void* ctx);

private:
  BoutReal shift; ///< Shift (alpha) parameter from TS
  Vec state;
  BoutReal ts_time; ///< Internal PETSc timestepper time

  PetscLib lib; ///< Handles initialising, finalising PETSc

  Vec u{nullptr};                    ///< PETSc solution vector
  TS ts{nullptr};                    ///< PETSc timestepper object
  SNES snes{nullptr};                ///< PETSc nonlinear solver object
  KSP ksp{nullptr};                  ///< PETSc linear solver
  Mat Jmf{nullptr};                  ///< Matrix Free Jacobian
  Mat Jfd{nullptr};                  ///< Finite Difference Jacobian
  MatFDColoring fdcoloring{nullptr}; ///< Matrix coloring context

  BoutReal next_output; ///< When the monitor should be called next

  bool interpolate; ///< Interpolate to regular times?

  bool diagnose;    ///< If true, print some information about current stage
  bool user_precon; ///< Use user-supplied preconditioning function?

  BoutReal atol; ///< Absolute tolerance
  BoutReal rtol; ///< Relative tolerance
  BoutReal stol; ///< Convergence tolerance

  int maxnl; ///< Maximum nonlinear iterations per SNES solve
  int maxf; ///< Maximum number of function evaluations allowed in the solver (default: 10000)
  int maxl; ///< Maximum linear iterations

  std::string ts_type;          ///< PETSc TS time solver type
  std::string adapt_type;       ///< TSAdaptType timestep adaptation
  std::string snes_type;        ///< PETSc SNES nonlinear solver type
  std::string ksp_type;         ///< PETSc KSP linear solver type
  std::string pc_type;          ///< Preconditioner type
  std::string pc_hypre_type;    ///< Hypre preconditioner type
  std::string line_search_type; ///< Line search type

  bool matrix_free;          ///< Use matrix free Jacobian
  bool matrix_free_operator; ///< Use matrix free Jacobian in the operator?
  int lag_jacobian;          ///< Re-use Jacobian
  bool use_coloring;         ///< Use matrix coloring
  void updateColoring();     ///< Updates the coloring using Jfd

  bool kspsetinitialguessnonzero; ///< Set initial guess to non-zero

  BoutReal start_timestep;
  int mxstep;
};

#endif // BOUT_HAS_PETSC

#endif // BOUT_PETSC_SOLVER_H
