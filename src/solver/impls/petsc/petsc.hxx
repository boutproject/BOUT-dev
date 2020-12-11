/**************************************************************************
 * Interface to PETSc solver
 * NOTE: This class needs tidying, generalising to use FieldData interface
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

#ifndef __PETSC_SOLVER_H__
#define __PETSC_SOLVER_H__

#include "bout/build_config.hxx"
#include "bout/solver.hxx"

#if not BOUT_HAS_PETSC

namespace {
RegisterUnavailableSolver registerunavailablepetsc("petsc",
                                                   "BOUT++ was not configured with PETSc");
}

#else

class PetscSolver;

#include <field2d.hxx>
#include <field3d.hxx>
#include <vector2d.hxx>
#include <vector3d.hxx>

#include <petsc.h>
// PETSc creates macros for MPI calls, which interfere with the MpiWrapper class
#undef MPI_Allreduce

#include <bout/petsclib.hxx>

#include <vector>

namespace {
RegisterSolver<PetscSolver> registersolverpetsc("petsc");
}

using BoutReal = PetscScalar;
#define OPT_SIZE 40

using rhsfunc = int (*)(BoutReal);

extern BoutReal simtime;

/// Monitor function called on every internal timestep
extern PetscErrorCode PetscMonitor(TS, PetscInt, PetscReal, Vec, void *ctx);
/// Monitor function for SNES
extern PetscErrorCode PetscSNESMonitor(SNES, PetscInt, PetscReal, void *ctx);

/// Compute IJacobian = dF/dU + a dF/dUdot  - a dummy matrix used for pc=none
#if PETSC_VERSION_GE(3, 5, 0)
extern PetscErrorCode solver_ijacobian(TS, PetscReal, Vec, Vec, PetscReal, Mat, Mat,
                                       void *);
#else
extern PetscErrorCode solver_ijacobian(TS, PetscReal, Vec, Vec, PetscReal, Mat *, Mat *,
                                       MatStructure *, void *);
#endif

/// Data for SNES
struct snes_info {
  PetscInt it;
  PetscInt linear_its;
  PetscReal time;
  PetscReal norm;
};

class PetscSolver : public Solver {
public:
  PetscSolver(Options *opts = nullptr);
  ~PetscSolver();

  // Can be called from physics initialisation to supply callbacks
  void setPrecon(PhysicsPrecon f) { prefunc = f; }
  void setJacobian(Jacobian j) override { jacfunc = j; }

  int init(int NOUT, BoutReal TIMESTEP) override;

  int run() override;

  // These functions used internally (but need to be public)

  /// Wrapper for the RHS function
  PetscErrorCode rhs(TS ts, PetscReal t, Vec globalin, Vec globalout);
  /// Wrapper for the preconditioner
  PetscErrorCode pre(PC pc, Vec x, Vec y);
  /// Wrapper for the Jacobian function
  PetscErrorCode jac(Vec x, Vec y);

  // Call back functions that need to access internal state
  friend PetscErrorCode PetscMonitor(TS, PetscInt, PetscReal, Vec, void *ctx);
  friend PetscErrorCode PetscSNESMonitor(SNES, PetscInt, PetscReal, void *ctx);
#if PETSC_VERSION_GE(3, 5, 0)
  friend PetscErrorCode solver_ijacobian(TS, PetscReal, Vec, Vec, PetscReal, Mat, Mat,
                                         void *);
#else
  friend PetscErrorCode solver_ijacobian(TS, PetscReal, Vec, Vec, PetscReal, Mat *, Mat *,
                                         MatStructure *, void *);
#endif

  PetscLogEvent solver_event, loop_event, init_event;

private:
  PhysicsPrecon prefunc; ///< Preconditioner
  Jacobian jacfunc;      ///< Jacobian - vector function

  BoutReal shift; ///< Shift (alpha) parameter from TS
  Vec state;
  BoutReal ts_time; ///< Internal PETSc timestepper time

  PetscLib lib; ///< Handles initialising, finalising PETSc

  Vec u;      ///< PETSc solution vector
  TS ts;      ///< PETSc timestepper object
  Mat J, Jmf; ///< RHS Jacobian
  MatFDColoring matfdcoloring;

  int nout;       ///< The number of outputs
  BoutReal tstep; ///< Time between outputs

  bool diagnose; ///< If true, print some information about current stage

  BoutReal next_output; ///< When the monitor should be called next

  PetscBool interpolate; ///< Whether to interpolate or not

  char output_name[PETSC_MAX_PATH_LEN];
  PetscBool output_flag;
  PetscInt prev_linear_its;
  BoutReal bout_snes_time;
  std::vector<snes_info> snes_list;

  bool adaptive; ///< Use adaptive timestepping
};

#endif // BOUT_HAS_PETSC

#endif // __PETSC_SOLVER_H__
