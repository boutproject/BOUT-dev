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

//Adding support for 3.4 based on BOUT++ 3.3 code <DD>

#ifdef BOUT_HAS_PETSC_3_4
class PetscSolver;

#ifndef __PETSC_SOLVER_H__
#define __PETSC_SOLVER_H__

#include <petsc.h>

#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/vector2d.hxx>
#include <bout/vector3d.hxx>

#include <bout/solver.hxx>

#include <bout/petsclib.hxx>

#include <vector>

typedef PetscScalar BoutReal;
typedef PetscInt integer;
typedef PetscBool boole;
#define OPT_SIZE 40

using std::vector;

extern BoutReal simtime;
extern PetscErrorCode PetscMonitor(TS,PetscInt,PetscReal,Vec,void *ctx);
extern PetscErrorCode PetscSNESMonitor(SNES,PetscInt,PetscReal,void *ctx);
extern int jstruc(int NVARS, int NXPE, int MXSUB, int NYPE, int MYSUB, int MZ, int MYG, int MXG);

extern PetscErrorCode solver_ijacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);

typedef struct snes_info {
  PetscInt it;
  PetscInt linear_its;
  PetscReal time;
  PetscReal norm;
} snes_info;

class PetscSolver : public Solver {
 public:
  PetscSolver(Options *opt = NULL);
  ~PetscSolver();

  // Can be called from physics initialisation to supply callbacks
  void setJacobian(Jacobian j) {jacfunc = j; }

  int init(int NOUT, BoutReal TIMESTEP) override;

  int run() override;
  PetscErrorCode run(MonitorFunc mon);

  // These functions used internally (but need to be public)

  PetscErrorCode rhs(TS ts,PetscReal t,Vec globalin,Vec globalout);
  PetscErrorCode pre(PC pc, Vec x, Vec y);
  PetscErrorCode jac(Vec x, Vec y);
  friend PetscErrorCode PetscMonitor(TS,PetscInt,PetscReal,Vec,void *ctx);
  friend PetscErrorCode PetscSNESMonitor(SNES,PetscInt,PetscReal,void *ctx);
  friend PetscErrorCode solver_ijacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);

  PetscLogEvent solver_event, loop_event, init_event;
 private:
  PetscLib lib; // Handles initialize / finalize
  
  Jacobian jacfunc; // Jacobian - vector function

  BoutReal shift;   // Shift (alpha) parameter from TS
  Vec state;
  BoutReal ts_time;

  Vec           u;
  TS            ts;
  Mat           J,Jmf;
  MatFDColoring matfdcoloring;

  Options *options;

  int nout;   // The number of outputs
  BoutReal tstep; // Time between outputs

  bool diagnose;

  BoutReal next_output;  // When the monitor should be called next

  PetscBool interpolate; // Whether to interpolate or not

  char output_name[PETSC_MAX_PATH_LEN];
  PetscBool output_flag;
  PetscInt prev_linear_its;
  BoutReal bout_snes_time;
  vector<snes_info> snes_list;
};


#endif // __PETSC_SOLVER_H__

#endif
