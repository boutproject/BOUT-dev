/**************************************************************************
 * Interface to PETSc 3.1 solver
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

#ifndef BOUT_HAS_PETSC_3_1

#include "emptysolver.hxx"
typedef EmptySolver Petsc31Solver;

#else
class Petsc31Solver;

#ifndef __PETSC31_SOLVER_H__
#define __PETSC31_SOLVER_H__

#include <petscts.h>

#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

#include "solver.hxx"

#include <vector>

typedef PetscScalar BoutReal;
typedef PetscInt integer;
typedef PetscTruth boole;
#define OPT_SIZE 40

using std::vector;

typedef int (*rhsfunc)(BoutReal);

EXTERN PetscErrorCode PreStep(TS);
EXTERN PetscErrorCode PostStep(TS);
EXTERN int jstruc(int NVARS, int NXPE, int MXSUB, int NYPE, int MYSUB, int MZ, int MYG, int MXG);

class Petsc31Solver : public Solver {
 public:
  Petsc31Solver();
  ~Petsc31Solver();
  
  int setup(int argc, char **argv);
  
  int init(rhsfunc f, int argc, char **argv, bool restarting, int NOUT, BoutReal TIMESTEP);

  int run(MonitorFunc f);

  // These functions used internally (but need to be public)
  PetscErrorCode rhs(TS ts,PetscReal t,Vec globalin,Vec globalout);
  friend PetscErrorCode PreStep(TS);
  friend PetscErrorCode PostStep(TS);

 private:
  Vec           u;
  TS            ts;
  Mat           J;
  MatFDColoring matfdcoloring;

  int nout;   // The number of outputs
  BoutReal tstep; // Time between outputs
  MonitorFunc monitor; // Monitor function to call regularly

  BoutReal next_time;  // When the monitor should be called next
  bool outputnext; // true if the monitor should be called next time
  
  PetscLogEvent USER_EVENT;
};


#endif // __PETSC31_SOLVER_H__

#endif
