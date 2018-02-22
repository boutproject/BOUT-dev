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

#ifdef BOUT_HAS_PETSC_3_1
class PetscSolver;

#ifndef __PETSC_SOLVER_H__
#define __PETSC_SOLVER_H__

#include <petscts.h>

#include <field2d.hxx>
#include <field3d.hxx>
#include <vector2d.hxx>
#include <vector3d.hxx>

#include <bout/solver.hxx>

#include <bout/petsclib.hxx>

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

class PetscSolver : public Solver {
 public:
  PetscSolver();
  ~PetscSolver();
  
  int init(int NOUT, BoutReal TIMESTEP) override;

  int run() override;

  // These functions used internally (but need to be public)
  PetscErrorCode rhs(TS ts,PetscReal t,Vec globalin,Vec globalout);
  friend PetscErrorCode PreStep(TS);
  friend PetscErrorCode PostStep(TS);

 private:
  PetscLib lib; // Handles initialisation and finalisation
  
  Vec           u;
  TS            ts;
  Mat           J;
  MatFDColoring matfdcoloring;

  int nout;   // The number of outputs
  BoutReal tstep; // Time between outputs

  BoutReal next_time;  // When the monitor should be called next
  bool outputnext; // true if the monitor should be called next time
};


#endif // __PETSC31_SOLVER_H__

#endif
