/**************************************************************************
 * Interface to PVODE solver
 * NOTE: This class needs tidying, generalising to use FieldData interface
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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

#ifndef BOUT_HAS_PVODE

#include "emptysolver.hxx"
typedef EmptySolver PvodeSolver;
 
#else
class PvodeSolver;

#ifndef __PVODE_SOLVER_H__
#define __PVODE_SOLVER_H__

#include "globals.hxx"

#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

#include "solver.hxx"
#include "bout_types.hxx"
#include "interpolation.hxx"   // Cell interpolation

#include <mpi.h>             // MPI data types and prototypes
#include <pvode/nvector.h>
#include <pvode/cvode.h>     // main CVODE header file
#include <pvode/iterativ.h>  // contains the enum for types of preconditioning
#include <pvode/cvspgmr.h>   // use CVSPGMR linear solver each internal step
#include <pvode/pvbbdpre.h>  // band preconditioner function prototypes

#include <stdlib.h>
#include <vector>

using std::vector;


class PvodeSolver : public Solver {
 public:
  PvodeSolver();
  ~PvodeSolver();

  void setPrecon(PhysicsPrecon f) {} // Doesn't do much yet
  
  BoutReal getCurrentTimestep() { return hcur; }
  
  int init(rhsfunc f, int argc, char **argv, bool restarting, int nout, BoutReal tstep);
  
  int run(MonitorFunc f);
  BoutReal run(BoutReal tout, int &ncalls, BoutReal &rhstime);

  // These functions used internally (but need to be public)
  void rhs(int N, BoutReal t, BoutReal *udata, BoutReal *dudata);
  void gloc(int N, BoutReal t, BoutReal *udata, BoutReal *dudata);

 private:
  int NOUT; // Number of outputs. Specified in init, needed in run
  BoutReal TIMESTEP; // Time between outputs
  BoutReal hcur; // Current internal timestep
  
  pvode::N_Vector u;
  pvode::machEnvType machEnv;
  void *cvode_mem;
  
  rhsfunc func; // RHS function
  rhsfunc gfunc; // Preconditioner function
};

#endif // __PVODE_SOLVER_H__

#endif
