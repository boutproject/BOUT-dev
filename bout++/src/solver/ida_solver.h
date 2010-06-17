/**************************************************************************
 * Interface to SUNDIALS IDA
 * 
 * Solver for DAE systems (so can handle constraints)
 *
 * NOTE: Only one solver can currently be compiled in
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

class Solver;

#ifndef __IDA_SOLVER_H__
#define __IDA_SOLVER_H__

#include "bout_types.h"
#include "field2d.h"
#include "field3d.h"
#include "vector2d.h"
#include "vector3d.h"

#include "generic_solver.h"

// NOTE: MPI must be included before SUNDIALS, otherwise complains
#include "mpi.h"

#include <nvector/nvector_parallel.h>

#include <vector>
using std::vector;

enum SOLVER_VAR_OP {LOAD_VARS, LOAD_DERIVS, SET_ID, SAVE_VARS, SAVE_DERIVS};

class Solver : public GenericSolver {
 public:
  Solver();
  ~Solver();

  void setPrecon(PhysicsPrecon f) {prefunc = f;}
  
  int init(rhsfunc f, int argc, char **argv, bool restarting, int nout, real tstep);
  
  int run(MonitorFunc f);
  real run(real tout, int &ncalls, real &rhstime);

  // These functions used internally (but need to be public)
  void res(real t, real *udata, real *dudata, real *rdata);
  void pre(real t, real cj, real delta, real *udata, real *rvec, real *zvec);
 private:
  int NOUT; // Number of outputs. Specified in init, needed in run
  real TIMESTEP; // Time between outputs
  
  rhsfunc func; // RHS function
  PhysicsPrecon prefunc; // Preconditioner
  
  N_Vector uvec, duvec, id; // Values, time-derivatives, and equation type
  void *idamem;

  // Loading data from BOUT++ to/from IDA
  void loop_vars_op(int jx, int jy, real *udata, int &p, SOLVER_VAR_OP op);
  void loop_vars(real *udata, SOLVER_VAR_OP op);

  void load_vars(real *udata);
  void load_derivs(real *udata);
  void set_id(real *udata);
  int save_vars(real *udata);
  void save_derivs(real *dudata);

  real pre_Wtime; // Time in preconditioner
  real pre_ncalls; // Number of calls to preconditioner
};

#endif // __IDA_SOLVER_H__

