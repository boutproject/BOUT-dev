/**************************************************************************
 * 4th-order Runge Kutta explicit method with adaptive timestepping
 * 
 * Always available, since doesn't depend on external library
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

class RK4Solver;

#ifndef __RK4_SOLVER_H__
#define __RK4_SOLVER_H__

#include "mpi.h"

#include "bout_types.h"
#include "field2d.h"
#include "field3d.h"
#include "vector2d.h"
#include "vector3d.h"

#include "solver.h"

class RK4Solver : public Solver {
 public:
  RK4Solver();
  ~RK4Solver();

  int init(rhsfunc f, int argc, char **argv, bool restarting, int nout, BoutReal tstep);
  
  int run(MonitorFunc f);
 private:
  BoutReal atol, rtol; // Tolerances for adaptive timestepping
  BoutReal max_timestep; // Maximum timestep
  BoutReal start_timestep; // Starting timestep
  int mxstep; // Maximum number of internal steps between outputs
  
  BoutReal *f0, *f1, *f2;
  
  BoutReal out_timestep; // The output timestep
  int nsteps; // Number of output steps
  
  BoutReal timestep; // The internal timestep
  
  int nlocal; // Number of variables on local processor
  
  void take_step(BoutReal curtime, BoutReal dt, 
                 BoutReal *start, BoutReal *result); // Take a single step to calculate f1
};

#endif // __KARNIADAKIS_SOLVER_H__

