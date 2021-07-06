/**************************************************************************
 * Power method for eigenvalue
 * 
 * Finds the largest (fastest growing) eigenvector and eigenvalue
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

class PowerSolver;

#ifndef __POWER_SOLVER_H__
#define __POWER_SOLVER_H__

#include <bout_types.hxx>
#include <bout/solver.hxx>

namespace {
RegisterSolver<PowerSolver> registersolverpower("power");
}

class Options;

class PowerSolver : public Solver {
 public:
  PowerSolver() : Solver() {}
  PowerSolver(Options*) : Solver() {}
  ~PowerSolver(){};
  
  int init(int nout, BoutReal tstep) override;
  
  int run() override;

  void outputVars(Options& output_options, bool save_repeat = true) override;
 private:

  BoutReal curtime; // Current simulation time (fixed)
  
  BoutReal eigenvalue;

  int nlocal, nglobal; // Number of variables
  Array<BoutReal> f0;  // The system state
  
  int nsteps; // Number of output steps
  
  BoutReal norm(Array<BoutReal> &state);
  void divide(Array<BoutReal> &in, BoutReal value);
};

#endif // __KARNIADAKIS_SOLVER_H__

