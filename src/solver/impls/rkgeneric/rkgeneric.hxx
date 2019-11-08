/**************************************************************************
 * Generic Runge Kutta explicit method with adaptive timestepping
 * 
 * Always available, since doesn't depend on external library
 * 
 **************************************************************************
 * Written by D Dickinson 2015
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

class RKGenericSolver;

#ifndef __RKGENERIC_SOLVER_H__
#define __RKGENERIC_SOLVER_H__

#include "mpi.h"

#include <bout_types.hxx>
#include <bout/solver.hxx>
#include <bout/rkscheme.hxx>

namespace {
RegisterSolver<RKGenericSolver> registersolverrkgeneric("rkgeneric");
}

class RKGenericSolver : public Solver {
 public:
  RKGenericSolver(Options *options);
  ~RKGenericSolver() = default;
  
  void resetInternalFields() override;
  
  //Utilities only used by the CTU bracket approach
  void setMaxTimestep(BoutReal dt) override;
  BoutReal getCurrentTimestep() override {return timestep; }

  //Setup solver and scheme
  int init(int nout, BoutReal tstep) override;

  //Actually evolve
  int run() override;

 private:
  //Take a step using the scheme
  BoutReal take_step(BoutReal timeIn,BoutReal dt, const Array<BoutReal> &start, 
		     Array<BoutReal> &resultFollow);

  //Used for storing current state and next step
  Array<BoutReal> f0, f2, tmpState;

  //Inputs
  BoutReal atol, rtol;   // Tolerances for adaptive timestepping
  BoutReal max_timestep; // Maximum timestep
  int mxstep; // Maximum number of internal steps between outputs
  bool adaptive;   // Adapt timestep?

  //Internal vars
  BoutReal out_timestep; // The output timestep
  int nsteps; // Number of output steps
  BoutReal timestep; // The internal timestep
  int nlocal, neq; // Number of variables on local processor and in total
  
  //Pointer to the actual scheme used
  std::unique_ptr<RKScheme> scheme{nullptr};

};

#endif // __RKGENERIC_SOLVER_H__

