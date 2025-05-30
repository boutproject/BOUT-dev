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

#ifndef BOUT_RK4_SOLVER_H
#define BOUT_RK4_SOLVER_H

#include "mpi.h"

#include <bout/bout_types.hxx>
#include <bout/solver.hxx>

namespace {
RegisterSolver<RK4Solver> registersolverrk4("rk4");
}

class RK4Solver : public Solver {
public:
  explicit RK4Solver(Options* opts = nullptr);

  void resetInternalFields() override;
  void setMaxTimestep(BoutReal dt) override;
  BoutReal getCurrentTimestep() override { return timestep; }

  int init() override;
  int run() override;

private:
  BoutReal atol, rtol;   //< Tolerances for adaptive timestepping
  BoutReal max_timestep; //< Maximum timestep
  BoutReal timestep;     //< The internal timestep
  int mxstep;            //< Maximum number of internal steps between outputs
  bool adaptive;         //< Adapt timestep?

  Array<BoutReal> f0, f1, f2;

  int nlocal, neq; //< Number of variables on local processor and in total

  /// Take a single step to calculate f1
  void take_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start,
                 Array<BoutReal>& result);

  Array<BoutReal> k1, k2, k3, k4, k5; //< Time-stepping arrays
};

#endif // BOUT_RK4_SOLVER_H
