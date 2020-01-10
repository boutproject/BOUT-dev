/**************************************************************************
 * 3rd-order Runge Kutta Strong Stability Preserving (SSP)
 * 
 * http://www.cscamm.umd.edu/tadmor/pub/linear-stability/Gottlieb-Shu-Tadmor.SIREV-01.pdf
 * 
 * Sigal Gottlieb, Chi-Wang Shu, and Eitan Tadmor, Strong stability-preserving 
 * high-order time discretization methods,
 * SIAM Rev. 43 (2001), no. 1, 89-112 (electronic). MR 2002f:65132
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

class RK3SSP;

#ifndef __RK3SSP_SOLVER_H__
#define __RK3SSP_SOLVER_H__

#include "mpi.h"

#include <bout_types.hxx>
#include <bout/solver.hxx>

namespace {
RegisterSolver<RK3SSP> registersolverrk3ssp("rk3ssp");
}

class RK3SSP : public Solver {
 public:
  RK3SSP(Options *opt = nullptr);
  ~RK3SSP(){};
  
  void setMaxTimestep(BoutReal dt) override;
  BoutReal getCurrentTimestep() override {return timestep; }
  
  int init(int nout, BoutReal tstep) override;
  
  int run() override;
 private:

  BoutReal max_timestep; // Maximum timestep
  int mxstep; // Maximum number of internal steps between outputs
  
  Array<BoutReal> f;
  
  BoutReal out_timestep; // The output timestep
  int nsteps; // Number of output steps
  
  BoutReal timestep; // The internal timestep

  int nlocal, neq; // Number of variables on local processor and in total
  
  void take_step(BoutReal curtime, BoutReal dt, 
                 Array<BoutReal> &start, Array<BoutReal> &result); // Take a single step to calculate f1
  
  Array<BoutReal> u1, u2, u3, L; // Time-stepping arrays
  
};

#endif // __RK4_SOLVER_H__

