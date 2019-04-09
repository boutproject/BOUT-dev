/**************************************************************************
 * 
 *
 * Always available, since doesn't depend on external library
 * 
 **************************************************************************
 * Copyright 2019 Ben Dudson
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

class SplitRK;

#pragma once

#ifndef SPLITRK_HXX
#define SPLITRK_HXX

#include <bout_types.hxx>
#include <bout/solver.hxx>

#include <bout/solverfactory.hxx>
namespace {
RegisterSolver<SplitRK> registersolversplitrk("splitrk");
}

class SplitRK : public Solver {
public:
  SplitRK(Options *opt = nullptr) : Solver(opt) {}
  ~SplitRK() = default;

  int init(int nout, BoutReal tstep) override;

  int run() override;
private:
  int nstages; ///< Number of stages in the RKL 
  
  BoutReal out_timestep; ///< The output timestep
  int nsteps; ///< Number of output steps

  int ninternal_steps; ///< Number of internal timesteps
  BoutReal timestep; ///< The internal timestep

  /// System state
  Array<BoutReal> state;
  
  // Temporary time-stepping arrays
  Array<BoutReal> u1, u2, u3, L;

  /// Take a combined step
  /// Uses 2nd order Strang splitting
  Array<BoutReal> take_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start);

  /// Take a step of the diffusion terms
  /// Uses the Runge-Kutta-Legendre 2nd order method
  Array<BoutReal> take_diffusion_step(BoutReal curtime, BoutReal dt,
                                      Array<BoutReal>& start);

  /// Take a step of the advection terms
  /// Uses the Strong Stability Preserving Runge-Kutta 3rd order method
  Array<BoutReal> take_advection_step(BoutReal curtime, BoutReal dt,
                                      Array<BoutReal>& start);
};

#endif
