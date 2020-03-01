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

namespace {
RegisterSolver<SplitRK> registersolversplitrk("splitrk");
}

class SplitRK : public Solver {
public:
  explicit SplitRK(Options *opt = nullptr) : Solver(opt) {}
  ~SplitRK() = default;

  int init(int nout, BoutReal tstep) override;

  int run() override;
private:
  int nstages{2}; ///< Number of stages in the RKL 
  
  BoutReal out_timestep{0.0}; ///< The output timestep
  int nsteps{0}; ///< Number of output steps
  
  BoutReal timestep{0.0}; ///< The internal timestep

  bool adaptive{true};   ///< Adapt timestep using tolerances?
  BoutReal atol{1e-10};  ///< Absolute tolerance
  BoutReal rtol{1e-5};   ///< Relative tolerance
  BoutReal max_timestep{1.0}; ///< Maximum timestep
  BoutReal max_timestep_change{2.0};  ///< Maximum factor by which the timestep should be changed
  int mxstep{1000};      ///< Maximum number of internal steps between outputs
  int adapt_period{1};   ///< Number of steps between checks

  bool diagnose{false};  ///< Turn on diagnostic output
  
  int nlocal{0}, neq{0}; ///< Number of variables on local processor and in total
  
  /// System state
  Array<BoutReal> state;
  
  /// Temporary time-stepping arrays
  /// These are used by both diffusion and advection time-step routines 
  Array<BoutReal> u1, u2, u3, dydt;

  /// Arrays used for adaptive timestepping
  Array<BoutReal> state1, state2;
  
  /// Take a combined step
  /// Uses 2nd order Strang splitting
  ///
  /// Note: start and result can be the same
  void take_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start,
                 Array<BoutReal>& result);

  /// Take a step of the diffusion terms
  /// Uses the Runge-Kutta-Legendre 2nd order method
  ///
  /// Note: start and result can be the same
  void take_diffusion_step(BoutReal curtime, BoutReal dt,
                                      Array<BoutReal>& start, Array<BoutReal>& result);

  /// Take a step of the advection terms
  /// Uses the Strong Stability Preserving Runge-Kutta 3rd order method
  ///
  /// Note: start and result can be the same
  void take_advection_step(BoutReal curtime, BoutReal dt, Array<BoutReal>& start,
                           Array<BoutReal>& result);
};

#endif
