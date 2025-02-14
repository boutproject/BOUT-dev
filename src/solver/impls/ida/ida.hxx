/**************************************************************************
 * Interface to SUNDIALS IDA
 *
 * IdaSolver for DAE systems (so can handle constraints)
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

#ifndef BOUT_IDA_SOLVER_H
#define BOUT_IDA_SOLVER_H

#include "bout/build_defines.hxx"
#include "bout/solver.hxx"

#if not BOUT_HAS_IDA

namespace {
RegisterUnavailableSolver
    registerunavailableida("ida", "BOUT++ was not configured with IDA/SUNDIALS");
}

#else

#include "bout/bout_types.hxx"
#include "bout/sundials_backports.hxx"

#include <string>

class IdaSolver;
class Options;

namespace {
RegisterSolver<IdaSolver> registersolverida("ida");
}

class IdaSolver : public Solver {
public:
  explicit IdaSolver(Options* opts = nullptr);
  ~IdaSolver() override;

  int init() override;
  int run() override;
  BoutReal run(BoutReal tout);

  // These functions used internally (but need to be public)
  void res(BoutReal t, BoutReal* udata, BoutReal* dudata, BoutReal* rdata);
  void pre(BoutReal t, BoutReal cj, BoutReal delta, BoutReal* udata, BoutReal* rvec,
           BoutReal* zvec);

private:
  /// Absolute tolerance
  BoutReal abstol;
  /// Relative tolerance
  BoutReal reltol;
  /// Maximum number of steps to take between outputs
  int mxsteps;
  /// Use user-supplied preconditioner
  bool use_precon;
  /// Correct the initial values
  bool correct_start;

  N_Vector uvec{nullptr};  // Values
  N_Vector duvec{nullptr}; // Time-derivatives
  N_Vector id{nullptr};    // Equation type
  void* idamem{nullptr};   // IDA internal memory block

  BoutReal pre_Wtime{0.0}; // Time in preconditioner
  int pre_ncalls{0};       // Number of calls to preconditioner

  /// SPGMR solver structure
  SUNLinearSolver sun_solver{nullptr};
  /// Context for SUNDIALS memory allocations
  sundials::Context suncontext;
};

#endif // BOUT_HAS_IDA
#endif // BOUT_IDA_SOLVER_H
