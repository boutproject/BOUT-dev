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

#ifndef __IDA_SOLVER_H__
#define __IDA_SOLVER_H__

#include "bout/build_config.hxx"
#include "bout/solver.hxx"

#if not BOUT_HAS_IDA

namespace {
RegisterUnavailableSolver registerunavailableida("ida",
                                                 "BOUT++ was not configured with IDA/SUNDIALS");
}

#else

#include "bout_types.hxx"

#include <sundials/sundials_config.h>
#if SUNDIALS_VERSION_MAJOR >= 3
#include <sunlinsol/sunlinsol_spgmr.h>
#endif

#include <nvector/nvector_parallel.h>

class IdaSolver;
class Options;

namespace {
RegisterSolver<IdaSolver> registersolverida("ida");
}

class IdaSolver : public Solver {
public:
  IdaSolver(Options* opts = nullptr);
  ~IdaSolver();

  int init(int nout, BoutReal tstep) override;

  int run() override;
  BoutReal run(BoutReal tout);

  // These functions used internally (but need to be public)
  void res(BoutReal t, BoutReal* udata, BoutReal* dudata, BoutReal* rdata);
  void pre(BoutReal t, BoutReal cj, BoutReal delta, BoutReal* udata, BoutReal* rvec,
           BoutReal* zvec);

private:
  int NOUT;          // Number of outputs. Specified in init, needed in run
  BoutReal TIMESTEP; // Time between outputs

  N_Vector uvec{nullptr};  // Values
  N_Vector duvec{nullptr}; // Time-derivatives
  N_Vector id{nullptr};    // Equation type
  void* idamem{nullptr};   // IDA internal memory block

  BoutReal pre_Wtime{0.0}; // Time in preconditioner
  int pre_ncalls{0};       // Number of calls to preconditioner

#if SUNDIALS_VERSION_MAJOR >= 3
  /// SPGMR solver structure
  SUNLinearSolver sun_solver{nullptr};
#endif
};

#endif // BOUT_HAS_IDA
#endif // __IDA_SOLVER_H__
