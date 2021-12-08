/**************************************************************************
 * Interface to ARKODE solver
 * NOTE: ARKode is currently in beta testing so use with cautious optimism
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

#ifndef __ARKODE_SOLVER_H__
#define __ARKODE_SOLVER_H__

#include "bout/build_config.hxx"
#include "bout/solver.hxx"

#if not BOUT_HAS_ARKODE

namespace {
RegisterUnavailableSolver registerunavailablearkode("arkode",
                                                    "BOUT++ was not configured with ARKODE/SUNDIALS");
}

#else

#include "bout_types.hxx"

#include <sundials/sundials_config.h>
#if SUNDIALS_VERSION_MAJOR >= 3
#include <sunlinsol/sunlinsol_spgmr.h>
#endif

#if SUNDIALS_VERSION_MAJOR >= 4
#include <sundials/sundials_nonlinearsolver.h>
#endif

#include <nvector/nvector_parallel.h>

#include <vector>

class ArkodeSolver;
class Options;

namespace {
RegisterSolver<ArkodeSolver> registersolverarkode("arkode");
}

class ArkodeSolver : public Solver {
public:
  ArkodeSolver(Options* opts = nullptr);
  ~ArkodeSolver();

  BoutReal getCurrentTimestep() override { return hcur; }

  int init(int nout, BoutReal tstep) override;

  int run() override;
  BoutReal run(BoutReal tout);

  // These functions used internally (but need to be public)
  void rhs_e(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void rhs_i(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void rhs(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void pre(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal* udata, BoutReal* rvec,
           BoutReal* zvec);
  void jac(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata);

private:
  int NOUT;          // Number of outputs. Specified in init, needed in run
  BoutReal TIMESTEP; // Time between outputs
  BoutReal hcur;     // Current internal timestep

  bool diagnose{false};      // Output additional diagnostics

  N_Vector uvec{nullptr};    // Values
  void* arkode_mem{nullptr}; // ARKODE internal memory block

  BoutReal pre_Wtime{0.0}; // Time in preconditioner
  int pre_ncalls{0};       // Number of calls to preconditioner

  // Diagnostics from ARKODE
  int nsteps{0};
  int nfe_evals{0};
  int nfi_evals{0};
  int nniters{0};
  int npevals{0};
  int nliters{0};

  void set_abstol_values(BoutReal* abstolvec_data, std::vector<BoutReal>& f2dtols,
                         std::vector<BoutReal>& f3dtols);
  void loop_abstol_values_op(Ind2D i2d, BoutReal* abstolvec_data, int& p,
                             std::vector<BoutReal>& f2dtols,
                             std::vector<BoutReal>& f3dtols, bool bndry);

#if SUNDIALS_VERSION_MAJOR >= 3
  /// SPGMR solver structure
  SUNLinearSolver sun_solver{nullptr};
#endif
#if SUNDIALS_VERSION_MAJOR >= 4
  /// Solver for functional iterations for Adams-Moulton
  SUNNonlinearSolver nonlinear_solver{nullptr};
#endif
};

#endif // BOUT_HAS_ARKODE
#endif // __ARKODE_SOLVER_H__
