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
RegisterUnavailableSolver
    registerunavailablearkode("arkode", "BOUT++ was not configured with ARKODE/SUNDIALS");
}

#else

#include "bout/bout_types.hxx"
#include "bout/sundials_backports.hxx"

#include <nvector/nvector_parallel.h>
#include <sundials/sundials_config.h>

#include <vector>

class ArkodeSolver;
class Options;

namespace {
RegisterSolver<ArkodeSolver> registersolverarkode("arkode");
}

class ArkodeSolver : public Solver {
public:
  explicit ArkodeSolver(Options* opts = nullptr);
  ~ArkodeSolver();

  BoutReal getCurrentTimestep() override { return hcur; }

  int init() override;

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
  BoutReal hcur; //< Current internal timestep

  bool diagnose{false}; //< Output additional diagnostics

  N_Vector uvec{nullptr};    //< Values
  void* arkode_mem{nullptr}; //< ARKODE internal memory block

  BoutReal pre_Wtime{0.0}; //< Time in preconditioner
  int pre_ncalls{0};       //< Number of calls to preconditioner

  /// Maximum number of steps to take between outputs
  int mxsteps;
  /// Use ImEx capability
  bool imex;
  /// Solve only explicit part
  bool solve_explicit;
  /// Solve only implicit part
  bool solve_implicit;
  /// Use linear implicit solver (only evaluates jacobian inversion once)
  bool set_linear;
  /// Solve explicit portion in fixed timestep mode. NOTE: This is not recommended except
  /// for code comparison
  bool fixed_step;
  /// Order of internal step
  int order;
  /// Fraction of the estimated explicitly stable step to use
  BoutReal cfl_frac;
  /// Set timestep adaptivity function:
  /// - 0: PID adaptivity (default)
  /// - 1: PI
  /// - 2: I
  /// - 3: explicit Gustafsson
  /// - 4: implicit Gustafsson
  /// - 5: ImEx Gustafsson
  int adap_method;
  /// Absolute tolerance
  BoutReal abstol;
  /// Relative tolerance
  BoutReal reltol;
  /// Use separate absolute tolerance for each field
  bool use_vector_abstol;
  /// Maximum timestep (only used if greater than zero)
  BoutReal max_timestep;
  /// Minimum timestep (only used if greater than zero)
  BoutReal min_timestep;
  /// Initial timestep (only used if greater than zero)
  BoutReal start_timestep;
  /// Use accelerated fixed point solver instead of Newton iterative
  bool fixed_point;
  /// Use user-supplied preconditioner function
  bool use_precon;
  /// Number of Krylov basis vectors to use
  int maxl;
  /// Use right preconditioning instead of left preconditioning
  bool rightprec;
  /// Use user-supplied Jacobian function
  bool use_jacobian;
  /// Use ARKode optimal parameters
  bool optimize;

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

  /// SPGMR solver structure
  SUNLinearSolver sun_solver{nullptr};
  /// Solver for functional iterations for Adams-Moulton
  SUNNonlinearSolver nonlinear_solver{nullptr};
  /// Context for SUNDIALS memory allocations
  sundials::Context suncontext;
};

#endif // BOUT_HAS_ARKODE
#endif // __ARKODE_SOLVER_H__
