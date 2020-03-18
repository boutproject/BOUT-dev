/**************************************************************************
 * Interface to SUNDIALS CVODE
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

#ifndef __SUNDIAL_SOLVER_H__
#define __SUNDIAL_SOLVER_H__

#ifdef BOUT_HAS_CVODE

#include "bout_types.hxx"
#include "bout/solver.hxx"
#include "bout/solverfactory.hxx"

#include <sundials/sundials_config.h>
#if SUNDIALS_VERSION_MAJOR >= 3
#include <sunlinsol/sunlinsol_spgmr.h>
#endif

#if SUNDIALS_VERSION_MAJOR >= 4
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#endif

#include <nvector/nvector_parallel.h>

#include <vector>

class CvodeSolver;
class Options;

namespace {
RegisterSolver<CvodeSolver> registersolvercvode("cvode");
}

class CvodeSolver : public Solver {
public:
  CvodeSolver(Options* opts = nullptr);
  ~CvodeSolver();

  void setJacobian(Jacobian j) override { jacfunc = j; }

  BoutReal getCurrentTimestep() override { return hcur; }

  int init(int nout, BoutReal tstep) override;

  int run() override;
  BoutReal run(BoutReal tout);

  void resetInternalFields() override;

  // These functions used internally (but need to be public)
  void rhs(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void pre(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal* udata, BoutReal* rvec,
           BoutReal* zvec);
  void jac(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata);

private:
  int NOUT;          // Number of outputs. Specified in init, needed in run
  BoutReal TIMESTEP; // Time between outputs
  BoutReal hcur;     // Current internal timestep

  Jacobian jacfunc{nullptr}; // Jacobian - vector function
  bool diagnose{false};      // Output additional diagnostics

  N_Vector uvec{nullptr};   // Values
  void* cvode_mem{nullptr}; // CVODE internal memory block

  BoutReal pre_Wtime{0.0}; // Time in preconditioner
  int pre_ncalls{0};       // Number of calls to preconditioner

  bool cvode_initialised = false;

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

#endif // BOUT_HAS_CVODE
#endif // __SUNDIAL_SOLVER_H__
