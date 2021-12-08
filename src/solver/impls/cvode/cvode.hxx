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

#include "bout/build_config.hxx"
#include "bout/solver.hxx"

#if not BOUT_HAS_CVODE

namespace {
RegisterUnavailableSolver registerunavailablecvode("cvode",
                                                   "BOUT++ was not configured with CVODE/SUNDIALS");
}

#else

#include "bout_types.hxx"

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

  bool diagnose{false};      // Output additional diagnostics

  N_Vector uvec{nullptr};   // Values
  void* cvode_mem{nullptr}; // CVODE internal memory block

  BoutReal pre_Wtime{0.0}; // Time in preconditioner
  int pre_ncalls{0};       // Number of calls to preconditioner

  // Diagnostics from CVODE
  int nsteps{0};
  int nfevals{0};
  int nniters{0};
  int npevals{0};
  int nliters{0};
  BoutReal last_step{0.0};
  int last_order{0};
  int num_fails{0};
  int nonlin_fails{0};
  int stab_lims{0};

  bool cvode_initialised = false;

  void set_vector_option_values(BoutReal* option_data, std::vector<BoutReal>& f2dtols,
                                std::vector<BoutReal>& f3dtols);
  void loop_vector_option_values_op(Ind2D i2d, BoutReal* option_data, int& p,
                                    std::vector<BoutReal>& f2dtols,
                                    std::vector<BoutReal>& f3dtols, bool bndry);
  template<class FieldType>
  std::vector<BoutReal> create_constraints(const std::vector<VarStr<FieldType>>& fields);
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
