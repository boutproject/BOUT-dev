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

#ifndef BOUT_SUNDIAL_SOLVER_H
#define BOUT_SUNDIAL_SOLVER_H

#include "bout/build_defines.hxx"
#include "bout/solver.hxx"

#if not BOUT_HAS_CVODE

namespace {
RegisterUnavailableSolver
    registerunavailablecvode("cvode", "BOUT++ was not configured with CVODE/SUNDIALS");
}

#else

#include "bout/bout_types.hxx"
#include "bout/region.hxx"
#include "bout/sundials_backports.hxx"

#include <string>
#include <vector>

class CvodeSolver;
class Options;

namespace {
RegisterSolver<CvodeSolver> registersolvercvode("cvode");
}

class CvodeSolver : public Solver {
public:
  explicit CvodeSolver(Options* opts = nullptr);
  ~CvodeSolver() override;

  BoutReal getCurrentTimestep() override { return hcur; }

  int init() override;
  int run() override;
  BoutReal run(BoutReal tout);

  void resetInternalFields() override;

  // These functions are used internally (but need to be public)
  void rhs(BoutReal t, BoutReal* udata, BoutReal* dudata, bool linear);
  void pre(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal* udata, BoutReal* rvec,
           BoutReal* zvec);
  void jac(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata);

private:
  BoutReal hcur; //< Current internal timestep

  bool diagnose{false}; //< Output additional diagnostics

  N_Vector uvec{nullptr};   //< Values
  void* cvode_mem{nullptr}; //< CVODE internal memory block

  BoutReal pre_Wtime{0.0}; //< Time in preconditioner
  int pre_ncalls{0};       //< Number of calls to preconditioner

  /// Use Adams Moulton implicit multistep. Otherwise BDF method
  bool adams_moulton;
  /// Use functional iteration instead of Newton
  bool func_iter;
  /// Maximum order of method to use. < 0 means no limit
  int max_order;
  bool stablimdet;
  /// Absolute tolerance
  BoutReal abstol;
  /// Relative tolerance
  BoutReal reltol;
  /// Use separate absolute tolerance for each field
  bool use_vector_abstol;
  /// Maximum number of internal steps between outputs.
  int mxsteps;
  /// Maximum time step size
  BoutReal max_timestep;
  /// Minimum time step size
  BoutReal min_timestep;
  /// Starting time step. < 0 then chosen by CVODE.
  BoutReal start_timestep;
  /// Maximum order
  int mxorder;
  /// Maximum number of nonlinear iterations allowed by CVODE before
  /// reducing timestep. CVODE default (used if this option is
  /// negative) is 3
  int max_nonlinear_iterations;
  /// Use CVODE function CVodeSetConstraints to constrain variables -
  /// the constraint to be applied is set by the positivity_constraint
  /// option in the subsection for each variable
  bool apply_positivity_constraints;
  /// Maximum number of linear iterations
  int maxl;
  /// Use preconditioner?
  bool use_precon;
  /// Use right preconditioner? Otherwise use left.
  bool rightprec;
  bool use_jacobian;
  BoutReal cvode_nonlinear_convergence_coef;
  BoutReal cvode_linear_convergence_coef;

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

  bool cvode_initialised{false};

  void set_vector_option_values(BoutReal* option_data, std::vector<BoutReal>& f2dtols,
                                std::vector<BoutReal>& f3dtols);
  void loop_vector_option_values_op(Ind2D i2d, BoutReal* option_data, int& p,
                                    std::vector<BoutReal>& f2dtols,
                                    std::vector<BoutReal>& f3dtols, bool bndry);
  template <class FieldType>
  std::vector<BoutReal> create_constraints(const std::vector<VarStr<FieldType>>& fields);

  /// SPGMR solver structure
  SUNLinearSolver sun_solver{nullptr};
  /// Solver for functional iterations for Adams-Moulton
  SUNNonlinearSolver nonlinear_solver{nullptr};
  /// Context for SUNDIALS memory allocations
  sundials::Context suncontext;
};

#endif // BOUT_HAS_CVODE
#endif // BOUT_SUNDIAL_SOLVER_H
