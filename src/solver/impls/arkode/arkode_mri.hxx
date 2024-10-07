/**************************************************************************
 * Interface to ARKODE MRI solver
 * NOTE: ARKode is currently in beta testing so use with cautious optimism
 *
 * NOTE: Only one solver can currently be compiled in
 *
 **************************************************************************
 * Copyright 2010-2024 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#ifndef BOUT_ARKODE_MRI_SOLVER_H
#define BOUT_ARKODE_MRI_SOLVER_H

#include "bout/build_config.hxx"
#include "bout/solver.hxx"

#if not BOUT_HAS_ARKODE

namespace {
RegisterUnavailableSolver
    registerunavailablearkodemri("arkode_mri", "BOUT++ was not configured with ARKODE/SUNDIALS");
}

#else

#include "bout/bout_enum_class.hxx"
#include "bout/bout_types.hxx"
#include "bout/sundials_backports.hxx"

#include <nvector/nvector_parallel.h>
#include <sundials/sundials_config.h>
#include <arkode/arkode_mristep.h>

#if SUNDIALS_CONTROLLER_SUPPORT
#include <sundials/sundials_adaptcontroller.h>
#endif

#include <vector>

class ArkodeMRISolver;
class Options;

namespace {
RegisterSolver<ArkodeMRISolver> registersolverarkodemri("arkode_mri");
}

// enum describing treatment of equations
// Note: Capitalized because `explicit` is a C++ reserved keyword
BOUT_ENUM_CLASS(MRI_Treatment, ImEx, Implicit, Explicit);

// Adaptivity method
BOUT_ENUM_CLASS(MRI_AdapMethod, PID, PI, I, Explicit_Gustafsson, Implicit_Gustafsson,
                ImEx_Gustafsson);

class ArkodeMRISolver : public Solver {
public:
  explicit ArkodeMRISolver(Options* opts = nullptr);
  ~ArkodeMRISolver();

  BoutReal getCurrentTimestep() override { return hcur; }

  int init() override;

  int run() override;
  BoutReal run(BoutReal tout);

  // These functions used internally (but need to be public)
  void rhs_se(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void rhs_si(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void rhs_fe(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void rhs_fi(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void rhs_s(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void rhs_f(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void rhs(BoutReal t, BoutReal* udata, BoutReal* dudata);
  void pre_s(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal* udata, BoutReal* rvec,
           BoutReal* zvec);
  void pre_f(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal* udata, BoutReal* rvec,
           BoutReal* zvec);
  void jac_s(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata);
  void jac_f(BoutReal t, BoutReal* ydata, BoutReal* vdata, BoutReal* Jvdata);

private:
  BoutReal hcur; //< Current internal timestep

  bool diagnose{false}; //< Output additional diagnostics

  N_Vector uvec{nullptr};    //< Values
  void* arkode_mem{nullptr}; //< ARKODE internal memory block
  void* inner_arkode_mem{nullptr}; //< ARKODE internal memory block
  MRIStepInnerStepper inner_stepper{nullptr}; //< inner stepper

  BoutReal pre_Wtime_s{0.0}; //< Time in preconditioner
  BoutReal pre_Wtime_f{0.0}; //< Time in preconditioner
  int pre_ncalls_s{0};       //< Number of calls to preconditioner
  int pre_ncalls_f{0};       //< Number of calls to preconditioner

  /// Maximum number of steps to take between outputs
  int mxsteps;
  /// Integrator treatment enum: IMEX, Implicit or Explicit
  MRI_Treatment treatment;
  MRI_Treatment inner_treatment;
  /// Use linear implicit solver (only evaluates jacobian inversion once)
  bool set_linear;
  bool inner_set_linear;
  /// Solve explicit portion in fixed timestep mode. NOTE: This is not recommended except
  /// for code comparison
  bool fixed_step;
  /// Order of the internal step
  int order;
  /// Timestep adaptivity function
  MRI_AdapMethod adap_method;
  MRI_AdapMethod inner_adap_method;
  /// Absolute tolerance
  BoutReal abstol;
  /// Relative tolerance
  BoutReal reltol;
  /// Use separate absolute tolerance for each field
  bool use_vector_abstol;
  /// Maximum timestep (only used if greater than zero)
  bool use_precon;
  bool inner_use_precon;
  /// Number of Krylov basis vectors to use
  int maxl;
  int inner_maxl;
  /// Use right preconditioning instead of left preconditioning
  bool rightprec;

  // Diagnostics from ARKODE MRI
  int nsteps{0};
  int nfe_evals{0};
  int nfi_evals{0};
  int nniters{0};
  int npevals{0};
  int nliters{0};
  int inner_nsteps{0};
  int inner_nfe_evals{0};
  int inner_nfi_evals{0};
  int inner_nniters{0};
  int inner_npevals{0};
  int inner_nliters{0};

  void set_abstol_values(BoutReal* abstolvec_data, std::vector<BoutReal>& f2dtols,
                         std::vector<BoutReal>& f3dtols);
  void loop_abstol_values_op(Ind2D i2d, BoutReal* abstolvec_data, int& p,
                             std::vector<BoutReal>& f2dtols,
                             std::vector<BoutReal>& f3dtols, bool bndry);

  /// SPGMR solver structure
  SUNLinearSolver sun_solver{nullptr};
  SUNLinearSolver inner_sun_solver{nullptr};
  /// Solver for implicit stages
  SUNNonlinearSolver nonlinear_solver{nullptr};
  SUNNonlinearSolver inner_nonlinear_solver{nullptr};
#if SUNDIALS_CONTROLLER_SUPPORT
  /// Timestep controller
  SUNAdaptController controller{nullptr};
  SUNAdaptController inner_controller{nullptr};
#endif
  /// Context for SUNDIALS memory allocations
  sundials::Context suncontext;
};

#endif // BOUT_HAS_ARKODE
#endif // BOUT_ARKODE_MRI_SOLVER_H
