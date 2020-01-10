/*!*************************************************************************
 * \file imex-bdf2.cxx
 *
 * 2nd order IMEX-BDF scheme
 *
 * Scheme taken from this paper: http://homepages.cwi.nl/~willem/DOCART/JCP07.pdf
 * W.Hundsdorfer, S.J.Ruuth "IMEX extensions of linear multistep methods with general
 * monotonicity and boundedness properties" JCP 225 (2007) 2016-2042
 *
 *
 * Uses PETSc for the SNES interface
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

#ifdef BOUT_HAS_PETSC

class IMEXBDF2;

#ifndef __IMEXBDF2_SOLVER_H__
#define __IMEXBDF2_SOLVER_H__

#include "mpi.h"

#include <bout_types.hxx>
#include <bout/solver.hxx>

#include <bout/petsclib.hxx>

#include <petsc.h>
#include <petscsnes.h>
// PETSc creates macros for MPI calls, which interfere with the MpiWrapper class
#undef MPI_Allreduce

namespace {
RegisterSolver<IMEXBDF2> registersolverimexbdf2("imexbdf2");
}

/// IMEX-BDF2 time integration solver
///
/// Scheme taken from this paper: http://homepages.cwi.nl/~willem/DOCART/JCP07.pdf
/// W.Hundsdorfer, S.J.Ruuth "IMEX extensions of linear multistep methods with general
/// monotonicity and boundedness properties" JCP 225 (2007) 2016-2042
///
/// The method has been extended to variable order, variable timestep,
/// and includes some adaptive capabilities
///
class IMEXBDF2 : public Solver {
 public:
  IMEXBDF2(Options *opt = nullptr);
  ~IMEXBDF2();

  /// Returns the current internal timestep
  BoutReal getCurrentTimestep() override {return timestep; }

  /// Initialise solver. Must be called once and only once
  ///
  /// @param[in] nout         Number of outputs
  /// @param[in] tstep        Time between outputs. NB: Not internal timestep
  int init(int nout, BoutReal tstep) override;

  /// Run the simulation
  int run() override;

  /// Nonlinear function. This is called by PETSc SNES object
  /// via a static C-style function. For implicit
  /// time integration this function calculates:
  ///
  ///     f = (x - gamma*G(x)) - rhs
  ///
  ///
  /// @param[in] x  The state vector
  /// @param[out] f  The vector for the result f(x)
  /// @param[in] linear   Specifies that the SNES solver is in a linear (KSP) inner loop, so the operator should be linearised if possible
  PetscErrorCode snes_function(Vec x, Vec f, bool linear);

  /// Preconditioner. Called by PCapply
  /// via a C-style static function.
  ///
  /// @param[in] x  The vector to be operated on
  /// @param[out] f  The result of the operation
  PetscErrorCode precon(Vec x, Vec f);
 private:
  static constexpr int MAX_SUPPORTED_ORDER = 4; //Should this be #defined instead?

  int maxOrder; ///< Specify the maximum order of the scheme to use (1/2/3)

  BoutReal out_timestep; ///< The output timestep
  int nsteps; ///< Number of output steps
  BoutReal timestep; ///< The internal timestep
  int ninternal;     ///< Number of internal steps per output
  int mxstep; ///< Maximum number of internal steps between outputs

  //Adaptivity

  /// Use adaptive timestepping?
  bool adaptive; // Do we want to do an error check to enable adaptivity?
  int nadapt; ///< How often do we check the error
  int mxstepAdapt;  ///< Maximum no. consecutive times we try to reduce timestep
  BoutReal scaleCushUp; ///< Don't increase timestep if scale factor < 1.0+scaleCushUp
  BoutReal scaleCushDown; ///< Don't decrease timestep if scale factor > 1.0-scaleCushDown
  BoutReal adaptRtol; ///< Target relative error for adaptivity.
  BoutReal dtMin; ///< Minimum timestep we want to use
  BoutReal dtMax; ///< Maximum timestep we want to use
  BoutReal dtMinFatal; ///< If timestep wants to drop below this we abort. Set -ve to deactivate

  //Scheme coefficients
  std::vector<BoutReal> uFac, fFac, gFac;
  BoutReal dtImp;

  int nlocal, neq; ///< Number of variables on local processor and in total

  /// Take a full step at requested order
  ///
  /// @param[in] curtime  The current simulation time
  /// @param[in] dt       The time step to take
  /// @param[in] order    The order of accuracy
  void take_step(BoutReal curtime, BoutReal dt, int order=2);

  /// Setup a SNES object
  /// This includes creating, setting functions, options,
  /// and internal (KSP) solver, and Jacobian options
  /// including coloring.
  ///
  void constructSNES(SNES *snesIn);

  /// Shuffle state along one step
  void shuffleState();

  /// Populate the *Fac vectors and dtImp with appropriate coefficients for this order
  void calculateCoeffs(int order);

  // Working memory
  Array<BoutReal> u ; ///< System state at current time
  std::vector<Array<BoutReal>> uV; ///< The solution history
  std::vector<Array<BoutReal>> fV; ///< The non-stiff solution history
  std::vector<BoutReal> timesteps; ///< Timestep history
  Array<BoutReal> rhs;
  Array<BoutReal> err;

  // Implicit solver
  PetscErrorCode solve_implicit(BoutReal curtime, BoutReal gamma);
  BoutReal implicit_gamma;
  BoutReal implicit_curtime;
  int predictor;    ///< Predictor method
  PetscLib lib; ///< Handles initialising, finalising PETSc
  Vec      snes_f;  ///< Used by SNES to store function
  Vec      snes_x;  ///< Result of SNES
  SNES     snes;    ///< SNES context
  SNES     snesAlt; ///< Alternative SNES object for adaptive checks
  SNES     snesUse; ///< The snes object to use in solve stage. Allows easy switching.
  Mat      Jmf;     ///< Matrix-free Jacobian

  // Diagnostics
  bool diagnose;  ///< Output diagnostics every timestep
  bool verbose;  ///< Gives a more verbose output for each timestep
  int linear_fails;   ///< Number of linear (KSP) convergence failures
  int nonlinear_fails;  ///< Numbef of nonlinear (SNES) convergence failures

  bool have_constraints; ///< Are there any constraint variables?
  Array<BoutReal> is_dae; ///< If using constraints, 1 -> DAE, 0 -> AE

  MatFDColoring fdcoloring; ///< Matrix coloring context, used for finite difference Jacobian evaluation

  template< class Op >
  void loopVars(BoutReal *u);

  /// Save variables from BOUT++ fields into a
  /// pre-allocated array \p u
  void saveVars(BoutReal *u);

  /// Load variables from input vector u into BOUT++ fields
  void loadVars(BoutReal *u);

  /// Save time derivatives from ddt() fields into
  /// a preallocated array \p u.
  void saveDerivs(BoutReal *u);
};

#endif // __IMEXBDF2_SOLVER_H__

#endif // BOUT_HAS_PETSC
