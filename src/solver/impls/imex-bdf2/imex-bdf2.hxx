/**************************************************************************
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

class IMEXBDF2 : public Solver {
 public:
  IMEXBDF2(Options *opt = NULL);
  ~IMEXBDF2();
  
  BoutReal getCurrentTimestep() {return timestep; }
  
  int init(bool restarting, int nout, BoutReal tstep);
  
  int run();
  
  PetscErrorCode snes_function(Vec x, Vec f, bool linear); // Nonlinear function
  PetscErrorCode precon(Vec x, Vec f); // Preconditioner
 private:
  static const int MAX_SUPPORTED_ORDER = 4; //Should this be #defined instead?
  
  int maxOrder; //Specify the maximum order of the scheme to use (1/2/3)

  BoutReal out_timestep; // The output timestep
  int nsteps; // Number of output steps
  BoutReal timestep; // The internal timestep
  int ninternal;     // Number of internal steps per output
  int mxstep; // Maximum number of internal steps between outputs

  //Adaptivity
  bool adaptive; //Do we want to do an error check to enable adaptivity?
  int nadapt; //How often do we check the error
  int mxstepAdapt;  //Maximum no. consecutive times we try to reduce timestep
  BoutReal scaleCushUp; //Don't increase timestep if scale factor < 1.0+scaleCushUp
  BoutReal scaleCushDown; //Don't decrease timestep if scale factor > 1.0-scaleCushDown
  BoutReal adaptRtol; //Target relative error for adaptivity.
  BoutReal dtMin; //Minimum timestep we want to use
  BoutReal dtMax; //Maximum timestep we want to use
  BoutReal dtMinFatal; //If timestep wants to drop below this we abort. Set -ve to deactivate
  
  //Scheme coefficients
  vector<BoutReal> uFac, fFac, gFac;
  BoutReal dtImp;

  int nlocal, neq; // Number of variables on local processor and in total
  
  // Take a full step at requested order
  void take_step(BoutReal curtime, BoutReal dt, int order=2); 

  void constructSNES(SNES *snesIn); //Setup a SNES object

  //Shuffle state along one step
  void shuffleState();

  //Populate the *Fac vectors and dtImp with appropriate coefficients for this order
  void calculateCoeffs(int order);

  // Working memory
  BoutReal *u ; // System state at current time 
  vector<BoutReal*> uV; //The solution history
  vector<BoutReal*> fV; //The non-stiff solution history
  //vector<BoutReal*> gV; //The stiff solution history
  vector<BoutReal> timesteps; //Timestep history
  BoutReal *rhs;
  BoutReal *err;

  // Implicit solver
  PetscErrorCode solve_implicit(BoutReal curtime, BoutReal gamma);
  BoutReal implicit_gamma;
  BoutReal implicit_curtime;
  int predictor;    // Predictor method
  PetscLib lib; // Handles initialising, finalising PETSc
  Vec      snes_f;  // Used by SNES to store function
  Vec      snes_x;  // Result of SNES
  SNES     snes;    // SNES context
  SNES     snesAlt; // Alternative SNES object for adaptive checks
  SNES     snesUse; // The snes object to use in solve stage. Allows easy switching.
  Mat      Jmf;     // Matrix-free Jacobian

  bool have_constraints; // Are there any constraint variables?
  BoutReal *is_dae; // 1 -> DAE, 0 -> AE
  
  MatFDColoring fdcoloring;
  
  template< class Op >
  void loopVars(BoutReal *u);
  
  void saveVars(BoutReal *u);
  void loadVars(BoutReal *u);
  void saveDerivs(BoutReal *u);
};

#endif // __IMEXBDF2_SOLVER_H__

#else // BOUT_HAS_PETSC

#include "../emptysolver.hxx"
typedef EmptySolver IMEXBDF2;

#endif // BOUT_HAS_PETSC
