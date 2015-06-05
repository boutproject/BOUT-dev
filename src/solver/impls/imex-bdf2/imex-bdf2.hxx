/**************************************************************************
 * 2nd order IMEX-BDF scheme
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
  
  PetscErrorCode snes_function(Vec x, Vec f); // Nonlinear function
 private:
  int mxstep; // Maximum number of internal steps between outputs
  
  BoutReal out_timestep; // The output timestep
  int nsteps; // Number of output steps
  
  BoutReal timestep; // The internal timestep
  int ninternal;     // Number of internal steps per output
  
  int nlocal, neq; // Number of variables on local processor and in total
  
  // Startup step
  void startup(BoutReal curtime, BoutReal dt);
  
  // Take a full step
  void take_step(BoutReal curtime, BoutReal dt); 
  
  // Working memory
  BoutReal *u, *u_1, *u_2; // System state at n, n-1 and n-2
  BoutReal *f_1, *f_2;     // F(u_1) and F(u_2)
  BoutReal *rhs;
  
  // Implicit solver
  PetscErrorCode solve_implicit(BoutReal curtime, BoutReal gamma);
  BoutReal implicit_gamma;
  BoutReal implicit_curtime;
  PetscLib lib; // Handles initialising, finalising PETSc
  Vec      snes_f;  // Used by SNES to store function
  Vec      snes_x;  // Result of SNES
  SNES     snes;    // SNES context
  Mat      Jmf;     // Matrix-free Jacobian
};

#endif // __IMEXBDF2_SOLVER_H__

#else // BOUT_HAS_PETSC

#include "../emptysolver.hxx"
typedef EmptySolver IMEXBDF2;

#endif // BOUT_HAS_PETSC
