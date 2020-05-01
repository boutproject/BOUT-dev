/**************************************************************************
 * 
 * Finds the steady-state solution of a set of equations
 * using PETSc for the SNES interface
 * 
 **************************************************************************
 * Copyright 2015 B.D.Dudson
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

#include "bout/build_config.hxx"

#if BOUT_HAS_PETSC

class SNESSolver;

#ifndef __SNES_SOLVER_H__
#define __SNES_SOLVER_H__

#include "mpi.h"

#include <bout_types.hxx>
#include <bout/solver.hxx>

#include <bout/petsclib.hxx>

#include <petsc.h>
#include <petscsnes.h>
// PETSc creates macros for MPI calls, which interfere with the MpiWrapper class
#undef MPI_Allreduce

namespace {
RegisterSolver<SNESSolver> registersolversnes("snes");
}

/// Uses PETSc's SNES interface to find a steady state solution to a
/// nonlinear ODE
class SNESSolver : public Solver {
 public:
  SNESSolver(Options *opt = nullptr) : Solver(opt) {}
  ~SNESSolver() {}
  
  int init(int nout, BoutReal tstep) override;
  
  int run() override;
  
  PetscErrorCode snes_function(Vec x, Vec f); ///< Nonlinear function
 private:
  int mxstep; ///< Maximum number of internal steps between outputs
  
  int nlocal; ///< Number of variables on local processor
  int neq; ///< Number of variables in total
  
  PetscLib lib;     ///< Handles initialising, finalising PETSc
  Vec      snes_f;  ///< Used by SNES to store function
  Vec      snes_x;  ///< Result of SNES
  SNES     snes;    ///< SNES context
  Mat      Jmf;     ///< Matrix-free Jacobian
  
};

#endif // __SNES_SOLVER_H__

#endif // BOUT_HAS_PETSC
