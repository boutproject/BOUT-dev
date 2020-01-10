/**************************************************************************
 * Interface to PVODE solver
 * NOTE: This class needs tidying, generalising to use FieldData interface
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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

#ifdef BOUT_HAS_PVODE

class PvodeSolver;

#ifndef __PVODE_SOLVER_H__
#define __PVODE_SOLVER_H__

#include <bout/solver.hxx>
#include <bout_types.hxx>

#include <pvode/nvector.h>
#include <pvode/cvode.h>     // main CVODE header file
#include <pvode/pvbbdpre.h>  // band preconditioner function prototypes

namespace {
RegisterSolver<PvodeSolver> registersolverpvode("pvode");
}

class PvodeSolver : public Solver {
 public:
  PvodeSolver(Options *opts);
  ~PvodeSolver();
  
  BoutReal getCurrentTimestep() override { return hcur; }
  
  int init(int nout, BoutReal tstep) override;
  
  int run() override;
  BoutReal run(BoutReal tout);

  // These functions used internally (but need to be public)
  void rhs(int N, BoutReal t, BoutReal *udata, BoutReal *dudata);
  void gloc(int N, BoutReal t, BoutReal *udata, BoutReal *dudata);

 private:
  int NOUT; // Number of outputs. Specified in init, needed in run
  BoutReal TIMESTEP; // Time between outputs
  BoutReal hcur; // Current internal timestep

  pvode::N_Vector u;
  pvode::machEnvType machEnv;
  void *cvode_mem;

  BoutReal abstol, reltol; // addresses passed in init must be preserved
  pvode::PVBBDData pdata;
};

#endif // __PVODE_SOLVER_H__

#endif
