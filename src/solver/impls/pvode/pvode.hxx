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

#include "bout/build_defines.hxx"

#if BOUT_HAS_PVODE

class PvodeSolver;

#ifndef BOUT_PVODE_SOLVER_H
#define BOUT_PVODE_SOLVER_H

#include <bout/bout_types.hxx>
#include <bout/solver.hxx>

#include <pvode/cvode.h> // main CVODE header file
#include <pvode/nvector.h>
#include <pvode/pvbbdpre.h> // band preconditioner function prototypes

namespace {
RegisterSolver<PvodeSolver> registersolverpvode("pvode");
}

class PvodeSolver : public Solver {
public:
  explicit PvodeSolver(Options* opts = nullptr);
  ~PvodeSolver();

  BoutReal getCurrentTimestep() override { return hcur; }

  int init() override;
  int run() override;
  BoutReal run(BoutReal tout);

  // These functions used internally (but need to be public)
  void rhs(int N, BoutReal t, BoutReal* udata, BoutReal* dudata);
  void gloc(int N, BoutReal t, BoutReal* udata, BoutReal* dudata);

private:
  /// Current internal timestep
  BoutReal hcur;
  /// Use preconditioner
  bool use_precon;
  /// Maximum Krylov dimension
  int precon_dimens;
  /// Preconditioner tolerance
  BoutReal precon_tol;
  /// Maximum number of steps
  int pvode_mxstep;

  pvode::N_Vector u{nullptr};
  pvode::machEnvType machEnv{nullptr};
  void* cvode_mem{nullptr};

  BoutReal abstol, reltol; // addresses passed in init must be preserved
  pvode::PVBBDData pdata{nullptr};

  bool pvode_initialised = false;
};

#endif // BOUT_PVODE_SOLVER_H

#endif
