/**************************************************************************
 * Interface to ARKODE solver
 * NOTE: ARKode is currently in beta testing so use with cautious optimism
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

#ifndef BOUT_HAS_ARKODE

#include "../emptysolver.hxx"
typedef EmptySolver ArkodeSolver;
 
#else
class ArkodeSolver;

#ifndef __ARKODE_SOLVER_H__
#define __ARKODE_SOLVER_H__

// NOTE: MPI must be included before SUNDIALS, otherwise complains
#include "mpi.h"

#include "bout_types.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

#include "bout/solver.hxx"

#include <arkode/arkode_spgmr.h>
#include <arkode/arkode_bbdpre.h>
#include <nvector/nvector_parallel.h>

#include <vector>
using std::vector;

#include <bout/solverfactory.hxx>
namespace {
RegisterSolver<ArkodeSolver> registersolverarkode("arkode");
}

class ArkodeSolver : public Solver {
  public:
    ArkodeSolver(Options *opts = NULL);
    ~ArkodeSolver();

    void setJacobian(Jacobian j) {jacfunc = j; }
    
    BoutReal getCurrentTimestep() { return hcur; }

    int init(int nout, BoutReal tstep) override;

    int run() override;
    BoutReal run(BoutReal tout);

    // These functions used internally (but need to be public)
    void rhs_e(BoutReal t, BoutReal *udata, BoutReal *dudata);
    void rhs_i(BoutReal t, BoutReal *udata, BoutReal *dudata);
    void rhs(BoutReal t, BoutReal *udata, BoutReal *dudata);
    void pre(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal *udata, BoutReal *rvec, BoutReal *zvec);
    void jac(BoutReal t, BoutReal *ydata, BoutReal *vdata, BoutReal *Jvdata);
  private:
    int NOUT; // Number of outputs. Specified in init, needed in run
    BoutReal TIMESTEP; // Time between outputs
    BoutReal hcur; // Current internal timestep
  
    Jacobian jacfunc; // Jacobian - vector function
    bool diagnose; // Output additional diagnostics
  
    N_Vector uvec; // Values
    void *arkode_mem;

    BoutReal pre_Wtime; // Time in preconditioner
    BoutReal pre_ncalls; // Number of calls to preconditioner
    
    void set_abstol_values(BoutReal* abstolvec_data, vector<BoutReal> &f2dtols, vector<BoutReal> &f3dtols);
    void loop_abstol_values_op(int jx, int jy, BoutReal* abstolvec_data, int &p, vector<BoutReal> &f2dtols, vector<BoutReal> &f3dtols, bool bndry);
};

#endif // __ARKODE_SOLVER_H__

#endif

