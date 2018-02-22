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

#ifndef BOUT_HAS_CVODE

#include "../emptysolver.hxx"
typedef EmptySolver CvodeSolver;
 
#else
class CvodeSolver;

#ifndef __SUNDIAL_SOLVER_H__
#define __SUNDIAL_SOLVER_H__

// NOTE: MPI must be included before SUNDIALS, otherwise complains
#include "mpi.h"

#include "bout_types.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

#include "bout/solver.hxx"

#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_bbdpre.h>
#include <nvector/nvector_parallel.h>

#include <vector>
using std::vector;

class CvodeSolver : public Solver {
  public:
    CvodeSolver(Options *opts = NULL);
    ~CvodeSolver();

    void setJacobian(Jacobian j) {jacfunc = j; }
    
    BoutReal getCurrentTimestep() { return hcur; }

    int init(int nout, BoutReal tstep) override;

    int run() override;
    BoutReal run(BoutReal tout);
    
    void resetInternalFields();

    // These functions used internally (but need to be public)
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
    void *cvode_mem;

    BoutReal pre_Wtime; // Time in preconditioner
    BoutReal pre_ncalls; // Number of calls to preconditioner
    
    void set_abstol_values(BoutReal* abstolvec_data, vector<BoutReal> &f2dtols, vector<BoutReal> &f3dtols);
    void loop_abstol_values_op(int jx, int jy, BoutReal* abstolvec_data, int &p, vector<BoutReal> &f2dtols, vector<BoutReal> &f3dtols, bool bndry);
};

#endif // __SUNDIAL_SOLVER_H__

#endif

