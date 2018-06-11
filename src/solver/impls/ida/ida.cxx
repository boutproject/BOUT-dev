/**************************************************************************
 * Interface to SUNDIALS IDA
 * 
 * IdaSolver for DAE systems (so can handle constraints)
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

#include "ida.hxx"

#ifdef BOUT_HAS_IDA

#include <boutcomm.hxx>
#include <interpolation.hxx> // Cell interpolation
#include <msg_stack.hxx>
#include <boutexception.hxx>

#include <ida/ida.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_bbdpre.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <output.hxx>

#include "unused.hxx"

#define ZERO        RCONST(0.)
#define ONE         RCONST(1.0)

#ifndef IDAINT
typedef int IDAINT;
#endif

static int idares(BoutReal t, N_Vector u, N_Vector du, N_Vector rr, void *user_data);
static int ida_bbd_res(IDAINT Nlocal, BoutReal t, 
		       N_Vector u, N_Vector du, N_Vector rr, void *user_data);
static int ida_pre(BoutReal t, N_Vector yy, 	 
		   N_Vector yp, N_Vector rr, 	 
		   N_Vector rvec, N_Vector zvec, 	 
		   BoutReal cj, BoutReal delta, 
		   void *user_data, N_Vector tmp);

IdaSolver::IdaSolver(Options *opts) : Solver(opts) {
  has_constraints = true; ///< This solver has constraints
}

IdaSolver::~IdaSolver() { }

/**************************************************************************
 * Initialise
 **************************************************************************/

 int IdaSolver::init(int nout, BoutReal tstep) {

  TRACE("Initialising IDA solver");

  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;
  
  // Save nout and tstep for use in run
  NOUT = nout;
  TIMESTEP = tstep;
  
  output.write("Initialising IDA solver\n");

  // Calculate number of variables
  int n2d = f2d.size();
  int n3d = f3d.size();
  int local_N = getLocalN();

  // Get total problem size
  int neq;
  if(MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    output_error.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3d, n2d, neq, local_N);

  // Allocate memory
  
  if((uvec = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == NULL)
    throw BoutException("ERROR: SUNDIALS memory allocation failed\n");
  if((duvec = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == NULL)
    throw BoutException("ERROR: SUNDIALS memory allocation failed\n");
  if((id = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == NULL)
    throw BoutException("ERROR: SUNDIALS memory allocation failed\n");
  
  // Put the variables into uvec
  save_vars(NV_DATA_P(uvec));
  
  // Get the starting time derivative
  run_rhs(simtime);
  
  // Put the time-derivatives into duvec
  save_derivs(NV_DATA_P(duvec));
  
  // Set the equation type in id(Differential or Algebraic. This is optional)
  set_id(NV_DATA_P(id));
  
  /// Get options
  int MXSUB = mesh->xend - mesh->xstart + 1;

  BoutReal abstol, reltol;
  int maxl;
  int mudq, mldq;
  int mukeep, mlkeep;
  bool use_precon;
  bool correct_start;

  OPTION(options, mudq, n3d*(MXSUB+2));
  OPTION(options, mldq, n3d*(MXSUB+2));
  OPTION(options, mukeep, n3d);
  OPTION(options, mlkeep, n3d);
  options->get("ATOL", abstol, 1.0e-12);
  options->get("RTOL", reltol, 1.0e-5);
  OPTION(options, maxl, 6*n3d);
  OPTION(options, use_precon, false);
  OPTION(options, correct_start, true);
  int mxsteps; // Maximum number of steps to take between outputs
  options->get("mxstep", mxsteps, 500);

  // Call IDACreate and IDAMalloc to initialise

  if((idamem = IDACreate()) == NULL)
    throw BoutException("ERROR: IDACreate failed\n");
  
  if( IDASetUserData(idamem, this) < 0 ) // For callbacks, need pointer to solver object
    throw BoutException("ERROR: IDASetUserData failed\n");

  if( IDASetId(idamem, id) < 0)
    throw BoutException("ERROR: IDASetID failed\n");

  if( IDAInit(idamem, idares, simtime, uvec, duvec) < 0 )
    throw BoutException("ERROR: IDAInit failed\n");
  
  if( IDASStolerances(idamem, reltol, abstol) < 0 )
    throw BoutException("ERROR: IDASStolerances failed\n");

  IDASetMaxNumSteps(idamem, mxsteps);

  // Call IDASpgmr to specify the IDA linear solver IDASPGMR
  if( IDASpgmr(idamem, maxl) )
    throw BoutException("ERROR: IDASpgmr failed\n");

  if(use_precon) {
    if(!have_user_precon()) {
      output.write("\tUsing BBD preconditioner\n");
      if( IDABBDPrecInit(idamem, local_N, mudq, mldq, mukeep, mlkeep, 
			 ZERO, ida_bbd_res, NULL) )
	throw BoutException("ERROR: IDABBDPrecInit failed\n");
    }else {
      output.write("\tUsing user-supplied preconditioner\n");
      if( IDASpilsSetPreconditioner(idamem, NULL, ida_pre) )
	throw BoutException("ERROR: IDASpilsSetPreconditioner failed\n");
    }
  }

  // Call IDACalcIC (with default options) to correct the initial values
  if(correct_start) {
    if( IDACalcIC(idamem, IDA_YA_YDP_INIT, 1e-6) )
      throw BoutException("ERROR: IDACalcIC failed\n");
  }

  return 0;
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int IdaSolver::run() {
  TRACE("IDA IdaSolver::run()");
  
  if(!initialised)
    throw BoutException("IdaSolver not initialised\n");

  for(int i=0;i<NOUT;i++) {
    
    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP);
    iteration++;

    /// Check if the run succeeded
    if(simtime < 0.0) {
      // Step failed
      throw BoutException("SUNDIALS IDA timestep failed\n");
    }
    
    /// Call the monitor function
    
    if(call_monitors(simtime, i, NOUT)) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

BoutReal IdaSolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run(%e)", tout);

  if(!initialised)
    throw BoutException("ERROR: Running IDA solver without initialisation\n");

  pre_Wtime = 0.0;
  pre_ncalls = 0.0;

  int flag = IDASolve(idamem, tout, &simtime, uvec, duvec, IDA_NORMAL);

  // Copy variables
  load_vars(NV_DATA_P(uvec));

  // Call rhs function to get extra variables at this time
  run_rhs(simtime);
  
  if(flag < 0) {
    output_error.write("ERROR IDA solve failed at t = %e, flag = %d\n", simtime, flag);
    return -1.0;
  }

  return simtime;
}

/**************************************************************************
 * Residual function F(t, u, du)
 **************************************************************************/

void IdaSolver::res(BoutReal t, BoutReal *udata, BoutReal *dudata, BoutReal *rdata)
{
  TRACE("Running RHS: IdaSolver::res(%e)", t);
  
  // Load state from udata
  load_vars(udata);
  
  // Call RHS function
  run_rhs(t);
  
  // Save derivatives to rdata (residual)
  save_derivs(rdata);
  
  // If a differential equation, subtract dudata
  int N = NV_LOCLENGTH_P(id);
  BoutReal *idd = NV_DATA_P(id);
  for(int i=0;i<N;i++) {
    if(idd[i] > 0.5) // 1 -> differential, 0 -> algebraic
      rdata[i] -= dudata[i];
  }
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void IdaSolver::pre(BoutReal t, BoutReal cj, BoutReal delta, BoutReal *udata, BoutReal *rvec, BoutReal *zvec) {
  TRACE("Running preconditioner: IdaSolver::pre(%e)", t);

  BoutReal tstart = MPI_Wtime();

  int N = NV_LOCLENGTH_P(id);
  
  if(!have_user_precon()) {
    // Identity (but should never happen)
    for(int i=0;i<N;i++)
      zvec[i] = rvec[i];
    return;
  }

  // Load state from udata (as with res function)
  load_vars(udata);

  // Load vector to be inverted into F_vars
  load_derivs(rvec);
  
  run_precon(t, cj, delta);

  // Save the solution from F_vars
  save_derivs(zvec);

  pre_Wtime += MPI_Wtime() - tstart;
  pre_ncalls++;
}

/**************************************************************************
 * IDA res function
 **************************************************************************/

static int idares(BoutReal t, 
                  N_Vector u, N_Vector du, N_Vector rr, 
                  void *user_data)
{
  BoutReal *udata = NV_DATA_P(u);
  BoutReal *dudata = NV_DATA_P(du);
  BoutReal *rdata = NV_DATA_P(rr);

  IdaSolver *s = static_cast<IdaSolver *>(user_data);

  // Calculate residuals
  s->res(t, udata, dudata, rdata);

  return 0;
}

/// Residual function for BBD preconditioner
static int ida_bbd_res(IDAINT UNUSED(Nlocal), BoutReal t, N_Vector u, N_Vector du,
                       N_Vector rr, void *user_data) {
  return idares(t, u, du, rr, user_data);
}

// Preconditioner function
static int ida_pre(BoutReal t, N_Vector yy, N_Vector UNUSED(yp), N_Vector UNUSED(rr),
                   N_Vector rvec, N_Vector zvec, BoutReal cj, BoutReal delta,
                   void *user_data, N_Vector UNUSED(tmp)) {
  BoutReal *udata = NV_DATA_P(yy);
  BoutReal *rdata = NV_DATA_P(rvec);
  BoutReal *zdata = NV_DATA_P(zvec);

  IdaSolver *s = static_cast<IdaSolver *>(user_data);

  // Calculate residuals
  s->pre(t, cj, delta, udata, rdata, zdata);

  return 0;
}

#endif
