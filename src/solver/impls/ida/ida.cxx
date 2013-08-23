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

#include <ida/ida.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_bbdpre.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <output.hxx>

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

IdaSolver::IdaSolver() : Solver()
{
  has_constraints = true; ///< This solver has constraints
  
  prefunc = NULL;
}

IdaSolver::~IdaSolver()
{
  if(initialised) {
    // Free IDA memory
    
    
    
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

 int IdaSolver::init(bool restarting, int nout, BoutReal tstep) {

  int msg_point = msg_stack.push("Initialising IDA solver");

  /// Call the generic initialisation first
  if(Solver::init(restarting, nout, tstep))
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
    output.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3d, n2d, neq, local_N);

  // Allocate memory
  
  if((uvec = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == NULL)
    bout_error("ERROR: SUNDIALS memory allocation failed\n");
  if((duvec = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == NULL)
    bout_error("ERROR: SUNDIALS memory allocation failed\n");
  if((id = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == NULL)
    bout_error("ERROR: SUNDIALS memory allocation failed\n");
  
  // Put the variables into uvec
  if(save_vars(NV_DATA_P(uvec)))
    bout_error("\tERROR: Initial variable value not set\n");
  
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
  Options *options = Options::getRoot();
  options = options->getSection("solver");
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
    bout_error("ERROR: IDACreate failed\n");
  
  if( IDASetUserData(idamem, this) < 0 ) // For callbacks, need pointer to solver object
    bout_error("ERROR: IDASetUserData failed\n");

  if( IDASetId(idamem, id) < 0)
    bout_error("ERROR: IDASetID failed\n");

  if( IDAInit(idamem, idares, simtime, uvec, duvec) < 0 )
    bout_error("ERROR: IDAInit failed\n");
  
  if( IDASStolerances(idamem, reltol, abstol) < 0 )
    bout_error("ERROR: IDASStolerances failed\n");

  IDASetMaxNumSteps(idamem, mxsteps);

  // Call IDASpgmr to specify the IDA linear solver IDASPGMR
  if( IDASpgmr(idamem, maxl) )
    bout_error("ERROR: IDASpgmr failed\n");

  if(use_precon) {
    if(prefunc == NULL) {
      output.write("\tUsing BBD preconditioner\n");
      if( IDABBDPrecInit(idamem, local_N, mudq, mldq, mukeep, mlkeep, 
			 ZERO, ida_bbd_res, NULL) )
	bout_error("ERROR: IDABBDPrecInit failed\n");
    }else {
      output.write("\tUsing user-supplied preconditioner\n");
      if( IDASpilsSetPreconditioner(idamem, NULL, ida_pre) )
	bout_error("ERROR: IDASpilsSetPreconditioner failed\n");
    }
  }

  // Call IDACalcIC (with default options) to correct the initial values
  if(correct_start) {
    if( IDACalcIC(idamem, IDA_YA_YDP_INIT, 1e-6) )
      bout_error("ERROR: IDACalcIC failed\n");
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return(0);
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int IdaSolver::run() {
#ifdef CHECK
  int msg_point = msg_stack.push("IDA IdaSolver::run()");
#endif
  
  if(!initialised)
    bout_error("IdaSolver not initialised\n");

  for(int i=0;i<NOUT;i++) {
    
    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP);
    iteration++;

    /// Check if the run succeeded
    if(simtime < 0.0) {
      // Step failed
      output.write("Timestep failed. Aborting\n");

      // Write restart to a different file
      restart.write("%s/BOUT.failed.%s", restartdir.c_str(), restartext.c_str());

      bout_error("SUNDIALS IDA timestep failed\n");
    }
    
    /// Write the restart file
    restart.write();
    
    if((archive_restart > 0) && (iteration % archive_restart == 0)) {
      restart.write("%s/BOUT.restart_%04d.%s", restartdir.c_str(), iteration, restartext.c_str());
    }
    
    /// Call the monitor function
    
    if(call_monitors(simtime, i, NOUT)) {
      // User signalled to quit
      
      // Write restart to a different file
      restart.write("%s/BOUT.final.%s", restartdir.c_str(), restartext.c_str());
      
      output.write("Monitor signalled to quit. Returning\n");
      break;
    }
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}

BoutReal IdaSolver::run(BoutReal tout) {
  if(!initialised)
    bout_error("ERROR: Running IDA solver without initialisation\n");

#ifdef CHECK
  int msg_point = msg_stack.push("Running solver: solver::run(%e)", tout);
#endif
  
  rhs_ncalls = 0;

  pre_Wtime = 0.0;
  pre_ncalls = 0.0;

  int flag = IDASolve(idamem, tout, &simtime, uvec, duvec, IDA_NORMAL);

  // Copy variables
  load_vars(NV_DATA_P(uvec));

  // Call rhs function to get extra variables at this time
  run_rhs(simtime);
  
  if(flag < 0) {
    output.write("ERROR IDA solve failed at t = %e, flag = %d\n", simtime, flag);
    return -1.0;
  }
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return simtime;
}

/**************************************************************************
 * Residual function F(t, u, du)
 **************************************************************************/

void IdaSolver::res(BoutReal t, BoutReal *udata, BoutReal *dudata, BoutReal *rdata)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: IdaSolver::res(%e)", t);
#endif
  
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

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void IdaSolver::pre(BoutReal t, BoutReal cj, BoutReal delta, BoutReal *udata, BoutReal *rvec, BoutReal *zvec) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running preconditioner: IdaSolver::pre(%e)", t);
#endif

  BoutReal tstart = MPI_Wtime();

  int N = NV_LOCLENGTH_P(id);
  
  if(prefunc == NULL) {
    // Identity (but should never happen)
    for(int i=0;i<N;i++)
      zvec[i] = rvec[i];
    return;
  }

  // Load state from udata (as with res function)
  load_vars(udata);

  // Load vector to be inverted into F_vars
  load_derivs(rvec);
  
  (*prefunc)(t, cj, delta);

  // Save the solution from F_vars
  save_derivs(zvec);

  pre_Wtime += MPI_Wtime() - tstart;
  pre_ncalls++;

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * PRIVATE FUNCTIONS
 **************************************************************************/

/// Perform an operation at a given (jx,jy) location, moving data between BOUT++ and CVODE
void IdaSolver::loop_vars_op(int jx, int jy, BoutReal *udata, int &p, SOLVER_VAR_OP op)
{
  BoutReal **d2d, ***d3d;
  int i;
  int jz;
 
  int n2d = f2d.size();
  int n3d = f3d.size();

  switch(op) {
  case LOAD_VARS: {
    /// Load variables from IDA into BOUT++
    
    // Loop over 2D variables
    for(i=0;i<n2d;i++) {
      d2d = f2d[i].var->getData(); // Get pointer to data
      d2d[jx][jy] = udata[p];
      p++;
    }
    
    for (jz=0; jz < mesh->ngz-1; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	d3d = f3d[i].var->getData(); // Get pointer to data
	d3d[jx][jy][jz] = udata[p];
	p++;
      }  
    }
    break;
  }
  case LOAD_DERIVS: {
    /// Load derivatives from IDA into BOUT++
    /// Used for preconditioner
    
    // Loop over 2D variables
    for(i=0;i<n2d;i++) {
      d2d = f2d[i].F_var->getData(); // Get pointer to data
      d2d[jx][jy] = udata[p];
      p++;
    }
    
    for (jz=0; jz < mesh->ngz-1; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	d3d = f3d[i].F_var->getData(); // Get pointer to data
	d3d[jx][jy][jz] = udata[p];
	p++;
      }  
    }
    
    break;
  }
  case SET_ID: {
    /// Set the type of equation (Differential or Algebraic)
    
    // Loop over 2D variables
    for(i=0;i<n2d;i++) {
      if(f2d[i].constraint) {
	udata[p] = ZERO;
      }else {
	udata[p] = ONE;
      }
      p++;
    }
    
    for (jz=0; jz < mesh->ngz-1; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	if(f3d[i].constraint) {
	  udata[p] = ZERO;
	}else {
	  udata[p] = ONE;
	}
	p++;
      }
    }
    
    break;
  }
  case SAVE_VARS: {
    /// Save variables from BOUT++ into IDA (only used at start of simulation)
    
    // Loop over 2D variables
    for(i=0;i<n2d;i++) {
      d2d = f2d[i].var->getData(); // Get pointer to data
      udata[p] = d2d[jx][jy];
      p++;
    }
    
    for (jz=0; jz < mesh->ngz-1; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	d3d = f3d[i].var->getData(); // Get pointer to data
	udata[p] = d3d[jx][jy][jz];
	p++;
      }  
    }
    break;
  }
    /// Save time-derivatives from BOUT++ into CVODE (returning RHS result)
  case SAVE_DERIVS: {
    
    // Loop over 2D variables
    for(i=0;i<n2d;i++) {
      d2d = f2d[i].F_var->getData(); // Get pointer to data
      udata[p] = d2d[jx][jy];
      p++;
    }
    
    for (jz=0; jz < mesh->ngz-1; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	d3d = f3d[i].F_var->getData(); // Get pointer to data
	udata[p] = d3d[jx][jy][jz];
	p++;
      }  
    }
    break;
  }
  }
}

/// Loop over variables and domain. Used for all data operations for consistency
void IdaSolver::loop_vars(BoutReal *udata, SOLVER_VAR_OP op)
{
  int jx, jy;
  int p = 0; // Counter for location in udata array

  int MYSUB = mesh->yend - mesh->ystart + 1;

  // Inner X boundary
  if(mesh->firstX()) {
    for(jx=0;jx<mesh->xstart;jx++)
      for(jy=0;jy<MYSUB;jy++)
	loop_vars_op(jx, jy+mesh->ystart, udata, p, op);
  }

  // Lower Y boundary region
  for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); xi++) {
    for(jy=0;jy<mesh->ystart;jy++)
      loop_vars_op(*xi, jy, udata, p, op);
  }
  
  // Upper Y boundary condition
  for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); xi++) {
    for(jy=mesh->yend+1;jy<mesh->ngy;jy++)
      loop_vars_op(*xi, jy, udata, p, op);
  }

  // Bulk of points
  for (jx=mesh->xstart; jx <= mesh->xend; jx++)
    for (jy=mesh->ystart; jy <= mesh->yend; jy++)
      loop_vars_op(jx, jy, udata, p, op);

  // Outer X boundary
  if(mesh->lastX()) {
    for(jx=mesh->xend+1;jx<mesh->ngx;jx++)
      for(jy=mesh->ystart;jy<=mesh->yend;jy++)
	loop_vars_op(jx, jy, udata, p, op);
  }
}

void IdaSolver::load_vars(BoutReal *udata)
{
  unsigned int i;
  
  // Make sure data is allocated
  for(i=0;i<f2d.size();i++)
    f2d[i].var->allocate();
  for(i=0;i<f3d.size();i++) {
    f3d[i].var->allocate();
    f3d[i].var->setLocation(f3d[i].location);
  }

  loop_vars(udata, LOAD_VARS);

  // Mark each vector as either co- or contra-variant

  for(i=0;i<v2d.size();i++)
    v2d[i].var->covariant = v2d[i].covariant;
  for(i=0;i<v3d.size();i++)
    v3d[i].var->covariant = v3d[i].covariant;
}

void IdaSolver::load_derivs(BoutReal *udata)
{
  unsigned int i;
  
  // Make sure data is allocated
  for(i=0;i<f2d.size();i++)
    f2d[i].F_var->allocate();
  for(i=0;i<f3d.size();i++) {
    f3d[i].F_var->allocate();
    f3d[i].F_var->setLocation(f3d[i].location);
  }

  loop_vars(udata, LOAD_DERIVS);

  // Mark each vector as either co- or contra-variant

  for(i=0;i<v2d.size();i++)
    v2d[i].F_var->covariant = v2d[i].covariant;
  for(i=0;i<v3d.size();i++)
    v3d[i].F_var->covariant = v3d[i].covariant;
}

void IdaSolver::set_id(BoutReal *udata)
{
  loop_vars(udata, SET_ID);
}

// This function only called during initialisation
int IdaSolver::save_vars(BoutReal *udata)
{
  unsigned int i;

  for(i=0;i<f2d.size();i++)
    if(f2d[i].var->getData() == (BoutReal**) NULL)
      return(1);

  for(i=0;i<f3d.size();i++)
    if(f3d[i].var->getData() == (BoutReal***) NULL)
      return(1);
  
  // Make sure vectors in correct basis
  for(i=0;i<v2d.size();i++) {
    if(v2d[i].covariant) {
      v2d[i].var->toCovariant();
    }else
      v2d[i].var->toContravariant();
  }
  for(i=0;i<v3d.size();i++) {
    if(v3d[i].covariant) {
      v3d[i].var->toCovariant();
    }else
      v3d[i].var->toContravariant();
  }

  loop_vars(udata, SAVE_VARS);

  return(0);
}

void IdaSolver::save_derivs(BoutReal *dudata)
{
  unsigned int i;

  // Make sure vectors in correct basis
  for(i=0;i<v2d.size();i++) {
    if(v2d[i].covariant) {
      v2d[i].F_var->toCovariant();
    }else
      v2d[i].F_var->toContravariant();
  }
  for(i=0;i<v3d.size();i++) {
    if(v3d[i].covariant) {
      v3d[i].F_var->toCovariant();
    }else
      v3d[i].F_var->toContravariant();
  }

  // Make sure 3D fields are at the correct cell location
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    if((*it).location != ((*it).F_var)->getLocation()) {
      //output.write("SOLVER: Interpolating\n");
      *((*it).F_var) = interp_to(*((*it).F_var), (*it).location);
    }
  }

  loop_vars(dudata, SAVE_DERIVS);
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
  
  IdaSolver *s = (IdaSolver*) user_data;

  // Calculate residuals
  s->res(t, udata, dudata, rdata);

  return 0;
}

/// Residual function for BBD preconditioner
static int ida_bbd_res(IDAINT Nlocal, BoutReal t, 
		       N_Vector u, N_Vector du, N_Vector rr, 
		       void *user_data)
{
  return idares(t, u, du, rr, user_data);
}

// Preconditioner function
static int ida_pre(BoutReal t, N_Vector yy, 	 
		   N_Vector yp, N_Vector rr, 	 
		   N_Vector rvec, N_Vector zvec, 	 
		   BoutReal cj, BoutReal delta, 
		   void *user_data, N_Vector tmp)
{
  BoutReal *udata = NV_DATA_P(yy);
  BoutReal *rdata = NV_DATA_P(rvec);
  BoutReal *zdata = NV_DATA_P(zvec);
  
  IdaSolver *s = (IdaSolver*) user_data;

  // Calculate residuals
  s->pre(t, cj, delta, udata, rdata, zdata);

  return 0;
}

#endif
