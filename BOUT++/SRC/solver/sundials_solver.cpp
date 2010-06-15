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

#include "sundials_solver.h"

#include "globals.h"
#include "boundary.h"
#include "interpolation.h" // Cell interpolation
#include "communicator.h"

#include <cvode/cvode.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#define ZERO        RCONST(0.)
#define ONE         RCONST(1.0)

static int cvode_rhs(real t, N_Vector u, N_Vector du, void *user_data);
static int cvode_bbd_rhs(int Nlocal, real t, N_Vector u, N_Vector du, 
			 void *user_data);

static int cvode_pre(real t, N_Vector yy, N_Vector yp,
		     N_Vector rvec, N_Vector zvec,
		     real gamma, real delta, int lr,
		     void *user_data, N_Vector tmp);

static int cvode_jac(N_Vector v, N_Vector Jv,
		     realtype t, N_Vector y, N_Vector fy,
		     void *user_data, N_Vector tmp);

Solver::Solver() : GenericSolver()
{
  has_constraints = false; ///< This solver doesn't have constraints
  
  prefunc = NULL;
  jacfunc = NULL;
}

Solver::~Solver()
{
  
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int Solver::init(rhsfunc f, int argc, char **argv, bool restarting, int nout, real tstep)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Initialising CVODE solver");
#endif

  /// Call the generic initialisation first
  if(GenericSolver::init(f, argc, argv, restarting, nout, tstep))
    return 1;

  // Save nout and tstep for use in run
  NOUT = nout;
  TIMESTEP = tstep;

  output.write("Initialising SUNDIALS' CVODE solver\n");

  // Set the rhs solver function
  func = f;

  // Calculate number of variables (in generic_solver)
  int local_N = getLocalN();
  
  // Get total problem size
  int neq;
  if(MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD)) {
    output.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3Dvars(), n2Dvars(), neq, local_N);

  // Allocate memory
  
  if((uvec = N_VNew_Parallel(MPI_COMM_WORLD, local_N, neq)) == NULL)
    bout_error("ERROR: SUNDIALS memory allocation failed\n");
  
  // Put the variables into uvec
  if(save_vars(NV_DATA_P(uvec)))
    bout_error("\tERROR: Initial variable value not set\n");
  
  /// Get options

  real abstol, reltol;
  int maxl;
  int mudq, mldq;
  int mukeep, mlkeep;
  bool use_precon, use_jacobian;
  real max_timestep;
  bool adams_moulton, func_iter; // Time-integration method
  options.setSection("solver");
  options.get("mudq", mudq, n3Dvars()*(MXSUB+2));
  options.get("mldq", mldq, n3Dvars()*(MXSUB+2));
  options.get("mukeep", mukeep, n3Dvars()+n2Dvars());
  options.get("mlkeep", mlkeep, n3Dvars()+n2Dvars());
  options.get("ATOL", abstol, 1.0e-12);
  options.get("RTOL", reltol, 1.0e-5);
  options.get("maxl", maxl, 5);
  options.get("use_precon", use_precon, false);
  options.get("use_jacobian", use_jacobian, false);
  options.get("max_timestep", max_timestep, -1.);
  
  int mxsteps; // Maximum number of steps to take between outputs
  options.get("pvode_mxstep", mxsteps, 500);
  
  options.get("adams_moulton", adams_moulton, false);
  
  int lmm = CV_BDF;
  if(adams_moulton) {
    // By default use functional iteration for Adams-Moulton
    lmm = CV_ADAMS;
    output.write("\tUsing Adams-Moulton implicit multistep method\n");
    options.get("func_iter", func_iter, true); 
  }else {
    output.write("\tUsing BDF method\n");
    // Use Newton iteration for BDF
    options.get("func_iter", func_iter, false); 
  }
  
  int iter = CV_NEWTON;
  if(func_iter)
    iter = CV_FUNCTIONAL;

  // Call CVodeCreate
  if((cvode_mem = CVodeCreate(lmm, iter)) == NULL)
    bout_error("ERROR: CVodeCreate failed\n");
  
  if( CVodeSetUserData(cvode_mem, this) < 0 ) // For callbacks, need pointer to solver object
    bout_error("ERROR: CVodeSetUserData failed\n");

  if( CVodeInit(cvode_mem, cvode_rhs, simtime, uvec) < 0 )
    bout_error("ERROR: CVodeInit failed\n");
  
  if( CVodeSStolerances(cvode_mem, reltol, abstol) < 0 )
    bout_error("ERROR: CVodeSStolerances failed\n");

  CVodeSetMaxNumSteps(cvode_mem, mxsteps);

  if(max_timestep > 0.0) {
    // Setting a maximum timestep
    CVodeSetMaxStep(cvode_mem, max_timestep);
  }
  
  /// Newton method can include Preconditioners and Jacobian function
  if(!func_iter) {
    output.write("\tUsing Newton iteration\n");
    /// Set Preconditioner
    if(use_precon) {
      
      if( CVSpgmr(cvode_mem, PREC_LEFT, maxl) != CVSPILS_SUCCESS )
	bout_error("ERROR: CVSpgmr failed\n");
      
      if(prefunc == NULL) {
	output.write("\tUsing BBD preconditioner\n");
	
	if( CVBBDPrecInit(cvode_mem, local_N, mudq, mldq, 
			  mukeep, mlkeep, ZERO, cvode_bbd_rhs, NULL) )
	  bout_error("ERROR: CVBBDPrecInit failed\n");
	
      }else {
	output.write("\tUsing user-supplied preconditioner\n");
	
	if( CVSpilsSetPreconditioner(cvode_mem, NULL, cvode_pre) )
	  bout_error("ERROR: CVSpilsSetPreconditioner failed\n");
      }
    }else {
      // Not using preconditioning
      
      output.write("\tNo preconditioning\n");
      
      if( CVSpgmr(cvode_mem, PREC_NONE, maxl) != CVSPILS_SUCCESS )
	bout_error("ERROR: CVSpgmr failed\n");
    }
    
    /// Set Jacobian-vector multiplication function
    
    if((use_jacobian) && (jacfunc != NULL)) {
      output.write("\tUsing user-supplied Jacobian function\n");
      
      if( CVSpilsSetJacTimesVecFn(cvode_mem, cvode_jac) != CVSPILS_SUCCESS )
	bout_error("ERROR: CVSpilsSetJacTimesVecFn failed\n");
    
    }else 
      output.write("\tUsing difference quotient approximation for Jacobian\n");
  }else {
    output.write("\tUsing Functional iteration\n");
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return(0);
}


/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int Solver::run(MonitorFunc monitor)
{
#ifdef CHECK
  int msg_point = msg_stack.push("CVODE Solver::run()");
#endif
  
  if(!initialised)
    bout_error("Solver not initialised\n");

  for(int i=0;i<NOUT;i++) {
    
    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP, rhs_ncalls, rhs_wtime);
    iteration++;

    /// Check if the run succeeded
    if(simtime < 0.0) {
      // Step failed
      output.write("Timestep failed. Aborting\n");

      // Write restart to a different file
      restart.write("%s/BOUT.final.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());

      bout_error("SUNDIALS timestep failed\n");
    }

    /// Write the restart file
    restart.write("%s/BOUT.restart.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());
    
    if((archive_restart > 0) && (iteration % archive_restart == 0)) {
      restart.write("%s/BOUT.restart_%04d.%d.%s", restartdir.c_str(), iteration, MYPE, restartext.c_str());
    }
    
    /// Call the monitor function
    
    if(monitor(simtime, i, NOUT)) {
      // User signalled to quit
      
      // Write restart to a different file
      restart.write("%s/BOUT.final.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());
      
      output.write("Monitor signalled to quit. Returning\n");
      break;
    }
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}

real Solver::run(real tout, int &ncalls, real &rhstime)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Running solver: solver::run(%e)", tout);
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  rhs_wtime = 0.0;
  rhs_ncalls = 0;

  pre_Wtime = 0.0;
  pre_ncalls = 0.0;

  int flag = CVode(cvode_mem, tout, uvec, &simtime, CV_NORMAL);
  
  ncalls = rhs_ncalls;
  rhstime = rhs_wtime;

  // Copy variables
  load_vars(NV_DATA_P(uvec));

  // Call rhs function to get extra variables at this time
  real tstart = MPI_Wtime();
  (*func)(simtime);
  rhs_wtime += MPI_Wtime() - tstart;
  rhs_ncalls++;
  
  if(flag < 0) {
    output.write("ERROR CVODE solve failed at t = %e, flag = %d\n", simtime, flag);
    return -1.0;
  }
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return simtime;
}

/**************************************************************************
 * RHS function du = F(t, u)
 **************************************************************************/

void Solver::rhs(real t, real *udata, real *dudata)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: Solver::res(%e)", t);
#endif

  real tstart = MPI_Wtime();
  
  // Load state from udata
  load_vars(udata);
  
  // Call RHS function
  (*func)(t);
  
  // Save derivatives to dudata
  save_derivs(dudata);
  
  rhs_wtime += MPI_Wtime() - tstart;
  rhs_ncalls++;

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void Solver::pre(real t, real gamma, real delta, real *udata, real *rvec, real *zvec)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Running preconditioner: Solver::pre(%e)", t);
#endif

  real tstart = MPI_Wtime();

  int N = NV_LOCLENGTH_P(uvec);
  
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
  
  (*prefunc)(t, gamma, delta);

  // Save the solution from vars
  save_vars(zvec);

  pre_Wtime += MPI_Wtime() - tstart;
  pre_ncalls++;

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * Jacobian-vector multiplication function
 **************************************************************************/

void Solver::jac(real t, real *ydata, real *vdata, real *Jvdata)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Running Jacobian: Solver::jac(%e)", t);
#endif
  
  if(jacfunc == NULL)
    bout_error("ERROR: No jacobian function supplied!\n");
  
  // Load state from ydate
  load_vars(ydata);
  
  // Load vector to be multiplied into F_vars
  load_derivs(vdata);
  
  // Call function
  (*jacfunc)(t);

  // Save Jv from vars
  save_vars(Jvdata);

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * PRIVATE FUNCTIONS
 **************************************************************************/

/// Perform an operation at a given (jx,jy) location, moving data between BOUT++ and CVODE
void Solver::loop_vars_op(int jx, int jy, real *udata, int &p, SOLVER_VAR_OP op)
{
  real **d2d, ***d3d;
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
    
    for (jz=0; jz < ncz; jz++) {
      
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
    
    for (jz=0; jz < ncz; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	d3d = f3d[i].F_var->getData(); // Get pointer to data
	d3d[jx][jy][jz] = udata[p];
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
    
    for (jz=0; jz < ncz; jz++) {
      
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
    
    for (jz=0; jz < ncz; jz++) {
      
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
void Solver::loop_vars(real *udata, SOLVER_VAR_OP op)
{
  int jx, jy;
  int p = 0; // Counter for location in udata array

  // Inner X boundary
  if(IDATA_DEST == -1) {
    for(jx=0;jx<MXG;jx++)
      for(jy=0;jy<MYSUB;jy++)
	loop_vars_op(jx, jy+MYG, udata, p, op);
  }

  for (jx=MXG; jx < MXSUB+MXG; jx++) {
    
    // Lower Y boundary region
    
    if( ((DDATA_INDEST == -1) && (jx < DDATA_XSPLIT)) ||
	((DDATA_OUTDEST == -1) && (jx >= DDATA_XSPLIT)) ) {
      for(jy=0;jy<MYG;jy++)
	loop_vars_op(jx, jy, udata, p, op);
    }
    
    for (jy=0; jy < MYSUB; jy++) {
      // Bulk of points
      loop_vars_op(jx, jy+MYG, udata, p, op);
    }

    // Upper Y boundary condition

    if( ((UDATA_INDEST == -1) && (jx < UDATA_XSPLIT)) ||
	((UDATA_OUTDEST == -1) && (jx >= UDATA_XSPLIT)) ) {
      for(jy=0;jy<MYG;jy++)
	loop_vars_op(jx, MYSUB+MYG+jy, udata, p, op);
    }
  }

  // Outer X boundary
  if(ODATA_DEST == -1) {
    for(jx=0;jx<MXG;jx++)
      for(jy=0;jy<MYSUB;jy++)
	loop_vars_op(MXG+MXSUB+jx, jy+MYG, udata, p, op);
  }
}

void Solver::load_vars(real *udata)
{
  unsigned int i;
  
  // Make sure data is allocated
  for(i=0;i<f2d.size();i++)
    f2d[i].var->Allocate();
  for(i=0;i<f3d.size();i++) {
    f3d[i].var->Allocate();
    f3d[i].var->setLocation(f3d[i].location);
  }

  loop_vars(udata, LOAD_VARS);

  // Mark each vector as either co- or contra-variant

  for(i=0;i<v2d.size();i++)
    v2d[i].var->covariant = v2d[i].covariant;
  for(i=0;i<v3d.size();i++)
    v3d[i].var->covariant = v3d[i].covariant;
}

void Solver::load_derivs(real *udata)
{
  unsigned int i;
  
  // Make sure data is allocated
  for(i=0;i<f2d.size();i++)
    f2d[i].F_var->Allocate();
  for(i=0;i<f3d.size();i++) {
    f3d[i].F_var->Allocate();
    f3d[i].F_var->setLocation(f3d[i].location);
  }

  loop_vars(udata, LOAD_DERIVS);

  // Mark each vector as either co- or contra-variant

  for(i=0;i<v2d.size();i++)
    v2d[i].F_var->covariant = v2d[i].covariant;
  for(i=0;i<v3d.size();i++)
    v3d[i].F_var->covariant = v3d[i].covariant;
}

// This function only called during initialisation
int Solver::save_vars(real *udata)
{
  unsigned int i;

  for(i=0;i<f2d.size();i++)
    if(f2d[i].var->getData() == (real**) NULL)
      return(1);

  for(i=0;i<f3d.size();i++)
    if(f3d[i].var->getData() == (real***) NULL)
      return(1);
  
  // Make sure vectors in correct basis
  for(i=0;i<v2d.size();i++) {
    if(v2d[i].covariant) {
      v2d[i].var->to_covariant();
    }else
      v2d[i].var->to_contravariant();
  }
  for(i=0;i<v3d.size();i++) {
    if(v3d[i].covariant) {
      v3d[i].var->to_covariant();
    }else
      v3d[i].var->to_contravariant();
  }

  loop_vars(udata, SAVE_VARS);

  return(0);
}

void Solver::save_derivs(real *dudata)
{
  unsigned int i;

  // Make sure vectors in correct basis
  for(i=0;i<v2d.size();i++) {
    if(v2d[i].covariant) {
      v2d[i].F_var->to_covariant();
    }else
      v2d[i].F_var->to_contravariant();
  }
  for(i=0;i<v3d.size();i++) {
    if(v3d[i].covariant) {
      v3d[i].F_var->to_covariant();
    }else
      v3d[i].F_var->to_contravariant();
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
 * CVODE RHS functions
 **************************************************************************/

static int cvode_rhs(real t, 
		     N_Vector u, N_Vector du, 
		     void *user_data)
{
  real *udata = NV_DATA_P(u);
  real *dudata = NV_DATA_P(du);
  
  Solver *s = (Solver*) user_data;

  // Calculate residuals
  s->rhs(t, udata, dudata);

  return 0;
}

/// RHS function for BBD preconditioner
static int cvode_bbd_rhs(int Nlocal, real t, 
			 N_Vector u, N_Vector du, 
			 void *user_data)
{
  return cvode_rhs(t, u, du, user_data);
}

/// Preconditioner function
static int cvode_pre(real t, N_Vector yy, N_Vector yp,
		     N_Vector rvec, N_Vector zvec,
		     real gamma, real delta, int lr,
		     void *user_data, N_Vector tmp)
{
  real *udata = NV_DATA_P(yy);
  real *rdata = NV_DATA_P(rvec);
  real *zdata = NV_DATA_P(zvec);
  
  Solver *s = (Solver*) user_data;

  // Calculate residuals
  s->pre(t, gamma, delta, udata, rdata, zdata);

  return 0;
}

/// Jacobian-vector multiplication function
static int cvode_jac(N_Vector v, N_Vector Jv,
		     realtype t, N_Vector y, N_Vector fy,
		     void *user_data, N_Vector tmp)
{
  real *ydata = NV_DATA_P(y);   ///< System state
  real *vdata = NV_DATA_P(v);   ///< Input vector
  real *Jvdata = NV_DATA_P(Jv);  ///< Jacobian*vector output
  
  Solver *s = (Solver*) user_data;
  
  s->jac(t, ydata, vdata, Jvdata);
  
  return 0;
}
