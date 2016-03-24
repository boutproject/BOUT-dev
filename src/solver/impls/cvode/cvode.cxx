/**************************************************************************
 * Interface to SUNDIALS CVODE
 * 
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

#include "cvode.hxx"

#ifdef BOUT_HAS_CVODE

#include <boutcomm.hxx>
#include <interpolation.hxx> // Cell interpolation
#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <cvode/cvode.h>
#include <cvode/cvode_bbdpre.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <output.hxx>

#define ZERO        RCONST(0.)
#define ONE         RCONST(1.0)

#ifndef CVODEINT
typedef int CVODEINT;
#endif

static int cvode_rhs(BoutReal t, N_Vector u, N_Vector du, void *user_data);
static int cvode_bbd_rhs(CVODEINT Nlocal, BoutReal t, N_Vector u, N_Vector du, 
			 void *user_data);

static int cvode_pre(BoutReal t, N_Vector yy, N_Vector yp,
		     N_Vector rvec, N_Vector zvec,
		     BoutReal gamma, BoutReal delta, int lr,
		     void *user_data, N_Vector tmp);

static int cvode_jac(N_Vector v, N_Vector Jv,
		     realtype t, N_Vector y, N_Vector fy,
		     void *user_data, N_Vector tmp);

CvodeSolver::CvodeSolver(Options *opts) : Solver(opts) {
  has_constraints = false; ///< This solver doesn't have constraints
  
  jacfunc = NULL;
  
  canReset = true;
}

CvodeSolver::~CvodeSolver() {
  if(initialised) {
    N_VDestroy_Parallel(uvec);
    CVodeFree(&cvode_mem);
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int CvodeSolver::init(bool restarting, int nout, BoutReal tstep) {
  int msg_point = msg_stack.push("Initialising CVODE solver");

  /// Call the generic initialisation first
  if(Solver::init(restarting, nout, tstep))
    return 1;

  // Save nout and tstep for use in run
  NOUT = nout;
  TIMESTEP = tstep;

  output.write("Initialising SUNDIALS' CVODE solver\n");

  // Calculate number of variables (in generic_solver)
  int local_N = getLocalN();

  // Get total problem size
  msg_stack.push("Allreduce localN -> GlobalN");
  int neq;
  if(MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    output.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }
  msg_stack.pop();

  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
                n3Dvars(), n2Dvars(), neq, local_N);

  // Allocate memory

  msg_stack.push("Allocating memory with N_VNew_Parallel");
  if((uvec = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == NULL)
    throw BoutException("ERROR: SUNDIALS memory allocation failed\n");
  msg_stack.pop();

  // Put the variables into uvec
  msg_stack.push("Saving variables into uvec");
  save_vars(NV_DATA_P(uvec));
  msg_stack.pop();

  /// Get options

  msg_stack.push("Getting options");
  BoutReal abstol, reltol;
  N_Vector abstolvec;
  int maxl;
  int mudq, mldq;
  int mukeep, mlkeep;
  int max_order;
  bool use_precon, use_jacobian, use_vector_abstol, stablimdet;
  BoutReal start_timestep, max_timestep;
  bool adams_moulton, func_iter; // Time-integration method
  int MXSUB = mesh->xend - mesh->xstart + 1;

  options->get("mudq", mudq, n3Dvars()*(MXSUB+2));
  options->get("mldq", mldq, n3Dvars()*(MXSUB+2));
  options->get("mukeep", mukeep, n3Dvars()+n2Dvars());
  options->get("mlkeep", mlkeep, n3Dvars()+n2Dvars());
  options->get("ATOL", abstol, 1.0e-12);
  options->get("RTOL", reltol, 1.0e-5);
  options->get("cvode_max_order", max_order, -1);
  options->get("cvode_stability_limit_detection", stablimdet, false);
  options->get("use_vector_abstol",use_vector_abstol,false);
  if (use_vector_abstol) {
    Options *abstol_options = Options::getRoot();
    BoutReal tempabstol;
    if((abstolvec = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == NULL)
      bout_error("ERROR: SUNDIALS memory allocation (abstol vector) failed\n");
    vector<BoutReal> f2dtols;
    vector<BoutReal> f3dtols;
    BoutReal* abstolvec_data = NV_DATA_P(abstolvec);
    for (int i=0; i<f2d.size(); i++) {
      abstol_options = Options::getRoot()->getSection(f2d[i].name);
      abstol_options->get("abstol", tempabstol, abstol);
      f2dtols.push_back(tempabstol);
    }
    for (int i=0; i<f3d.size(); i++) {
      abstol_options = Options::getRoot()->getSection(f3d[i].name);
      abstol_options->get("atol", tempabstol, abstol);
      f3dtols.push_back(tempabstol);
    }
    set_abstol_values(abstolvec_data, f2dtols, f3dtols);
  }

  options->get("maxl", maxl, 5);
  OPTION(options, use_precon,   false);
  OPTION(options, use_jacobian, false);
  OPTION(options, max_timestep, -1.);
  OPTION(options, start_timestep, -1);
  OPTION(options, diagnose,     false);

  int mxsteps; // Maximum number of steps to take between outputs
  options->get("mxstep", mxsteps, 500);

  int mxorder; // Maximum lmm order to be used by the solver
  options->get("mxorder", mxorder, -1);

  options->get("adams_moulton", adams_moulton, false);

  int lmm = CV_BDF;
  if(adams_moulton) {
    // By default use functional iteration for Adams-Moulton
    lmm = CV_ADAMS;
    output.write("\tUsing Adams-Moulton implicit multistep method\n");
    options->get("func_iter", func_iter, true); 
  }else {
    output.write("\tUsing BDF method\n");
    // Use Newton iteration for BDF
    options->get("func_iter", func_iter, false); 
  }

  int iter = CV_NEWTON;
  if(func_iter)
    iter = CV_FUNCTIONAL;
  msg_stack.pop();

  // Call CVodeCreate
  msg_stack.push("Calling CVodeCreate");
  if((cvode_mem = CVodeCreate(lmm, iter)) == NULL)
    throw BoutException("CVodeCreate failed\n");
  msg_stack.pop();

  msg_stack.push("Calling CVodeSetUserData");
  if( CVodeSetUserData(cvode_mem, this) < 0 ) // For callbacks, need pointer to solver object
    throw BoutException("CVodeSetUserData failed\n");
  msg_stack.pop();

  msg_stack.push("Calling CVodeInit");
  if( CVodeInit(cvode_mem, cvode_rhs, simtime, uvec) < 0 )
    throw BoutException("CVodeInit failed\n");
  msg_stack.pop();

  msg_stack.push("Calling CVodeSetMaxOrder");
  if (max_order>0) {
    if ( CVodeSetMaxOrd(cvode_mem, max_order) < 0)
      throw BoutException("CVodeSetMaxOrder failed\n");
  }
  
  msg_stack.push("Calling CVodeSetstabLimDet");
  if (stablimdet) {
    if ( CVodeSetStabLimDet(cvode_mem, stablimdet) < 0)
      throw BoutException("CVodeSetstabLimDet failed\n");
  }
  
  if (use_vector_abstol) {
    msg_stack.push("Calling CVodeSStolerances");
    if( CVodeSVtolerances(cvode_mem, reltol, abstolvec) < 0 )
      throw BoutException("CVodeSStolerances failed\n");
  msg_stack.pop();
  }
  else {
    if( CVodeSStolerances(cvode_mem, reltol, abstol) < 0 )
      throw BoutException("CVodeSStolerances failed\n");
    msg_stack.pop();
  }

  CVodeSetMaxNumSteps(cvode_mem, mxsteps);

  if(max_timestep > 0.0) {
    // Setting a maximum timestep
    CVodeSetMaxStep(cvode_mem, max_timestep);
  }

  if(start_timestep > 0.0) {
    // Setting a user-supplied initial guess for the appropriate timestep
    CVodeSetInitStep(cvode_mem, start_timestep);
  }
  
  if(start_timestep > 0.0) {
    CVodeSetInitStep(cvode_mem, start_timestep);
  }

  if(mxorder > 0) {
    // Setting the maximum solver order
    CVodeSetMaxOrd(cvode_mem, mxorder);
  }

  /// Newton method can include Preconditioners and Jacobian function
  if(!func_iter) {
    output.write("\tUsing Newton iteration\n");
    /// Set Preconditioner
    msg_stack.push("Setting preconditioner");
    if(use_precon) {

      int prectype = PREC_LEFT;
      bool rightprec;
      options->get("rightprec", rightprec, false);
      if(rightprec)
        prectype = PREC_RIGHT;
      
      if( CVSpgmr(cvode_mem, prectype, maxl) != CVSPILS_SUCCESS )
        bout_error("ERROR: CVSpgmr failed\n");

      if(!have_user_precon()) {
        output.write("\tUsing BBD preconditioner\n");

        if( CVBBDPrecInit(cvode_mem, local_N, mudq, mldq, 
              mukeep, mlkeep, ZERO, cvode_bbd_rhs, NULL) )
          bout_error("ERROR: CVBBDPrecInit failed\n");

      } else {
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
    msg_stack.pop();

    /// Set Jacobian-vector multiplication function

    if((use_jacobian) && (jacfunc != NULL)) {
      output.write("\tUsing user-supplied Jacobian function\n");

      msg_stack.push("Setting Jacobian-vector multiply");
      if( CVSpilsSetJacTimesVecFn(cvode_mem, cvode_jac) != CVSPILS_SUCCESS )
        bout_error("ERROR: CVSpilsSetJacTimesVecFn failed\n");

      msg_stack.pop();
    }else
      output.write("\tUsing difference quotient approximation for Jacobian\n");
  }else {
    output.write("\tUsing Functional iteration\n");
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}


/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int CvodeSolver::run() {
#ifdef CHECK
  int msg_point = msg_stack.push("CvodeSolver::run()");
#endif

  if(!initialised)
    throw BoutException("CvodeSolver not initialised\n");

  for(int i=0;i<NOUT;i++) {

    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP);
    iteration++;

    /// Check if the run succeeded
    if(simtime < 0.0) {
      // Step failed
      output.write("Timestep failed. Aborting\n");
      
      throw BoutException("SUNDIALS timestep failed\n");
    }
    
    if(diagnose) {
      // Print additional diagnostics
      long int nsteps, nfevals, nniters, npevals, nliters;
      
      CVodeGetNumSteps(cvode_mem, &nsteps);
      CVodeGetNumRhsEvals(cvode_mem, &nfevals);
      CVodeGetNumNonlinSolvIters(cvode_mem, &nniters);
      CVSpilsGetNumPrecSolves(cvode_mem, &npevals);
      CVSpilsGetNumLinIters(cvode_mem, &nliters);

      output.write("\nCVODE: nsteps %ld, nfevals %ld, nniters %ld, npevals %ld, nliters %ld\n", 
                   nsteps, nfevals, nniters, npevals, nliters);
      
      output.write("    -> Newton iterations per step: %e\n", 
                   ((double) nniters) / ((double) nsteps));
      output.write("    -> Linear iterations per Newton iteration: %e\n",
                   ((double) nliters) / ((double) nniters));
      output.write("    -> Preconditioner evaluations per Newton: %e\n",
                   ((double) npevals) / ((double) nniters));
    }

    /// Call the monitor function

    if(call_monitors(simtime, i, NOUT)) {
      // User signalled to quit
      break;
    }
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}

BoutReal CvodeSolver::run(BoutReal tout) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running solver: solver::run(%e)", tout);
#endif

  MPI_Barrier(BoutComm::get());
  
  rhs_ncalls = 0;

  pre_Wtime = 0.0;
  pre_ncalls = 0.0;

  int flag;
  if(!monitor_timestep) {
    // Run in normal mode
    flag = CVode(cvode_mem, tout, uvec, &simtime, CV_NORMAL);
  }else {
    // Run in single step mode, to call timestep monitors
    BoutReal internal_time;
    CVodeGetCurrentTime(cvode_mem, &internal_time);
    while(internal_time < tout) {
      // Run another step
      BoutReal last_time = internal_time;
      flag = CVode(cvode_mem, tout, uvec, &internal_time, CV_ONE_STEP);
      
      if(flag < 0) {
        output.write("ERROR CVODE solve failed at t = %e, flag = %d\n", internal_time, flag);
        return -1.0;
      }
      
      // Call timestep monitor
      call_timestep_monitors(internal_time, internal_time - last_time);
    }
    // Get output at the desired time
    flag = CVodeGetDky(cvode_mem, tout, 0, uvec);
    simtime = tout;
  }

  // Copy variables
  load_vars(NV_DATA_P(uvec));

  // Call rhs function to get extra variables at this time
  run_rhs(simtime);

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

void CvodeSolver::rhs(BoutReal t, BoutReal *udata, BoutReal *dudata) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: CvodeSolver::res(%e)", t);
#endif

  // Load state from udata
  load_vars(udata);

  // Get the current timestep
  // Note: CVodeGetCurrentStep updated too late in older versions
  CVodeGetLastStep(cvode_mem, &hcur);
  
  // Call RHS function
  run_rhs(t);

  // Save derivatives to dudata
  save_derivs(dudata);

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void CvodeSolver::pre(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal *udata, BoutReal *rvec, BoutReal *zvec) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running preconditioner: CvodeSolver::pre(%e)", t);
#endif

  BoutReal tstart = MPI_Wtime();

  int N = NV_LOCLENGTH_P(uvec);
  
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
  
  run_precon(t, gamma, delta);

  // Save the solution from F_vars
  save_derivs(zvec);

  pre_Wtime += MPI_Wtime() - tstart;
  pre_ncalls++;

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * Jacobian-vector multiplication function
 **************************************************************************/

void CvodeSolver::jac(BoutReal t, BoutReal *ydata, BoutReal *vdata, BoutReal *Jvdata) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running Jacobian: CvodeSolver::jac(%e)", t);
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
  save_derivs(Jvdata);

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * CVODE RHS functions
 **************************************************************************/

static int cvode_rhs(BoutReal t, 
		     N_Vector u, N_Vector du, 
		     void *user_data) {
  
  BoutReal *udata = NV_DATA_P(u);
  BoutReal *dudata = NV_DATA_P(du);
  
  CvodeSolver *s = (CvodeSolver*) user_data;
  
  // Calculate RHS function
  int rhs_status = 0;
  try {
    s->rhs(t, udata, dudata);
  }
  catch (BoutRhsFail error) {
    return 1;
  }
  return 0;
}

/// RHS function for BBD preconditioner
static int cvode_bbd_rhs(CVODEINT Nlocal, BoutReal t, 
			 N_Vector u, N_Vector du, 
			 void *user_data)
{
  return cvode_rhs(t, u, du, user_data);
}

/// Preconditioner function
static int cvode_pre(BoutReal t, N_Vector yy, N_Vector yp,
		     N_Vector rvec, N_Vector zvec,
		     BoutReal gamma, BoutReal delta, int lr,
		     void *user_data, N_Vector tmp)
{
  BoutReal *udata = NV_DATA_P(yy);
  BoutReal *rdata = NV_DATA_P(rvec);
  BoutReal *zdata = NV_DATA_P(zvec);
  
  CvodeSolver *s = (CvodeSolver*) user_data;

  // Calculate residuals
  s->pre(t, gamma, delta, udata, rdata, zdata);

  return 0;
}

/// Jacobian-vector multiplication function
static int cvode_jac(N_Vector v, N_Vector Jv,
		     realtype t, N_Vector y, N_Vector fy,
		     void *user_data, N_Vector tmp)
{
  BoutReal *ydata = NV_DATA_P(y);   ///< System state
  BoutReal *vdata = NV_DATA_P(v);   ///< Input vector
  BoutReal *Jvdata = NV_DATA_P(Jv);  ///< Jacobian*vector output
  
  CvodeSolver *s = (CvodeSolver*) user_data;
  
  s->jac(t, ydata, vdata, Jvdata);
  
  return 0;
}

/**************************************************************************
 * vector abstol functions
 **************************************************************************/

void CvodeSolver::set_abstol_values(BoutReal* abstolvec_data, vector<BoutReal> &f2dtols, vector<BoutReal> &f3dtols) {
  int jx, jy;
  int p = 0; // Counter for location in abstolvec_data array

  int MYSUB = mesh->yend - mesh->ystart + 1;

  // Inner X boundary
  if(mesh->firstX() && !mesh->periodicX) {
    for(jx=0;jx<mesh->xstart;jx++)
      for(jy=0;jy<MYSUB;jy++)
	loop_abstol_values_op(jx, jy+mesh->ystart, abstolvec_data, p, f2dtols, f3dtols, true);
  }

  // Lower Y boundary region
  for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); xi++) {
    for(jy=0;jy<mesh->ystart;jy++)
      loop_abstol_values_op(*xi, jy, abstolvec_data, p, f2dtols, f3dtols, true);
  }

  // Bulk of points
  for (jx=mesh->xstart; jx <= mesh->xend; jx++)
    for (jy=mesh->ystart; jy <= mesh->yend; jy++)
      loop_abstol_values_op(jx, jy, abstolvec_data, p, f2dtols, f3dtols, false);
  
  // Upper Y boundary condition
  for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); xi++) {
    for(jy=mesh->yend+1;jy<mesh->ngy;jy++)
      loop_abstol_values_op(*xi, jy, abstolvec_data, p, f2dtols, f3dtols, true);
  }

  // Outer X boundary
  if(mesh->lastX() && !mesh->periodicX) {
    for(jx=mesh->xend+1;jx<mesh->ngx;jx++)
      for(jy=mesh->ystart;jy<=mesh->yend;jy++)
	loop_abstol_values_op(jx, jy, abstolvec_data, p, f2dtols, f3dtols, true);
  }
}

void CvodeSolver::loop_abstol_values_op(int jx, int jy, BoutReal* abstolvec_data, int &p, vector<BoutReal> &f2dtols, vector<BoutReal> &f3dtols, bool bndry) {
  // Loop over 2D variables
  for(int i=0;i<f2dtols.size();i++) {
    if(bndry && !f2d[i].evolve_bndry)
      continue;
    abstolvec_data[p] = f2dtols[i];
    p++;
  }
  
  for (int jz=0; jz < mesh->ngz-1; jz++) {
    // Loop over 3D variables
    for(int i=0;i<f3dtols.size();i++) {
      if(bndry && !f3d[i].evolve_bndry)
	continue;
      abstolvec_data[p] = f3dtols[i];
      p++;
    }  
  }
}

void CvodeSolver::resetInternalFields() {
  TRACE("CvodeSolver::resetInternalFields");
  save_vars(NV_DATA_P(uvec));
  
  if ( CVodeReInit(cvode_mem, simtime, uvec) < 0 )
    throw BoutException("CVodeReInit failed\n");
  
}

#endif
