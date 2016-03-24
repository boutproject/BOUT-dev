/**************************************************************************
 * Experimental interface to SUNDIALS ARKode IMEX solver
 *
 * NOTE: ARKode is still in beta testing so use with cautious optimism 
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Nick Walkden, nick.walkden@ccfe.ac.uk
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

#include "arkode.hxx"

#ifdef BOUT_HAS_ARKODE

#include <boutcomm.hxx>
#include <interpolation.hxx> // Cell interpolation
#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <arkode/arkode.h>
#include <arkode/arkode_bbdpre.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <output.hxx>

#define ZERO        RCONST(0.)
#define ONE         RCONST(1.0)

#ifndef ARKODEINT
typedef int ARKODEINT;
#endif

static int arkode_rhs_e(BoutReal t, N_Vector u, N_Vector du, void *user_data);
static int arkode_rhs_i(BoutReal t, N_Vector u, N_Vector du, void *user_data);
static int arkode_rhs(BoutReal t, N_Vector u, N_Vector du, void *user_data);

static int arkode_bbd_rhs(ARKODEINT Nlocal, BoutReal t, N_Vector u, N_Vector du, 
			 void *user_data);

static int arkode_pre(BoutReal t, N_Vector yy, N_Vector yp,
		     N_Vector rvec, N_Vector zvec,
		     BoutReal gamma, BoutReal delta, int lr,
		     void *user_data, N_Vector tmp);

static int arkode_jac(N_Vector v, N_Vector Jv,
		     realtype t, N_Vector y, N_Vector fy,
		     void *user_data, N_Vector tmp);

ArkodeSolver::ArkodeSolver(Options *opts) : Solver(opts) {
  has_constraints = false; ///< This solver doesn't have constraints
  
  jacfunc = NULL;
}

ArkodeSolver::~ArkodeSolver() {
  if(initialised) {
    N_VDestroy_Parallel(uvec);
    ARKodeFree(&arkode_mem);
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int ArkodeSolver::init(bool restarting, int nout, BoutReal tstep) {
  int msg_point = msg_stack.push("Initialising ARKODE solver");

  /// Call the generic initialisation first
  if(Solver::init(restarting, nout, tstep))
    return 1;

  // Save nout and tstep for use in run
  NOUT = nout;
  TIMESTEP = tstep;

  output.write("Initialising SUNDIALS' ARKODE solver\n");

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
  int order;
  bool use_precon, use_jacobian, use_vector_abstol,set_linear;
  BoutReal start_timestep, max_timestep, min_timestep,fixed_timestep;
  bool imex,expl,impl; // Time-integration method
  int MXSUB = mesh->xend - mesh->xstart + 1;
  BoutReal cfl_frac;
  bool fixed_step;  
  
  options->get("mudq", mudq, n3Dvars()*(MXSUB+2));
  options->get("mldq", mldq, n3Dvars()*(MXSUB+2));
  options->get("mukeep", mukeep, n3Dvars()+n2Dvars());
  options->get("mlkeep", mlkeep, n3Dvars()+n2Dvars());
  options->get("ATOL", abstol, 1.0e-12);
  options->get("RTOL", reltol, 1.0e-5);
  options->get("order", order, 4);
  options->get("cfl_frac", cfl_frac, -1.0);
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

  options->get("maxl", maxl, 0);
  OPTION(options, use_precon,   false);
  OPTION(options, use_jacobian, false);
  OPTION(options, max_timestep, -1.);
  OPTION(options, min_timestep, -1.);
  OPTION(options, start_timestep, -1);
  OPTION(options, diagnose,     false);

  int mxsteps; // Maximum number of steps to take between outputs
  options->get("mxstep", mxsteps, 500);

  int mxorder; // Maximum lmm order to be used by the solver
  options->get("mxorder", mxorder, -1);

  options->get("imex", imex, true); //Use ImEx capability
  options->get("explicit",expl,true);//Solve only explicit part 
  options->get("implicit",impl,true);//Solve only implicit part

  msg_stack.push("Calling ARKodeCreate");
  if((arkode_mem = ARKodeCreate()) == NULL)
    throw BoutException("ARKodeCreate failed\n");
  msg_stack.pop();

  msg_stack.push("Calling ARKodeSetUserData");
  if( ARKodeSetUserData(arkode_mem, this) != ARK_SUCCESS ) // For callbacks, need pointer to solver object
    throw BoutException("ARKodeSetUserData failed\n");
  msg_stack.pop();

  if(imex) {   //Use ImEx solver 
  output.write("\tUsing ARKode ImEx solver \n"); 
  msg_stack.push("Calling ARKodeInit");
  if( ARKodeInit(arkode_mem, arkode_rhs_e, arkode_rhs_i, simtime, uvec) != ARK_SUCCESS  ) //arkode_rhs_e holds the explicit part, arkode_rhs_i holds the implicit part
    throw BoutException("ARKodeInit failed\n");
  msg_stack.pop();
  
  msg_stack.push("Calling ARKodeSetImEx");
    if(expl && impl){
      if( ARKodeSetImEx(arkode_mem) != ARK_SUCCESS  )
        throw BoutException("ARKodeSetImEx failed\n");
      
  } else if(expl){
        msg_stack.push("Calling ARKodeSetExplicit");
        if( ARKodeSetExplicit(arkode_mem) != ARK_SUCCESS )
          throw BoutException("ARKodeSetExplicit failed\n");
   } else {
     msg_stack.push("Calling ARKodeSetImplicit");
        if( ARKodeSetImplicit(arkode_mem) != ARK_SUCCESS )
          throw BoutException("ARKodeSetExplicit failed\n");
        }
  } else { 

    if(expl){	//Use purely explicit solver
	output.write("\tUsing ARKode Explicit solver \n");
      msg_stack.push("Calling ARKodeInit");
      if( ARKodeInit(arkode_mem, arkode_rhs, NULL, simtime, uvec) != ARK_SUCCESS  ) //arkode_rhs_e holds the explicit part, arkode_rhs_i holds the implicit part
        throw BoutException("ARKodeInit failed\n");
      msg_stack.pop();
      msg_stack.push("Calling ARKodeSetExplicit");
      if( ARKodeSetExplicit(arkode_mem) != ARK_SUCCESS )
        throw BoutException("ARKodeSetExplicit failed\n");
    } else {	//Use purely implicit solver
	output.write("\tUsing ARKode Implicit solver \n");
      msg_stack.push("Calling ARKodeInit");
      if( ARKodeInit(arkode_mem, NULL, arkode_rhs, simtime, uvec) != ARK_SUCCESS  ) //arkode_rhs_e holds the explicit part, arkode_rhs_i holds the implicit part
        throw BoutException("ARKodeInit failed\n");
      msg_stack.pop();
      msg_stack.push("Calling ARKodeSetImplicit");
      if( ARKodeSetImplicit(arkode_mem) != ARK_SUCCESS )
        throw BoutException("ARKodeSetExplicit failed\n");
   }
 }     
  msg_stack.pop();

  OPTION(options,set_linear,false);
  if(set_linear){  //Use linear implicit solver (only evaluates jacobian inversion once
	output.write("\tSetting ARKode implicit solver to Linear\n");
	if( ARKodeSetLinear(arkode_mem,1) != ARK_SUCCESS )
  		throw BoutException("ARKodeSetLinear failed\n");
  }

  OPTION(options,fixed_step,false);	//Solve explicit portion in fixed timestep mode
					//NOTE: This is not recommended except for code comparison 
  if(fixed_step){
	options->get("timestep",fixed_timestep,0.0);	//If not given, default to adaptive timestepping
	msg_stack.push("Calling ARKodeSetFixedStep");
	if( ARKodeSetFixedStep(arkode_mem,fixed_timestep) != ARK_SUCCESS )
		throw BoutException("ARKodeSetFixedStep failed\n");
	msg_stack.pop();
  }



  msg_stack.push("Calling ARKodeSetOrder");
    if ( ARKodeSetOrder(arkode_mem, order) != ARK_SUCCESS)
      throw BoutException("ARKodeSetOrder failed\n");
  msg_stack.pop();  

  msg_stack.push("Calling ARKodeSetCFLFraction");
    if( ARKodeSetCFLFraction(arkode_mem, cfl_frac) != ARK_SUCCESS)
        throw BoutException("ARKodeSetCFLFraction failed\n");
  msg_stack.pop();  
 
  //Set timestep adaptivity function
  int adap_method;
  OPTION(options,adap_method,0);
  // 0 -> PID adaptivity (default)
  // 1 -> PI 
  // 2 -> I
  // 3 -> explicit Gustafsson
  // 4 -> implicit Gustafsson
  // 5 -> ImEx Gustafsson
 
  msg_stack.push("Calling ARKodeSetAdaptivityMethod");
    if ( ARKodeSetAdaptivityMethod(arkode_mem, adap_method,1,1,NULL) != ARK_SUCCESS)
        throw BoutException("ARKodeSetAdaptivityMethod failed\n");
  msg_stack.pop();

  if (use_vector_abstol) {
    msg_stack.push("Calling ARKodeSVtolerances");
    if( ARKodeSVtolerances(arkode_mem, reltol, abstolvec) != ARK_SUCCESS )
      throw BoutException("ARKodeSVtolerances failed\n");
  msg_stack.pop();
  }
  else {
    if( ARKodeSStolerances(arkode_mem, reltol, abstol) != ARK_SUCCESS )
      throw BoutException("ARKodeSStolerances failed\n");
    msg_stack.pop();
  }

  msg_stack.push("Calling ARKodeSetMaxNumSteps");
  if( ARKodeSetMaxNumSteps(arkode_mem, mxsteps) != ARK_SUCCESS )
     throw BoutException("ARKodeSetMaxNumSteps failed\n");
  msg_stack.pop(); 

  if(max_timestep > 0.0) {
  msg_stack.push("Calling ARKodeSetMaxStep");
    if( ARKodeSetMaxStep(arkode_mem, max_timestep) != ARK_SUCCESS )
        throw BoutException("ARKodeSetMaxStep failed\n");
  msg_stack.pop();
  }

  if(min_timestep > 0.0) {
  msg_stack.push("Calling ARKodeSetMinStep");
    if( ARKodeSetMinStep(arkode_mem, min_timestep) != ARK_SUCCESS )
      throw BoutException("ARKodeSetMinStep failed\n");
  msg_stack.pop();
  }
 
  if(start_timestep > 0.0) {
  msg_stack.push("Calling ARKodeSetInitStep"); 
    if( ARKodeSetInitStep(arkode_mem, start_timestep) != ARK_SUCCESS )
       throw BoutException("ARKodeSetInitStep failed");
  msg_stack.pop();
  }

 
  //ARKodeSetPredictorMethod(arkode_mem,4); 

  /// Newton method can include Preconditioners and Jacobian function
  bool fixed_point;
  OPTION(options,fixed_point,false);

  if(fixed_point){	//Use accellerated fixed point
  output.write("\tUsing accellerated fixed point solver\n");
    if( ARKodeSetFixedPoint(arkode_mem, 3.0) )
        throw BoutException("ARKodeSetFixedPoint failed\n");
  }else{
    output.write("\tUsing Newton iteration\n");
    if( ARKodeSetNewton(arkode_mem) )
        throw BoutException("ARKodeSetNewton failed\n");
  }

    /// Set Preconditioner
    msg_stack.push("Setting preconditioner");
    if(use_precon) {

      int prectype = PREC_LEFT;
      bool rightprec;
      options->get("rightprec", rightprec, false);
      if(rightprec)
        prectype = PREC_RIGHT;
      
      if( ARKSpgmr(arkode_mem, prectype, maxl) != ARKSPILS_SUCCESS )
        bout_error("ERROR: ARKSpgmr failed\n");

      if(!have_user_precon()) {
        output.write("\tUsing BBD preconditioner\n");

        if( ARKBBDPrecInit(arkode_mem, local_N, mudq, mldq, 
              mukeep, mlkeep, ZERO, arkode_bbd_rhs, NULL) != ARKSPILS_SUCCESS )
          bout_error("ERROR: ARKBBDPrecInit failed\n");

      } else {
        output.write("\tUsing user-supplied preconditioner\n");

        if( ARKSpilsSetPreconditioner(arkode_mem, NULL, arkode_pre) != ARKSPILS_SUCCESS )
          bout_error("ERROR: ARKSpilsSetPreconditioner failed\n");
      }
    }else {
      // Not using preconditioning

      output.write("\tNo preconditioning\n");

      if( ARKSpgmr(arkode_mem, PREC_NONE, maxl) != ARKSPILS_SUCCESS )
        bout_error("ERROR: ARKSpgmr failed\n");
    }
    msg_stack.pop();
   
 
    /// Set Jacobian-vector multiplication function

    if((use_jacobian) && (jacfunc != NULL)) {
      output.write("\tUsing user-supplied Jacobian function\n");

      msg_stack.push("Setting Jacobian-vector multiply");
      if( ARKSpilsSetJacTimesVecFn(arkode_mem, arkode_jac) != ARKSPILS_SUCCESS )
        bout_error("ERROR: ARKSpilsSetJacTimesVecFn failed\n");

      msg_stack.pop();
    }else
      output.write("\tUsing difference quotient approximation for Jacobian\n"); 

  

 
//Use ARKode optimal parameters
  bool optimize;
  OPTION(options,optimize,false);
  msg_stack.push("Calling ARKodeSetOptimialParams");
  if(optimize){
	output.write("\tUsing ARKode inbuilt optimization\n");
	if( ARKodeSetOptimalParams(arkode_mem) != ARK_SUCCESS )
		throw BoutException("ARKodeSetOptimalParams failed");
	}
  msg_stack.pop();
 
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}


/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int ArkodeSolver::run() {
#ifdef CHECK
  int msg_point = msg_stack.push("ArkodeSolver::run()");
#endif

  

  if(!initialised)
    throw BoutException("ArkodeSolver not initialised\n");

  for(int i=0;i<NOUT;i++) {

    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP);
    iteration++;

    /// Check if the run succeeded
    if(simtime < 0.0) {
      // Step failed
      output.write("Timestep failed. Aborting\n");
      
      throw BoutException("ARKode timestep failed\n");
    }
    
    if(diagnose) {
      // Print additional diagnostics
      long int nsteps, nfe_evals, nfi_evals, nniters, npevals, nliters;
      
      ARKodeGetNumSteps(arkode_mem, &nsteps);
      ARKodeGetNumRhsEvals(arkode_mem, &nfe_evals, &nfi_evals);
      ARKodeGetNumNonlinSolvIters(arkode_mem, &nniters);
      ARKSpilsGetNumPrecEvals(arkode_mem, &npevals);
      ARKSpilsGetNumLinIters(arkode_mem, &nliters);

      output.write("\nARKODE: nsteps %ld, nfe_evals %ld, nfi_evals %ld, nniters %ld, npevals %ld, nliters %ld\n", 
                   nsteps, nfe_evals, nfi_evals, nniters, npevals, nliters);
      
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

BoutReal ArkodeSolver::run(BoutReal tout) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running solver: solver::run(%e)", tout);
#endif

  MPI_Barrier(BoutComm::get());

  rhs_ncalls = 0;
  rhs_ncalls_i = 0;
  rhs_ncalls_e = 0;

  pre_Wtime = 0.0;
  pre_ncalls = 0.0;

  int flag;
  if(!monitor_timestep) {
    // Run in normal mode
    flag = ARKode(arkode_mem, tout, uvec, &simtime, ARK_NORMAL);
  }else {
    // Run in single step mode, to call timestep monitors
    BoutReal internal_time;
    ARKodeGetCurrentTime(arkode_mem, &internal_time);
    while(internal_time < tout) {
      // Run another step
      BoutReal last_time = internal_time;
      flag = ARKode(arkode_mem, tout, uvec, &internal_time, ARK_ONE_STEP);
      
      if(flag != ARK_SUCCESS) {
        output.write("ERROR ARKODE solve failed at t = %e, flag = %d\n", internal_time, flag);
        return -1.0;
      }
      
      // Call timestep monitor
      call_timestep_monitors(internal_time, internal_time - last_time);
    }
    // Get output at the desired time
    flag = ARKodeGetDky(arkode_mem, tout, 0, uvec);
    simtime = tout;
  }

  // Copy variables
  load_vars(NV_DATA_P(uvec));
  // Call rhs function to get extra variables at this time
  run_rhs(simtime);
  //run_diffusive(simtime);
  if(flag != ARK_SUCCESS) {
    output.write("ERROR ARKODE solve failed at t = %e, flag = %d\n", simtime, flag);
    return -1.0;
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return simtime;
}

/**************************************************************************
 * Explicit RHS function du = F_E(t, u)
 **************************************************************************/

void ArkodeSolver::rhs_e(BoutReal t, BoutReal *udata, BoutReal *dudata) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: ArkodeSolver::rhs_e(%e)", t);
#endif

  // Load state from udata
  load_vars(udata);

  // Get the current timestep
  // Note: ARKodeGetCurrentStep updated too late in older versions
  ARKodeGetLastStep(arkode_mem, &hcur);
  
  // Call RHS function
  run_convective(t);

  // Save derivatives to dudata
  save_derivs(dudata);

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}


/**************************************************************************
 *   Implicit RHS function du = F_I(t, u)
 **************************************************************************/

void ArkodeSolver::rhs_i(BoutReal t, BoutReal *udata, BoutReal *dudata) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: ArkodeSolver::rhs_i(%e)", t);
#endif

  load_vars(udata);
  ARKodeGetLastStep(arkode_mem, &hcur);
  // Call Implicit RHS function
  run_diffusive(t);
  save_derivs(dudata);
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
  }
  
  
/**************************************************************************
 *   Full  RHS function du = F(t, u)
 **************************************************************************/
void ArkodeSolver::rhs(BoutReal t, BoutReal *udata, BoutReal *dudata) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: ArkodeSolver::rhs(%e)", t);
#endif

  load_vars(udata);
  ARKodeGetLastStep(arkode_mem, &hcur);
  // Call Implicit RHS function
  run_rhs(t);
  save_derivs(dudata);
  
  #ifdef CHECK
     msg_stack.pop(msg_point);
  #endif
   }
  



/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void ArkodeSolver::pre(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal *udata, BoutReal *rvec, BoutReal *zvec) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running preconditioner: ArkodeSolver::pre(%e)", t);
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

  // Load vector t:q
  // o be inverted into F_vars
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

void ArkodeSolver::jac(BoutReal t, BoutReal *ydata, BoutReal *vdata, BoutReal *Jvdata) {
#ifdef CHECK
  int msg_point = msg_stack.push("Running Jacobian: ArkodeSolver::jac(%e)", t);
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
 * ARKODE explicit RHS functions
 **************************************************************************/

static int arkode_rhs_e(BoutReal t, 
		     N_Vector u, N_Vector du, 
		     void *user_data) {
  
  BoutReal *udata = NV_DATA_P(u);
  BoutReal *dudata = NV_DATA_P(du);
  
  ArkodeSolver *s = (ArkodeSolver*) user_data;
  
  // Calculate RHS function
  int rhs_e_status = 0;
  try {
    s->rhs_e(t, udata, dudata);
  }
  catch (BoutRhsFail error) {
    return 1;
  }
  return 0;
}



static int arkode_rhs_i(BoutReal t,
                     N_Vector u, N_Vector du,
                     void *user_data) {

  BoutReal *udata = NV_DATA_P(u);
  BoutReal *dudata = NV_DATA_P(du);

  ArkodeSolver *s = (ArkodeSolver*) user_data;

  //Calculate RHS function
  int rhs_e_status = 0;
  try {
     s->rhs_i(t, udata, dudata);
    }
  catch (BoutRhsFail error) {               
     return 1;
    } 
  return 0;
  }
  

static int arkode_rhs(BoutReal t,
                     N_Vector u, N_Vector du,
                     void *user_data) {

  BoutReal *udata = NV_DATA_P(u);
  BoutReal *dudata = NV_DATA_P(du);

  ArkodeSolver *s = (ArkodeSolver*) user_data;

  //Calculate RHS function
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
static int arkode_bbd_rhs(ARKODEINT Nlocal, BoutReal t, 
			 N_Vector u, N_Vector du, 
			 void *user_data)
{
  return arkode_rhs_i(t, u, du, user_data);
}

/// Preconditioner function
static int arkode_pre(BoutReal t, N_Vector yy, N_Vector yp,
		     N_Vector rvec, N_Vector zvec,
		     BoutReal gamma, BoutReal delta, int lr,
		     void *user_data, N_Vector tmp)
{
  BoutReal *udata = NV_DATA_P(yy);
  BoutReal *rdata = NV_DATA_P(rvec);
  BoutReal *zdata = NV_DATA_P(zvec);
  
  ArkodeSolver *s = (ArkodeSolver*) user_data;

  // Calculate residuals
  s->pre(t, gamma, delta, udata, rdata, zdata);

  return 0;
}

/// Jacobian-vector multiplication function
static int arkode_jac(N_Vector v, N_Vector Jv,
		     realtype t, N_Vector y, N_Vector fy,
		     void *user_data, N_Vector tmp)
{
  BoutReal *ydata = NV_DATA_P(y);   ///< System state
  BoutReal *vdata = NV_DATA_P(v);   ///< Input vector
  BoutReal *Jvdata = NV_DATA_P(Jv);  ///< Jacobian*vector output
  
  ArkodeSolver *s = (ArkodeSolver*) user_data;
  
  s->jac(t, ydata, vdata, Jvdata);
  
  return 0;
}

/**************************************************************************
 * vector abstol functions
 **************************************************************************/

void ArkodeSolver::set_abstol_values(BoutReal* abstolvec_data, vector<BoutReal> &f2dtols, vector<BoutReal> &f3dtols) {
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

void ArkodeSolver::loop_abstol_values_op(int jx, int jy, BoutReal* abstolvec_data, int &p, vector<BoutReal> &f2dtols, vector<BoutReal> &f3dtols, bool bndry) {
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

#endif
