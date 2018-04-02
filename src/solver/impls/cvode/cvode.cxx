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

#include "unused.hxx"

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

int CvodeSolver::init(int nout, BoutReal tstep) {
  TRACE("Initialising CVODE solver");

  /// Call the generic initialisation first
  if(Solver::init(nout, tstep))
    return 1;

  // Save nout and tstep for use in run
  NOUT = nout;
  TIMESTEP = tstep;

  output_progress.write("Initialising SUNDIALS' CVODE solver\n");

  // Calculate number of variables (in generic_solver)
  int local_N = getLocalN();

  // Get total problem size
  int neq;
  {TRACE("Allreduce localN -> GlobalN");
    if (MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
      throw BoutException("ERROR: MPI_Allreduce failed!\n");
    }
  }

  output_info.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
                    n3Dvars(), n2Dvars(), neq, local_N);

  // Allocate memory
  {TRACE("Allocating memory with N_VNew_Parallel");
    if((uvec = N_VNew_Parallel(BoutComm::get(), local_N, neq)) == NULL)
      throw BoutException("ERROR: SUNDIALS memory allocation failed\n");
  }

  // Put the variables into uvec
  {TRACE("Saving variables into uvec");
    save_vars(NV_DATA_P(uvec));
  }

  /// Get options
  BoutReal abstol, reltol;
  // Initialise abstolvec to nullptr to avoid compiler maybed-uninitialised warning
  N_Vector abstolvec = nullptr;
  int maxl;
  int mudq, mldq;
  int mukeep, mlkeep;
  int max_order;
  bool use_precon, use_jacobian, use_vector_abstol, stablimdet;
  BoutReal start_timestep, max_timestep;
  bool adams_moulton, func_iter; // Time-integration method
  int MXSUB = mesh->xend - mesh->xstart + 1;
  int mxsteps; // Maximum number of steps to take between outputs
  int mxorder; // Maximum lmm order to be used by the solver
  int lmm = CV_BDF;
  int iter = CV_NEWTON;

  {TRACE("Getting options");
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
	throw BoutException("ERROR: SUNDIALS memory allocation (abstol vector) failed\n");
      vector<BoutReal> f2dtols;
      vector<BoutReal> f3dtols;
      BoutReal* abstolvec_data = NV_DATA_P(abstolvec);
      for (const auto& f : f2d) {
	abstol_options = Options::getRoot()->getSection(f.name);
	abstol_options->get("abstol", tempabstol, abstol);
	f2dtols.push_back(tempabstol);
      }
      for (const auto& f : f3d) {
	abstol_options = Options::getRoot()->getSection(f.name);
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

    options->get("mxstep", mxsteps, 500);
    options->get("mxorder", mxorder, -1);
    options->get("adams_moulton", adams_moulton, false);

    if(adams_moulton) {
      // By default use functional iteration for Adams-Moulton
      lmm = CV_ADAMS;
      output_info.write("\tUsing Adams-Moulton implicit multistep method\n");
      options->get("func_iter", func_iter, true); 
    }else {
      output_info.write("\tUsing BDF method\n");
      // Use Newton iteration for BDF
      options->get("func_iter", func_iter, false); 
    }

    if(func_iter)
      iter = CV_FUNCTIONAL;
  }//End of options TRACE

  // Call CVodeCreate
  {TRACE("Calling CVodeCreate");
    if((cvode_mem = CVodeCreate(lmm, iter)) == NULL)
      throw BoutException("CVodeCreate failed\n");
  }

  {TRACE("Calling CVodeSetUserData");
    if( CVodeSetUserData(cvode_mem, this) < 0 ) // For callbacks, need pointer to solver object
      throw BoutException("CVodeSetUserData failed\n");
  }

  {TRACE("Calling CVodeInit");
    if( CVodeInit(cvode_mem, cvode_rhs, simtime, uvec) < 0 )
      throw BoutException("CVodeInit failed\n");
  }

  
  if (max_order>0) {
    TRACE("Calling CVodeSetMaxOrder");
    if ( CVodeSetMaxOrd(cvode_mem, max_order) < 0)
      throw BoutException("CVodeSetMaxOrder failed\n");
  }
   
  if (stablimdet) {
    TRACE("Calling CVodeSetstabLimDet");
    if ( CVodeSetStabLimDet(cvode_mem, stablimdet) < 0)
      throw BoutException("CVodeSetstabLimDet failed\n");
  }
  
  if (use_vector_abstol) {
    TRACE("Calling CVodeSVtolerances");
    if( CVodeSVtolerances(cvode_mem, reltol, abstolvec) < 0 )
      throw BoutException("CVodeSStolerances failed\n");
  }
  else {
    TRACE("Calling CVodeSStolerances");
    if( CVodeSStolerances(cvode_mem, reltol, abstol) < 0 )
      throw BoutException("CVodeSStolerances failed\n");
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
  if (!func_iter) {
    output_info.write("\tUsing Newton iteration\n");
    /// Set Preconditioner
    TRACE("Setting preconditioner");
    if (use_precon) {

      int prectype = PREC_LEFT;
      bool rightprec;
      options->get("rightprec", rightprec, false);
      if (rightprec)
        prectype = PREC_RIGHT;
      
      if ( CVSpgmr(cvode_mem, prectype, maxl) != CVSPILS_SUCCESS )
        throw BoutException("ERROR: CVSpgmr failed\n");

      if (!have_user_precon()) {
        output_info.write("\tUsing BBD preconditioner\n");

        if( CVBBDPrecInit(cvode_mem, local_N, mudq, mldq, 
              mukeep, mlkeep, ZERO, cvode_bbd_rhs, NULL) )
          throw BoutException("ERROR: CVBBDPrecInit failed\n");

      } else {
        output_info.write("\tUsing user-supplied preconditioner\n");

        if( CVSpilsSetPreconditioner(cvode_mem, NULL, cvode_pre) )
          throw BoutException("ERROR: CVSpilsSetPreconditioner failed\n");
      }
    }else {
      // Not using preconditioning

      output_info.write("\tNo preconditioning\n");

      if( CVSpgmr(cvode_mem, PREC_NONE, maxl) != CVSPILS_SUCCESS )
        throw BoutException("ERROR: CVSpgmr failed\n");
    }

    /// Set Jacobian-vector multiplication function

    if((use_jacobian) && (jacfunc != NULL)) {
      output_info.write("\tUsing user-supplied Jacobian function\n");

      TRACE("Setting Jacobian-vector multiply");
      if( CVSpilsSetJacTimesVecFn(cvode_mem, cvode_jac) != CVSPILS_SUCCESS )
        throw BoutException("ERROR: CVSpilsSetJacTimesVecFn failed\n");
    }else
      output_info.write("\tUsing difference quotient approximation for Jacobian\n");
  }else {
    output_info.write("\tUsing Functional iteration\n");
  }

  return 0;
}


/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int CvodeSolver::run() {
  TRACE("CvodeSolver::run()");

  if(!initialised)
    throw BoutException("CvodeSolver not initialised\n");

  for(int i=0;i<NOUT;i++) {

    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP);
    iteration++;

    /// Check if the run succeeded
    if(simtime < 0.0) {
      // Step failed
      throw BoutException("SUNDIALS CVODE timestep failed\n");
    }
    
    if (diagnose) {
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
                   static_cast<BoutReal>(nniters) / static_cast<BoutReal>(nsteps));
      output.write("    -> Linear iterations per Newton iteration: %e\n",
                   static_cast<BoutReal>(nliters) / static_cast<BoutReal>(nniters));
      output.write("    -> Preconditioner evaluations per Newton: %e\n",
                   static_cast<BoutReal>(npevals) / static_cast<BoutReal>(nniters));

      // Last step size
      BoutReal last_step;
      CVodeGetLastStep(cvode_mem, &last_step);

      // Order used in last step
      int last_order;
      CVodeGetLastOrder(cvode_mem, &last_order);

      output.write("    -> Last step size: %e, order: %d\n", last_step, last_order);
      
      // Local error test failures
      long int num_fails;
      CVodeGetNumErrTestFails(cvode_mem, &num_fails);

      // Number of nonlinear convergence failures
      long int nonlin_fails;
      CVodeGetNumNonlinSolvConvFails(cvode_mem, &nonlin_fails);
      
      output.write("    -> Local error fails: %d, nonlinear convergence fails: %d\n", num_fails, nonlin_fails);

      // Stability limit order reductions
      long int stab_lims;
      CVodeGetNumStabLimOrderReds(cvode_mem, &stab_lims);
      
      output.write("    -> Stability limit order reductions: %d\n", stab_lims);
      
    }

    /// Call the monitor function

    if(call_monitors(simtime, i, NOUT)) {
      // User signalled to quit
      break;
    }
  }

  return 0;
}

BoutReal CvodeSolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run(%e)", tout);

  MPI_Barrier(BoutComm::get());
  
  rhs_ncalls = 0;

  pre_Wtime = 0.0;
  pre_ncalls = 0.0;

  int flag;
  if (!monitor_timestep) {
    // Run in normal mode
    flag = CVode(cvode_mem, tout, uvec, &simtime, CV_NORMAL);
  } else {
    // Run in single step mode, to call timestep monitors
    BoutReal internal_time;
    CVodeGetCurrentTime(cvode_mem, &internal_time);
    while(internal_time < tout) {
      // Run another step
      BoutReal last_time = internal_time;
      flag = CVode(cvode_mem, tout, uvec, &internal_time, CV_ONE_STEP);
      
      if (flag < 0) {
        throw BoutException("ERROR CVODE solve failed at t = %e, flag = %d\n", internal_time, flag);
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

  if (flag < 0) {
    throw BoutException("ERROR CVODE solve failed at t = %e, flag = %d\n", simtime, flag);
  }

  return simtime;
}

/**************************************************************************
 * RHS function du = F(t, u)
 **************************************************************************/

void CvodeSolver::rhs(BoutReal t, BoutReal *udata, BoutReal *dudata) {
  TRACE("Running RHS: CvodeSolver::res(%e)", t);

  // Load state from udata
  load_vars(udata);

  // Get the current timestep
  // Note: CVodeGetCurrentStep updated too late in older versions
  CVodeGetLastStep(cvode_mem, &hcur);
  
  // Call RHS function
  run_rhs(t);

  // Save derivatives to dudata
  save_derivs(dudata);
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

void CvodeSolver::pre(BoutReal t, BoutReal gamma, BoutReal delta, BoutReal *udata, BoutReal *rvec, BoutReal *zvec) {
  TRACE("Running preconditioner: CvodeSolver::pre(%e)", t);

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
}

/**************************************************************************
 * Jacobian-vector multiplication function
 **************************************************************************/

void CvodeSolver::jac(BoutReal t, BoutReal *ydata, BoutReal *vdata, BoutReal *Jvdata) {
  TRACE("Running Jacobian: CvodeSolver::jac(%e)", t);
  
  if(jacfunc == NULL)
    throw BoutException("ERROR: No jacobian function supplied!\n");
  
  // Load state from ydate
  load_vars(ydata);
  
  // Load vector to be multiplied into F_vars
  load_derivs(vdata);
  
  // Call function
  (*jacfunc)(t);

  // Save Jv from vars
  save_derivs(Jvdata);
}

/**************************************************************************
 * CVODE RHS functions
 **************************************************************************/

static int cvode_rhs(BoutReal t, 
		     N_Vector u, N_Vector du, 
		     void *user_data) {
  
  BoutReal *udata = NV_DATA_P(u);
  BoutReal *dudata = NV_DATA_P(du);

  CvodeSolver *s = static_cast<CvodeSolver *>(user_data);

  // Calculate RHS function
  try {
    s->rhs(t, udata, dudata);
  } catch (BoutRhsFail &error) {
    return 1;
  }
  return 0;
}

/// RHS function for BBD preconditioner
static int cvode_bbd_rhs(CVODEINT UNUSED(Nlocal), BoutReal t, N_Vector u, N_Vector du,
                         void *user_data) {
  return cvode_rhs(t, u, du, user_data);
}

/// Preconditioner function
static int cvode_pre(BoutReal t, N_Vector yy, N_Vector UNUSED(yp), N_Vector rvec,
                     N_Vector zvec, BoutReal gamma, BoutReal delta, int UNUSED(lr),
                     void *user_data, N_Vector UNUSED(tmp)) {
  BoutReal *udata = NV_DATA_P(yy);
  BoutReal *rdata = NV_DATA_P(rvec);
  BoutReal *zdata = NV_DATA_P(zvec);

  CvodeSolver *s = static_cast<CvodeSolver *>(user_data);

  // Calculate residuals
  s->pre(t, gamma, delta, udata, rdata, zdata);

  return 0;
}

/// Jacobian-vector multiplication function
static int cvode_jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector UNUSED(fy),
                     void *user_data, N_Vector UNUSED(tmp)) {
  BoutReal *ydata = NV_DATA_P(y);    ///< System state
  BoutReal *vdata = NV_DATA_P(v);    ///< Input vector
  BoutReal *Jvdata = NV_DATA_P(Jv);  ///< Jacobian*vector output

  CvodeSolver *s = static_cast<CvodeSolver *>(user_data);

  s->jac(t, ydata, vdata, Jvdata);
  
  return 0;
}

/**************************************************************************
 * vector abstol functions
 **************************************************************************/

void CvodeSolver::set_abstol_values(BoutReal* abstolvec_data, vector<BoutReal> &f2dtols, vector<BoutReal> &f3dtols) {
  int p = 0; // Counter for location in abstolvec_data array

  // All boundaries
  for (auto &i2d : mesh->getRegion2D("RGN_BNDRY")) {
    loop_abstol_values_op(i2d, abstolvec_data, p, f2dtols, f3dtols, true);
  }
  // Bulk of points
  for (auto &i2d : mesh->getRegion2D("RGN_NOBNDRY")) {
    loop_abstol_values_op(i2d, abstolvec_data, p, f2dtols, f3dtols, false);
  }
}

void CvodeSolver::loop_abstol_values_op(Ind2D UNUSED(i2d),
                                        BoutReal *abstolvec_data, int &p,
                                        vector<BoutReal> &f2dtols,
                                        vector<BoutReal> &f3dtols, bool bndry) {
  // Loop over 2D variables
  for(vector<BoutReal>::size_type i=0; i<f2dtols.size(); i++) {
    if(bndry && !f2d[i].evolve_bndry) {
      continue;
    }
    abstolvec_data[p] = f2dtols[i];
    p++;
  }
  
  for (int jz=0; jz < mesh->LocalNz; jz++) {
    // Loop over 3D variables
    for(vector<BoutReal>::size_type i=0; i<f3dtols.size(); i++) {
      if(bndry && !f3d[i].evolve_bndry) {
        continue;
      }
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
