/**************************************************************************
 * Interface to PVODE solver
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

#include "pvode.hxx"

#ifdef BOUT_HAS_PVODE

#include <bout/mesh.hxx>
#include <boutcomm.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <bout/sys/timer.hxx>
#include <boutexception.hxx>

#include "unused.hxx"

#include <pvode/iterativ.h>  // contains the enum for types of preconditioning
#include <pvode/cvspgmr.h>   // use CVSPGMR linear solver each internal step
#include <pvode/pvbbdpre.h>  // band preconditioner function prototypes

using namespace pvode;

void solver_f(integer N, BoutReal t, N_Vector u, N_Vector udot, void *f_data);
void solver_gloc(integer N, BoutReal t, BoutReal* u, BoutReal* udot, void *f_data);
void solver_cfn(integer N, BoutReal t, N_Vector u, void *f_data);

const BoutReal ZERO = 0.0;

long int iopt[OPT_SIZE];
BoutReal ropt[OPT_SIZE];

PvodeSolver::PvodeSolver(Options *options) : Solver(options) {
  has_constraints = false; ///< This solver doesn't have constraints
}

PvodeSolver::~PvodeSolver() {
  if(pvode_initialised) {
    // Free CVODE memory
    
    N_VFree(u);
    PVBBDFree(pdata);
    CVodeFree(cvode_mem);
    PVecFreeMPI(machEnv);
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int PvodeSolver::init(int nout, BoutReal tstep) {
  TRACE("Initialising PVODE solver");

  int mudq, mldq, mukeep, mlkeep;
  boole optIn;
  int i;
  bool use_precon;
  int precon_dimens;
  BoutReal precon_tol;

  int n2d = n2Dvars(); // Number of 2D variables
  int n3d = n3Dvars(); // Number of 3D variables

  /// Call the generic initialisation first
  if(Solver::init(nout, tstep))
    return 1;
  
  // Save nout and tstep for use in run
  NOUT = nout;
  TIMESTEP = tstep;

  output.write("Initialising PVODE solver\n");

  int local_N = getLocalN();

  if(local_N == 0) {
    throw BoutException("No local evolving variables");
  }
  
  // Get total problem size
  int neq;
  if(MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("\tERROR: MPI_Allreduce failed!\n");
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3d, n2d, neq, local_N);

  // Set machEnv block
  machEnv = static_cast<machEnvType>(PVecInitMPI(BoutComm::get(), local_N, neq, pargc, pargv));

  if (machEnv == nullptr) {
    throw BoutException("\tError: PVecInitMPI failed\n");
  }

  // Allocate memory, and set problem data, initial values, tolerances

  u = N_VNew(neq, machEnv);

  ///////////// GET OPTIONS /////////////

  int pvode_mxstep;
  // Compute band_width_default from actually added fields, to allow for multiple Mesh objects
  //
  // Previous implementation was equivalent to:
  //   int MXSUB = mesh->xend - mesh->xstart + 1;
  //   int band_width_default = n3Dvars()*(MXSUB+2);
  int band_width_default = 0;
  for (const auto& fvar : f3d) {
    Mesh* localmesh = fvar.var->getMesh();
    band_width_default += localmesh->xend - localmesh->xstart + 3;
  }

  
  options->get("mudq", mudq, band_width_default);
  options->get("mldq", mldq, band_width_default);
  options->get("mukeep", mukeep, 0);
  options->get("mlkeep", mlkeep, 0);
  options->get("ATOL", abstol, 1.0e-12);
  options->get("RTOL", reltol, 1.0e-5);
  options->get("use_precon", use_precon, false);
  options->get("precon_dimens", precon_dimens, 50);
  options->get("precon_tol", precon_tol, 1.0e-4);
  options->get("mxstep", pvode_mxstep, 500);

  pdata = PVBBDAlloc(local_N, mudq, mldq, mukeep, mlkeep, ZERO, 
                     solver_gloc, solver_cfn, static_cast<void*>(this));

  if (pdata == nullptr) {
    throw BoutException("\tError: PVBBDAlloc failed.\n");
  }

  ////////// SAVE DATA TO CVODE ///////////

  // Set pointer to data array in vector u.
  BoutReal *udata = N_VDATA(u);
  save_vars(udata);
  
  /* Call CVodeMalloc to initialize CVODE: 
     
     neq     is the problem size = number of equations
     f       is the user's right hand side function in y'=f(t,y)
     T0      is the initial time
     u       is the initial dependent variable vector
     BDF     specifies the Backward Differentiation Formula
     NEWTON  specifies a Newton iteration
     SS      specifies scalar relative and absolute tolerances
     &reltol and &abstol are pointers to the scalar tolerances
     data    is the pointer to the user-defined data block
     NULL    is the pointer to the error message file
     FALSE   indicates there are no optional inputs in iopt and ropt
     iopt, ropt  communicate optional integer and BoutReal input/output

     A pointer to CVODE problem memory is returned and stored in cvode_mem.  */

  optIn = TRUE; for(i=0;i<OPT_SIZE;i++)iopt[i]=0; 
                for(i=0;i<OPT_SIZE;i++)ropt[i]=ZERO;
		iopt[MXSTEP]=pvode_mxstep;

  cvode_mem = CVodeMalloc(neq, solver_f, simtime, u, BDF, NEWTON, SS, &reltol, &abstol,
                          this, nullptr, optIn, iopt, ropt, machEnv);

  if (cvode_mem == nullptr) {
    throw BoutException("\tError: CVodeMalloc failed.\n");
  }

  /* Call CVSpgmr to specify the CVODE linear solver CVSPGMR with
     left preconditioning, modified Gram-Schmidt orthogonalization,
     default values for the maximum Krylov dimension maxl and the tolerance
     parameter delt, preconditioner setup and solve routines from the
     PVBBDPRE module, and the pointer to the preconditioner data block.    */

  if(use_precon) {
    CVSpgmr(cvode_mem, LEFT, MODIFIED_GS, precon_dimens, precon_tol, PVBBDPrecon, PVBBDPSol, pdata);
  }else {
    CVSpgmr(cvode_mem, NONE, MODIFIED_GS, 10, ZERO, PVBBDPrecon, PVBBDPSol, pdata);
  }

  /*  CVSpgmr(cvode_mem, NONE, MODIFIED_GS, 10, 0.0, PVBBDPrecon, PVBBDPSol, pdata); */

  // PvodeSolver is now initialised fully
  pvode_initialised = true;
  
  return(0);
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int PvodeSolver::run() {
  TRACE("PvodeSolver::run()");
  
  if(!pvode_initialised)
    throw BoutException("PvodeSolver not initialised\n");
  
  for(int i=0;i<NOUT;i++) {
    
    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP);
    iteration++;

    /// Check if the run succeeded
    if(simtime < 0.0) {
      // Step failed
      output.write("Timestep failed. Aborting\n");
      
      throw BoutException("PVODE timestep failed\n");
    }
    
    /// Call the monitor function
    
    if(call_monitors(simtime, i, NOUT)) {
      // User signalled to quit
      break;
    }
  }
  
  return 0;
}

BoutReal PvodeSolver::run(BoutReal tout) {
  TRACE("Running solver: solver::run(%e)", tout);

  BoutReal *udata;
  
  // Set pointer to data array in vector u.
  udata = N_VDATA(u);

  // Run CVODE
  int flag;
  if(!monitor_timestep) {
    // Run in normal mode
    flag = CVode(cvode_mem, tout, u, &simtime, NORMAL);
  }else {
    // Run in single step mode, to call timestep monitors
    BoutReal internal_time = static_cast<CVodeMem>(cvode_mem)->cv_tn;
    //CvodeGetCurrentTime(cvode_mem, &internal_time);
    
    while(internal_time < tout) {
      // Run another step
      BoutReal last_time = internal_time;
      flag = CVode(cvode_mem, tout, u, &internal_time, ONE_STEP);
      if(flag < 0) {
        output_error.write("ERROR CVODE solve failed at t = %e, flag = %d\n", internal_time, flag);
        return -1.0;
      }
      
      // Call timestep monitor
      call_timestep_monitors(internal_time, internal_time - last_time);
    }
    // Get output at the desired time
    flag = CVodeDky(cvode_mem, tout, 0, u);
    simtime = tout;
  }

  // Copy variables
  load_vars(udata);
  
  // Call rhs function to get extra variables at this time
  run_rhs(simtime);

  // Check return flag
  if(flag != SUCCESS) {
    output_error.write("ERROR CVODE step failed, flag = %d\n", flag);
    return(-1.0);
  }

  return simtime;
}

/**************************************************************************
 * RHS function
 **************************************************************************/

void PvodeSolver::rhs(int UNUSED(N), BoutReal t, BoutReal *udata, BoutReal *dudata) {
  TRACE("Running RHS: PvodeSolver::rhs(%e)", t);

  // Get current timestep
  hcur = 0.0; //((CVodeMemRec*) cvode_mem)->cv_h;

  // Load state from CVODE
  load_vars(udata);

  // Call function
  run_rhs(t);

  // Save derivatives to CVODE
  save_derivs(dudata);
}

void PvodeSolver::gloc(int UNUSED(N), BoutReal t, BoutReal *udata, BoutReal *dudata) {
  TRACE("Running RHS: PvodeSolver::gloc(%e)", t);

  Timer timer("rhs");

  // Load state from CVODE
  load_vars(udata);

  // Call function
  run_rhs(t);

  // Save derivatives to CVODE
  save_derivs(dudata);
}

/**************************************************************************
 * CVODE rhs function
 **************************************************************************/

void solver_f(integer N, BoutReal t, N_Vector u, N_Vector udot, void *f_data) {
  BoutReal *udata, *dudata;
  PvodeSolver *s;

  udata = N_VDATA(u);
  dudata = N_VDATA(udot);

  s = static_cast<PvodeSolver *>(f_data);

  s->rhs(N, t, udata, dudata);
}

// Preconditioner RHS
void solver_gloc(integer N, BoutReal t, BoutReal *u, BoutReal *udot, void *f_data) {
  PvodeSolver *s;

  s = static_cast<PvodeSolver *>(f_data);

  s->gloc(N, t, u, udot);
}

// Preconditioner communication function
void solver_cfn(integer UNUSED(N), BoutReal UNUSED(t), N_Vector UNUSED(u),
                void *UNUSED(f_data)) {
  // doesn't do anything at the moment
}

#endif 
