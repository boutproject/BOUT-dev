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

using namespace pvode;

void solver_f(integer N, BoutReal t, N_Vector u, N_Vector udot, void *f_data);
void solver_gloc(integer N, BoutReal t, BoutReal* u, BoutReal* udot, void *f_data);
void solver_cfn(integer N, BoutReal t, N_Vector u, void *f_data);

const BoutReal ZERO = 0.0;

static BoutReal abstol, reltol; // addresses passed in init must be preserved
static PVBBDData pdata;

long int iopt[OPT_SIZE];
BoutReal ropt[OPT_SIZE];

PvodeSolver::PvodeSolver() : Solver()
{
  gfunc = (rhsfunc) NULL;

  has_constraints = false; ///< This solver doesn't have constraints
}

PvodeSolver::~PvodeSolver()
{
  if(initialised) {
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

int PvodeSolver::init(rhsfunc f, int argc, char **argv, bool restarting, int nout, BoutReal tstep)
{
  int mudq, mldq, mukeep, mlkeep;
  boole optIn;
  int i;
  bool use_precon;
  int precon_dimens;
  BoutReal precon_tol;

  int n2d = n2Dvars(); // Number of 2D variables
  int n3d = n3Dvars(); // Number of 3D variables

#ifdef CHECK
  int msg_point = msg_stack.push("Initialising PVODE solver");
#endif

  /// Call the generic initialisation first
  if(Solver::init(f, argc, argv, restarting, nout, tstep))
    return 1;
  
  // Save nout and tstep for use in run
  NOUT = nout;
  TIMESTEP = tstep;

  output.write("Initialising PVODE solver\n");
  
  // Set the rhs solver function
  func = f;
  if(gfunc == (rhsfunc) NULL)
    gfunc = f; // If preconditioner function not specified, use f

  int local_N = getLocalN();
  
  // Get total problem size
  int neq;
  if(MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    output.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3d, n2d, neq, local_N);

  // Set machEnv block
  machEnv = (machEnvType) PVecInitMPI(BoutComm::get(), local_N, neq, &argc, &argv);

  if (machEnv == NULL) {
    if(MYPE == 0)
      output.write("\tError: PVecInitMPI failed\n");
    return(1);
  }

  //.allocate memory, and set problem data, initial values, tolerances

  u = N_VNew(neq, machEnv);

  ///////////// GET OPTIONS /////////////

  int pvode_mxstep;
  int MXSUB = mesh->xend - mesh->xstart + 1;
  
  Options *options = Options::getRoot();
  options = options->getSection("solver");
  options->get("mudq", mudq, n3d*(MXSUB+2));
  options->get("mldq", mldq, n3d*(MXSUB+2));
  options->get("mukeep", mukeep, 0);
  options->get("mlkeep", mlkeep, 0);
  options->get("ATOL", abstol, 1.0e-12);
  options->get("RTOL", reltol, 1.0e-5);
  options->get("use_precon", use_precon, false);
  options->get("precon_dimens", precon_dimens, 50);
  options->get("precon_tol", precon_tol, 1.0e-4);
  options->get("mxstep", pvode_mxstep, 500);

  pdata = PVBBDAlloc(local_N, mudq, mldq, mukeep, mlkeep, ZERO, 
                     solver_gloc, solver_cfn, (void*) this);
  
  if (pdata == NULL) {
    output.write("\tError: PVBBDAlloc failed.\n");
    return(1);
  }

  ////////// SAVE DATA TO CVODE ///////////

  // Set pointer to data array in vector u.
  BoutReal *udata = N_VDATA(u);
  if(save_vars(udata)) {
    bout_error("\tError: Initial variable value not set\n");
    return(1);
  }
  
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

  cvode_mem = CVodeMalloc(neq, solver_f, simtime, u, BDF, NEWTON, SS, &reltol,
                          &abstol, this, NULL, optIn, iopt, ropt, machEnv);

  if(cvode_mem == NULL) {
    output.write("\tError: CVodeMalloc failed.\n");
    return(1);
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
  
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
  
  return(0);
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

int PvodeSolver::run(MonitorFunc monitor)
{
#ifdef CHECK
  int msg_point = msg_stack.push("PvodeSolver::run()");
#endif
  
  if(!initialised)
    bout_error("PvodeSolver not initialised\n");
  
  for(int i=0;i<NOUT;i++) {
    
    /// Run the solver for one output timestep
    simtime = run(simtime + TIMESTEP, rhs_ncalls, rhs_wtime);
    iteration++;

    /// Check if the run succeeded
    if(simtime < 0.0) {
      // Step failed
      output.write("Timestep failed. Aborting\n");

      // Write restart to a different file
      restart.write("%s/BOUT.failed.%d.%s", restartdir.c_str(), MYPE, restartext.c_str());

      bout_error("PVODE timestep failed\n");
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

BoutReal PvodeSolver::run(BoutReal tout, int &ncalls, BoutReal &rhstime)
{
  BoutReal *udata;
  int flag;

#ifdef CHECK
  int msg_point = msg_stack.push("Running solver: solver::run(%e)", tout);
#endif

  rhs_wtime = 0.0;
  rhs_ncalls = 0;

  // Set pointer to data array in vector u.
  udata = N_VDATA(u);

  // Run CVODE
  flag = CVode(cvode_mem, tout, u, &simtime, NORMAL);

  ncalls = rhs_ncalls;
  rhstime = rhs_wtime;

  // Copy variables
  load_vars(udata);
  
  // Call rhs function to get extra variables at this time
  BoutReal tstart = MPI_Wtime();
  (*func)(simtime);
  rhstime += MPI_Wtime() - tstart;
  ncalls++;

  // Check return flag
  if(flag != SUCCESS) {
    output.write("ERROR CVODE step failed, flag = %d\n", flag);
    return(-1.0);
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return simtime;
}

/**************************************************************************
 * RHS function
 **************************************************************************/

void PvodeSolver::rhs(int N, BoutReal t, BoutReal *udata, BoutReal *dudata)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: PvodeSolver::rhs(%e)", t);
#endif

  // Load state from CVODE
  load_vars(udata);

  // Call function
  int flag = run_rhs(t);

  // Save derivatives to CVODE
  save_derivs(dudata);

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void PvodeSolver::gloc(int N, BoutReal t, BoutReal *udata, BoutReal *dudata)
{
  int flag;
  BoutReal tstart;

#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: PvodeSolver::gloc(%e)", t);
#endif

  tstart = MPI_Wtime();

  // Load state from CVODE
  load_vars(udata);

  // Call function
  flag = (*gfunc)(t);

  // Save derivatives to CVODE
  save_derivs(dudata);

  rhs_wtime += MPI_Wtime() - tstart;
  rhs_ncalls++;

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * PRIVATE FUNCTIONS
 **************************************************************************/

/// Perform an operation at a given (jx,jy) location, moving data between BOUT++ and CVODE
void PvodeSolver::loop_vars_op(int jx, int jy, BoutReal *udata, int &p, SOLVER_VAR_OP op)
{
  BoutReal **d2d, ***d3d;
  unsigned int i;
  int jz;

#ifdef CHECK
  //msg_stack.push("loop_vars_op(jx=%d, jy=%d, op=%d)", jx, jy, op);
  
  if((jx < 0) || (jx >= mesh->ngx) ||
     (jy < 0) || (jy >= mesh->ngy)) {
    output.write("OUT OF DOMAIN: %d, %d", jx, jy);
    bout_error("ERROR: Out of range in loop_vars_op\n");
  }
#endif

  unsigned int n2d = n2Dvars();
  unsigned int n3d = n3Dvars();
 
  switch(op) {
  case LOAD_VARS: {
    /// Load variables from CVODE into BOUT++
    
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
  case SAVE_VARS: {
    /// Save variables from BOUT++ into CVODE (only used at start of simulation)
    
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
#ifdef CHECK
  //msg_stack.pop();
#endif
}

/// Loop over variables and domain. Used for all data operations for consistency
void PvodeSolver::loop_vars(BoutReal *udata, SOLVER_VAR_OP op)
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
  RangeIter *xi = mesh->iterateBndryLowerY();
  for(xi->first(); !xi->isDone(); xi->next()) {
    for(jy=0;jy<mesh->ystart;jy++)
      loop_vars_op(xi->ind, jy, udata, p, op);
  }
  delete xi;

  // Bulk of points
  for (jx=mesh->xstart; jx <= mesh->xend; jx++)
    for (jy=mesh->ystart; jy <= mesh->yend; jy++)
      loop_vars_op(jx, jy, udata, p, op);
  
  // Upper Y boundary condition
  xi = mesh->iterateBndryUpperY();
  for(xi->first(); !xi->isDone(); xi->next()) {
    for(jy=mesh->yend+1;jy<mesh->ngy;jy++)
      loop_vars_op(xi->ind, jy, udata, p, op);
  }
  delete xi;

  // Outer X boundary
  if(mesh->lastX()) {
    for(jx=mesh->xend+1;jx<mesh->ngx;jx++)
      for(jy=mesh->ystart;jy<=mesh->yend;jy++)
	loop_vars_op(jx, jy, udata, p, op);
  }
}

void PvodeSolver::load_vars(BoutReal *udata)
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

// This function only called during initialisation
int PvodeSolver::save_vars(BoutReal *udata)
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

void PvodeSolver::save_derivs(BoutReal *dudata)
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
 * CVODE rhs function
 **************************************************************************/

void solver_f(integer N, BoutReal t, N_Vector u, N_Vector udot, void *f_data)
{
  BoutReal *udata, *dudata;
  PvodeSolver *s;

  udata = N_VDATA(u);
  dudata = N_VDATA(udot);
  
  s = (PvodeSolver*) f_data;

  s->rhs(N, t, udata, dudata);
}

// Preconditioner RHS
void solver_gloc(integer N, BoutReal t, BoutReal* u, BoutReal* udot, void *f_data)
{
  PvodeSolver *s;
  
  s = (PvodeSolver*) f_data;

  s->gloc(N, t, u, udot);
}

// Preconditioner communication function
void solver_cfn(integer N, BoutReal t, N_Vector u, void *f_data)
{
  // doesn't do anything at the moment
}

#endif 
