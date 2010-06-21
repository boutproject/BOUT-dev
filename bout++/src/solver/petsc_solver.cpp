/**************************************************************************
 * Interface to PVODE solver
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

#include "petsc_solver.h"

#include "globals.h"

#include <stdlib.h>

#include "boundary.h"
#include "interpolation.h" // Cell interpolation


EXTERN PetscErrorCode solver_f(TS ts, real t, Vec globalin, Vec globalout, void *f_data);

Solver::Solver()
{
  has_constraints = false; // No constraints
}

Solver::~Solver()
{
  if(initialised) {
    // Free CVODE memory
    
    // VecDestroy(u);
    // TSDestroy(ts);

    initialised = false;
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int Solver::init(rhsfunc f, int argc, char **argv, bool restarting, int NOUT, real TIMESTEP)
{
  int neq;
  int mudq, mldq, mukeep, mlkeep;
  bool use_precon;
  int precon_dimens;
  real precon_tol;

  // Save NOUT and TIMESTEP for use later
  nout = NOUT;
  tstep = TIMESTEP;

#ifdef CHECK
  int msg_point = msg_stack.push("Initialising PETSc solver");
#endif

  /// Call the generic initialisation first
  GenericSolver::init(f, argc, argv, restarting, NOUT, TIMESTEP);

  output.write("Initialising PETSc solver\n");
  
  // Set the rhs solver function
  func = f;

  int n2d = n2Dvars();       // Number of 2D variables
  int n3d = n3Dvars();       // Number of 3D variables
  int local_N = getLocalN(); // Number of evolving variables on this processor

  /********** Get total problem size **********/

  if(MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD)) {
    output.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3d, n2d, neq, local_N);

  // Allocate memory, and set problem data, initial values, tolerances

  VecCreate(MPI_COMM_WORLD, &u);
  VecSetSizes(u, local_N, neq);
  VecSetFromOptions(u);

  ////////// SAVE INITIAL STATE TO PETSc VECTOR ///////////

  // Set pointer to data array in vector u.
  real *udata;
  
  VecGetArray(u,&udata);
  if(save_vars(udata)) {
    bout_error("\tError: Initial variable value not set\n");
    return(1);
  }
  VecRestoreArray(u,&udata);

  PetscErrorCode  ierr;
  PC              pc;
  const PCType    pctype;
  PetscTruth      pcnone,J_load,flg;
  Mat             J;
  MatStructure    J_structure; 
  PetscMPIInt     rank;
  char            load_file[PETSC_MAX_PATH_LEN];  /* jacobian input file name */
  PetscInt        M1,N1;
  PetscTruth      J_write=PETSC_FALSE,J_slowfd=PETSC_FALSE;
  
  // Create timestepper 
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  TSCreate(MPI_COMM_WORLD,&ts);
  TSSetProblemType(ts,TS_NONLINEAR);
  TSSetType(ts,TSSUNDIALS);
  TSSetApplicationContext(ts, this);

  // Set user provided RHSFunction
  TSSetRHSFunction(ts,solver_f,this);
  
  //Sets the general-purpose update function called at the beginning of every time step. 
  //This function can change the time step.
  // TSSetPreStep(ts, PreStep);

  //Sets the general-purpose update function called after every time step -- Copy variables?
  // TSSetPostStep(ts, PostStep);

  ///////////// GET OPTIONS /////////////
 	int MXSUB = mesh->xend - mesh->xstart + 1;

  options.setSection("solver");
  options.get("mudq", mudq, n3d*(MXSUB+2));
  options.get("mldq", mldq, n3d*(MXSUB+2));
  options.get("mukeep", mukeep, 0);
  options.get("mlkeep", mlkeep, 0);
  options.get("use_precon", use_precon, false);
  options.get("precon_dimens", precon_dimens, 50);
  options.get("precon_tol", precon_tol, 1.0e-4);
  
  // Set tolerances
  real abstol, reltol;
  options.get("ATOL", abstol, 1.0e-12);
  options.get("RTOL", reltol, 1.0e-5);
  TSSundialsSetTolerance(ts, abstol, reltol);

  // Select Adams-Moulton or BDF method
  bool adams_moulton;
  OPTION(adams_moulton, false);
  if(adams_moulton) {
    output.write("\tUsing Adams-Moulton implicit multistep method\n");
    TSSundialsSetType(ts, SUNDIALS_ADAMS);
  }else {
    output.write("\tUsing BDF method\n");
    TSSundialsSetType(ts, SUNDIALS_BDF);
  }

  // Initial timestep. By default just use TIMESTEP
  real initial_tstep;
  OPTION(initial_tstep, TIMESTEP);
  TSSetInitialTimeStep(ts,0.0,initial_tstep);

  // Maximum number of steps
  int pvode_mxstep;
  OPTION(pvode_mxstep, 500); // Number of steps between outputs
  pvode_mxstep *= NOUT; // Total number of steps
  PetscReal tfinal = NOUT*TIMESTEP; // Final output time
  TSSetDuration(ts,pvode_mxstep,tfinal);

  /////////////////////////////////////////////////////

  TSSetTime(ts, simtime); // Set the simulation time
  TSSetSolution(ts,u);    // Set the current solution
  
  TSSetFromOptions(ts);   // Apply PETSc options
  TSSetUp(ts);

  // Create RHSJacobian J
  TSSundialsGetPC(ts,&pc);
  PCGetType(pc,&pctype);
  PetscTypeCompare((PetscObject)pc,PCNONE,&pcnone);
  if (!pcnone){ 
    output.write("\t Set preconditoner .... tstart %g, J local size %d\n",simtime,local_N);
    PetscOptionsGetString(PETSC_NULL,"-J_load",load_file,PETSC_MAX_PATH_LEN-1,&J_load);
    if(J_load){
      PetscViewer     fd;   
      if (!rank){
        PetscPrintf(PETSC_COMM_SELF,"load Jmat ...\n");
      }
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,load_file,FILE_MODE_READ,&fd);
      MatLoad(fd,MATMPIAIJ,&J);
      PetscViewerDestroy(fd);
    } else { // create Jacobian matrix by slow fd
      MatCreate(PETSC_COMM_WORLD,&J);
      MatSetSizes(J,local_N,local_N,PETSC_DECIDE,PETSC_DECIDE);
      MatSetFromOptions(J);

      // Get nonzero pattern of J - color_none !!!
      MatSeqAIJSetPreallocation(J,10,PETSC_NULL);
      MatMPIAIJSetPreallocation(J,10,PETSC_NULL,10,PETSC_NULL);

      PetscOptionsHasName(PETSC_NULL,"-J_slowfd",&J_slowfd);
      if (J_slowfd){ // create Jacobian matrix by slow fd
        PetscPrintf(PETSC_COMM_SELF,"compute Jmat by slow fd...\n");
        TSDefaultComputeJacobian(ts,simtime,u,&J,&J,&J_structure,this);
      } else { // get sparse pattern of the Jacobian
        PetscPrintf(PETSC_COMM_SELF,"get sparse pattern of the Jacobian...\n");
        int stat;

	int MXSUB = mesh->xend - mesh->xstart + 1;	
	int NVARS=n3Dvars()+n2Dvars();
        printf("NVARS=%d\n",NVARS);  
        printf("mesh->NXPE=%d\n", mesh->NXPE);
        printf("MXSUB=%d\n",MXSUB); 
/*        printf("NYPE=%d\n",NYPE);   //-number of poloidal subdomains*/
/*        printf("MYSUB=%d\n",MYSUB);  //64; //-number of poloidal polongs per subdomain*/
        printf("MZ=%d\n",mesh->ngz);     //65;    // toroidal grid size +1
        printf("MYG=%d\n",mesh->ystart);    //-poloidal guard cells
        printf("MXG=%d\n",mesh->xstart);    //-radial guard cells

/*        stat=jstruc( NVARS,  mesh->NXPE,  MXSUB,  NYPE,  MYSUB,  MZ,  MYG,  MXG);*/
      }

      PetscInt diag;
      ierr = MatMissingDiagonal(J,&flg,&diag);CHKERRQ(ierr);
      if (flg) printf("missing %d -th diagonal!\n",diag);
    }

    
    MatGetSize(J,&M1,&N1);
    if (!rank){
      PetscPrintf(PETSC_COMM_SELF,"J size: %d,%d\n",M1,N1);
    }
    PetscOptionsHasName(PETSC_NULL,"-J_write",&J_write);
    if (J_write){ /* write J into a binary file for viewing its data structure */
      PetscViewer    viewer;
      PetscPrintf(PETSC_COMM_WORLD,"[%d] writing J in binary to data_petsc/J.dat, size %d...\n",rank,M1);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"data_petsc/J.dat",FILE_MODE_WRITE,&viewer);
      MatView(J,viewer);
      PetscViewerDestroy(viewer);
    }

    ISColoring      iscoloring;
    MatFDColoring   matfdcoloring = 0;

    // Create coloring context of J to be used during time stepping 
    ierr = MatGetColoring(J,MATCOLORING_SL,&iscoloring);CHKERRQ(ierr); // error if J is type of mpiaij???
    ierr = MatFDColoringCreate(J,iscoloring,&matfdcoloring);CHKERRQ(ierr);
    ierr = ISColoringDestroy(iscoloring);CHKERRQ(ierr);
    ierr = MatFDColoringSetFunction(matfdcoloring,(PetscErrorCode (*)(void))solver_f,this);CHKERRQ(ierr);
    ierr = MatFDColoringSetFromOptions(matfdcoloring);CHKERRQ(ierr);
    ierr = TSSetRHSJacobian(ts,J,J,TSDefaultComputeJacobianColor,matfdcoloring);CHKERRQ(ierr);

#ifdef MYDEBUG     
    // test TSDefaultComputeJacobianColor()
    ierr = TSDefaultComputeJacobianColor(ts,0.0,u_petsc,&J,&J,&J_structure,matfdcoloring);CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(data.comm, "[%d] TSDefaultComputeJacobianColor is done\n",rank);
    ierr = PetscSynchronizedFlush(data.comm);CHKERRQ(ierr);
#endif
  } // if (!pcnone)

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return(0);
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

PetscErrorCode Solver::run(MonitorFunc mon)
{
  integer steps;
  real ftime;
  
  // Set when the next call to monitor is desired
  next_time = simtime + tstep;
  monitor = mon; // Store the monitor function pointer
  outputnext = false;

  PetscFunctionBegin;
  PetscFunctionReturn(TSStep(ts,&steps,&ftime));
}

/**************************************************************************
 * RHS function
 **************************************************************************/

PetscErrorCode Solver::rhs(TS ts, real t, Vec udata, Vec dudata)
{
  int flag;
  real *udata_array, *dudata_array;

  PetscFunctionBegin;
#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: Solver::rhs(%e)", t);
#endif

  real tstart = MPI_Wtime();

  // Load state from PETSc
  VecGetArray(udata, &udata_array);
  load_vars(udata_array);
  VecRestoreArray(udata, &udata_array);

  // Call RHS function
  flag = (*func)(t);

  // Save derivatives to PETSc
  VecGetArray(dudata, &dudata_array);
  save_derivs(dudata_array);
  VecRestoreArray(dudata, &dudata_array);

  // Update the wall time and number of RHS calls
  rhs_wtime += MPI_Wtime() - tstart;
  rhs_ncalls++;

  simtime = t; // Update the simulation time

  // Decide whether to call the monitor
  if(t >= next_time) {
  //if(outputnext) {
    // NOTE: Not using if(t >= next_time) to avoid floating-point comparisons
    
    iteration++; // Increment the 'iteration' number. keeps track of outputs
    
    // Call the monitor function
    if(monitor(simtime, iteration, nout)) {
      // User signalled to quit
      
      // Write restart to a different file
      char restartname[512];
      sprintf(restartname, "data/BOUT.final.%d.pdb", MYPE);
      restart.write(restartname);
      
      output.write("Monitor signalled to quit. Returning\n");
      
      PetscFunctionReturn(1);
    }
    
    // Reset iteration and wall-time count
    rhs_ncalls = 0;
    rhs_wtime = 0.0;

    outputnext = false;
    next_time = simtime + tstep; // Set the next output time
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  PetscFunctionReturn(0);
}

/**************************************************************************
 * PRIVATE FUNCTIONS
 **************************************************************************/

/// Perform an operation at a given (jx,jy) location, moving data between BOUT++ and CVODE
void Solver::loop_vars_op(int jx, int jy, real *udata, int &p, SOLVER_VAR_OP op)
{
  real **d2d, ***d3d;
  unsigned int i;
  int jz;
 
  unsigned int n2d = f2d.size();
  unsigned int n3d = f3d.size();

  switch(op) {
  case LOAD_VARS: {
    /// Load variables from PETSc into BOUT++
    
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
}

/// Loop over variables and domain. Used for all data operations for consistency
void Solver::loop_vars(real *udata, SOLVER_VAR_OP op)
{
  int jx, jy;
  int p = 0; // Counter for location in udata array

  int MYSUB = mesh->yend - mesh->ystart + 1;

  // Inner X boundary
  if(mesh->first_x()) {
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
  if(mesh->last_x()) {
    for(jx=mesh->xend+1;jx<mesh->ngx;jx++)
      for(jy=mesh->ystart;jy<=mesh->yend;jy++)
	loop_vars_op(jx, jy, udata, p, op);
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
 * Static functions which can be used for PETSc callbacks
 **************************************************************************/
#undef __FUNCT__  
#define __FUNCT__ "Solver::solver_f"
PetscErrorCode solver_f(TS ts, real t, Vec globalin, Vec globalout, void *f_data)
{
  Solver *s;
  
  PetscFunctionBegin;
  s = (Solver*) f_data;
  PetscFunctionReturn(s->rhs(ts, t, globalin, globalout));
}

#undef __FUNCT__  
#define __FUNCT__ "Solver::PreUpdate"
PetscErrorCode PreStep(TS ts) 
{
  Solver *s;
	PetscReal t, dt;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
	ierr = TSGetTime(ts, &t);CHKERRQ(ierr);
	ierr = TSGetTimeStep(ts, &dt);CHKERRQ(ierr);
  ierr = TSGetApplicationContext(ts, (void **)&s);CHKERRQ(ierr);

  output.write("Pre-update %e\n", t);

  if((t + dt) >= s->next_time) {
    // Going past next output time
    
    dt = s->next_time - t;
    s->outputnext = true;
  }else
    s->outputnext = false;

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "Solver::PostUpdate"
PetscErrorCode PostStep(TS ts) 
{
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
