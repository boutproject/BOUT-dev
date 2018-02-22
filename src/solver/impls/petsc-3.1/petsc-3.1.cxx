/**************************************************************************
 * Interface to PETSc 3.1 solver
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

#ifdef BOUT_HAS_PETSC_3_1

#include "petsc-3.1.hxx"

#include <globals.hxx>
#include <boutcomm.hxx>

#include <stdlib.h>

#include <interpolation.hxx> // Cell interpolation
#include <msg_stack.hxx>

#include <petsc.h>

#include <output.hxx>

EXTERN PetscErrorCode solver_f(TS ts, BoutReal t, Vec globalin, Vec globalout, void *f_data);

PetscSolver::PetscSolver() {
  has_constraints = false; // No constraints
  this->J = 0;
  this->matfdcoloring = 0;
}

PetscSolver::~PetscSolver() {
  if(initialised) {
    // Free CVODE memory

    VecDestroy(u);
    if (J){MatDestroy(J);}
    if (matfdcoloring){MatFDColoringDestroy(matfdcoloring);}
    TSDestroy(ts);

    initialised = false;
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int PetscSolver::init(int NOUT, BoutReal TIMESTEP) {
  TRACE("Initialising PETSc-3.1 solver");

  PetscErrorCode  ierr;
  int             neq;
  int             mudq, mldq, mukeep, mlkeep;
  bool            use_precon;
  int             precon_dimens;
  BoutReal        precon_tol;
  MPI_Comm        comm = PETSC_COMM_WORLD;

  // Save NOUT and TIMESTEP for use later
  nout = NOUT;
  tstep = TIMESTEP;

  /// Call the generic initialisation first
  Solver::init(NOUT, TIMESTEP);

  output.write("Initialising PETSc-3.1 solver\n");

  PetscInt n2d = n2Dvars();       // Number of 2D variables
  PetscInt n3d = n3Dvars();       // Number of 3D variables
  PetscInt local_N = getLocalN(); // Number of evolving variables on this processor

  /********** Get total problem size **********/
  if(MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    output_error.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }

  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3d, n2d, neq, local_N);

  ierr = VecCreate(BoutComm::get(), &u);CHKERRQ(ierr);
  ierr = VecSetSizes(u, local_N, PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);

  ////////// SAVE INITIAL STATE TO PETSc VECTOR ///////////
  // Set pointer to data array in vector u.
  BoutReal *udata;
  
  ierr = VecGetArray(u,&udata);CHKERRQ(ierr);
  save_vars(udata);
  
  ierr = VecRestoreArray(u,&udata);CHKERRQ(ierr);

  PetscTruth      J_load;
  MatStructure    J_structure; 
  PetscMPIInt     rank;
  char            load_file[PETSC_MAX_PATH_LEN];  /* jacobian input file name */
  PetscTruth      J_write=PETSC_FALSE,J_slowfd=PETSC_FALSE;
  ISColoring      iscoloring;
  
  // Create timestepper 
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = TSCreate(BoutComm::get(),&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  ierr = TSSetApplicationContext(ts, this);CHKERRQ(ierr);

  // Set user provided RHSFunction
  ierr = TSSetRHSFunction(ts,solver_f,this);CHKERRQ(ierr);
  
  //Sets the general-purpose update function called at the beginning of every time step. 
  //This function can change the time step.
  // TSSetPreStep(ts, PreStep);

  //Sets the general-purpose update function called after every time step -- Copy variables?
  // TSSetPostStep(ts, PostStep);

  ///////////// GET OPTIONS /////////////
  int MXSUB = mesh->xend - mesh->xstart + 1;

  Options *options = Options::getRoot();
  options = options->getSection("solver");
  OPTION(options, mudq, n3d*(MXSUB+2));
  OPTION(options, mldq, n3d*(MXSUB+2));
  OPTION(options, mukeep, 0);
  OPTION(options, mlkeep, 0);
  OPTION(options, use_precon, false);
  OPTION(options, precon_dimens, 50);
  OPTION(options, precon_tol, 1.0e-4);
  
  // Set Sundials tolerances
  BoutReal abstol, reltol;
  options->get("ATOL", abstol, 1.0e-12);
  options->get("RTOL", reltol, 1.0e-5);
  ierr = TSSundialsSetTolerance(ts, abstol, reltol);CHKERRQ(ierr);

  // Select Sundials Adams-Moulton or BDF method
  bool adams_moulton;
  OPTION(options, adams_moulton, false);
  if (adams_moulton) {
    output.write("\tUsing Adams-Moulton implicit multistep method\n");
    ierr = TSSundialsSetType(ts, SUNDIALS_ADAMS);CHKERRQ(ierr);
  } else {
    output.write("\tUsing BDF method\n");
    ierr = TSSundialsSetType(ts, SUNDIALS_BDF);CHKERRQ(ierr);
  }
  
  // Initial time and timestep. By default just use TIMESTEP
  BoutReal initial_tstep;
  OPTION(options, initial_tstep, TIMESTEP);
  ierr = TSSetInitialTimeStep(ts,simtime,initial_tstep);CHKERRQ(ierr);

  // Maximum number of steps
  int mxstep;
  OPTION(options, mxstep, 500); // Number of steps between outputs
  mxstep *= NOUT; // Total number of steps
  PetscReal tfinal = NOUT*TIMESTEP; // Final output time'=
  output.write("\tSet mxstep %d, tfinal %g, simtime %g\n",mxstep,tfinal,simtime);
  ierr = TSSetDuration(ts,mxstep,tfinal);CHKERRQ(ierr);

  // Set the current solution
  ierr = TSSetSolution(ts,u);CHKERRQ(ierr);    
  
  // Create RHSJacobian J
  SNES            snes;
  KSP             ksp;
  PC              pc;
  const PCType    pctype;
  PetscTruth      pcnone,sundialstype;

  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  ierr = TSSetUp(ts);CHKERRQ(ierr); // enable correct query of pctype below

  ierr = PetscTypeCompare((PetscObject)ts,TSSUNDIALS,&sundialstype);CHKERRQ(ierr);
  if (sundialstype){   
    ierr = TSSundialsGetPC(ts,&pc);CHKERRQ(ierr); // Sundials does not use SNES!
  } else {
    ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  }
  ierr = PCGetType(pc,&pctype);CHKERRQ(ierr);
  ierr = PetscTypeCompare((PetscObject)pc,PCNONE,&pcnone);CHKERRQ(ierr);
  output.write("\tSundialstype %d, PCNONE %d\n",sundialstype,pcnone);

  if (sundialstype && pcnone) return(0);

  //Create Jacobian matrix to be used by preconditioner
  output.write("\tSet preconditoner .... tstart %g, J localsize %d\n",simtime,local_N);
  ierr = PetscOptionsGetString(PETSC_NULL,"-J_load",load_file,PETSC_MAX_PATH_LEN-1,&J_load);CHKERRQ(ierr);
  if(J_load){
    PetscViewer     fd;   
    if (!rank){
      ierr = PetscPrintf(PETSC_COMM_SELF,"load Jmat ...\n");CHKERRQ(ierr);
    }
    ierr = PetscViewerBinaryOpen(comm,load_file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = MatLoad(fd,MATAIJ,&J);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);
    
  } else { // create Jacobian matrix by slow fd

    PetscInt MXSUB = mesh->xend - mesh->xstart + 1;
    PetscInt MYSUB = mesh->yend - mesh->ystart + 1;

    PetscInt nx = mesh->xend;//MXSUB;
    PetscInt ny = mesh->yend;//MYSUB;

    /* number of z points */
    PetscInt nz  = mesh->LocalNz;

    /* number of degrees (variables) at each grid point */
    if(n2Dvars() != 0) {
      throw BoutException("PETSc solver can't handle 2D variables yet. Sorry\n");
    }
    
    PetscInt dof = n3Dvars();

    /* Stencil width. Hardcoded to 2 until there is a public method to get mesh->MXG */
    PetscInt n = local_N; //mesh->xend*mesh->yend*nz*dof; //<- that doesn't seem to work. Why is n3Dvars()*nz?
    // PetscInt n = MXSUB*MYSUB*nz*dof;
    PetscInt sw = 2;
    PetscInt dim = 3;
    PetscInt cols = sw*2*3+1;
    PetscInt prealloc = cols*dof;
    PetscInt preallocblock = prealloc*dof;

    ierr = MatCreate(comm,&J);CHKERRQ(ierr);
    ierr = MatSetType(J, MATBAIJ);CHKERRQ(ierr);
    output << "n: " << n << "\t\t local_N: " << local_N << endl;
    ierr = MatSetSizes(J,local_N, local_N, PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);CHKERRQ(ierr);

    // Get nonzero pattern of J - color_none !!!
    ierr = MatSeqAIJSetPreallocation(J,prealloc,PETSC_NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(J,prealloc,PETSC_NULL,prealloc,PETSC_NULL);CHKERRQ(ierr);
    
    ierr = MatSeqBAIJSetPreallocation(J,dof,prealloc,PETSC_NULL);CHKERRQ(ierr);   
    ierr = MatMPIBAIJSetPreallocation(J,dof,prealloc,PETSC_NULL,prealloc,PETSC_NULL);CHKERRQ(ierr);
    ierr = MatSeqSBAIJSetPreallocation(J,dof,prealloc,PETSC_NULL);CHKERRQ(ierr);
    ierr = MatMPISBAIJSetPreallocation(J,dof,prealloc,PETSC_NULL,prealloc,PETSC_NULL);CHKERRQ(ierr);
    
    ierr = PetscOptionsHasName(PETSC_NULL,"-J_slowfd",&J_slowfd);CHKERRQ(ierr);
    if (J_slowfd){ // create Jacobian matrix by slow fd
      ierr = PetscPrintf(PETSC_COMM_SELF,"compute Jmat by slow fd...\n");CHKERRQ(ierr);
      ierr = TSDefaultComputeJacobian(ts,simtime,u,&J,&J,&J_structure,this);CHKERRQ(ierr);
    } else { // get sparse pattern of the Jacobian
      throw BoutException("Path followed in PETSc solver not yet implemented in general "
                          "-- experimental hard coded values here. Sorry\n");

      ierr = PetscPrintf(PETSC_COMM_SELF,"get sparse pattern of the Jacobian...\n");CHKERRQ(ierr);
      
      ISLocalToGlobalMapping ltog, ltogb;
      PetscInt i, j, k, d, s;
      PetscInt gi, gj;

      MatStencil stencil[cols];
      
      PetscScalar one[preallocblock];
      for (i = 0; i < preallocblock; i++) {
        one[i] = 1.0;
      }

      // Change this block for parallel. Currently this is just the identity
      // map since advect1d has no branch cuts (and we are only testing
      // single processor now)
      PetscInt ltog_array[n];
      for (i = 0; i < n; i++) {
        ltog_array[i] = i;
      }

      // Also change this for parallel. This defines the 'global stencil'
      // where starts are the starting index in each dimension and dims
      // are the size
      PetscInt starts[3], dims[3];
      starts[0] = starts[1] = starts[2] = 0;
      dims[0] = nx;
      dims[1] = ny;
      dims[2] = nz;

      // This doesn't need to be changed for parallel if ltog_array, starts, and dims are all correctly set in parallel
      ierr = ISLocalToGlobalMappingCreate(comm, n, ltog_array, &ltog);CHKERRQ(ierr);
      ierr = ISLocalToGlobalMappingBlock(ltog, dof, &ltogb);CHKERRQ(ierr);

      ierr = MatSetBlockSize(J, dof);CHKERRQ(ierr);
      ierr = MatSetLocalToGlobalMapping(J, ltog);CHKERRQ(ierr);
      ierr = MatSetLocalToGlobalMappingBlock(J, ltogb);CHKERRQ(ierr);
      ierr = MatSetStencil(J, dim, dims, starts, dof);CHKERRQ(ierr);

      bool xperiodic = false;
      // Need to figure out how to tell if y is periodic
      bool yperiodic = true;
      
      for(k=0;k<nz;k++) {
        output << "----- " << k << " -----" << endl;
        for(j=mesh->ystart; j <= mesh->yend; j++) {
          // cout << "j " << mesh->YGLOBAL(j) << ": ";
          gj = mesh->YGLOBAL(j);
          
          // X range depends on whether there are X boundaries
          int xmin = mesh->xstart;
          if(mesh->firstX())
            xmin = 0; // This processor includes a boundary region
          int xmax = mesh->xend;
          if(mesh->lastX())
            xmax = mesh->LocalNx-1;
            
          for(i=xmin; i <= xmax; i++) {
            gi = mesh->XGLOBAL(i);
            
            // Check if X and Y are periodic
            yperiodic = mesh->periodicY(i);
            xperiodic = mesh->periodicX;
            
            d = 0;
            stencil[d].k = k;
            stencil[d].j = gj;
            stencil[d].i = gi;
            stencil[d].c = dof;
            for (s = 0; s < sw; s++) {
              d++;
              stencil[d].k = k;
              stencil[d].j = gj;
              if(xperiodic) {
                stencil[d].i = (gi+s+1 >= nx) ? s-(nx-1-gi) : gi+s+1;
              } else {
                stencil[d].i = (gi+s+1 >= nx) ? -1 : gi+s+1;
              }
              stencil[d].c = dof;

              d++;
              stencil[d].k = k;
              stencil[d].j = gj;
              if(xperiodic) {
                stencil[d].i = (gi-s-1 < 0) ? nx-1-(s-gi) : gi-s-1;
              } else {
                stencil[d].i = gi-s-1;
              }
              stencil[d].c = dof;
            }
            for (s = 0; s < sw; s++) {
              d++;
              stencil[d].k = k;
              if(yperiodic) {
                stencil[d].j = (gj+s+1 >= ny) ? s-(ny-1-gj) : gj+s+1;
              } else {
                stencil[d].j = (gj+s+1 >= ny) ? -1 : gj+s+1;
              }
              stencil[d].i = gi;
              stencil[d].c = dof;

              d++;
              stencil[d].k = k;
              if(yperiodic) {
                stencil[d].j = (gj-s-1 < 0) ? ny-1-(s-gj) : gj-s-1;
              } else {
                stencil[d].j = gj-s-1;
              }
              stencil[d].i = gi;
              stencil[d].c = dof;
            }
            for (s = 0; s < sw; s++) {
              d++;
              stencil[d].k = (k+s+1 >= nz) ? s-(nz-1-k) : k+s+1;
              stencil[d].j = gj;
              stencil[d].i = gi;
              stencil[d].c = dof;

              d++;
              stencil[d].k = (k-s-1 < 0) ? nz-1-(s-k) : k-s-1;
              stencil[d].j = gj;
              stencil[d].i = gi;
              stencil[d].c = dof;
            }
            ierr = MatSetValuesBlockedStencil(J, 1, stencil, cols, stencil, one, INSERT_VALUES);CHKERRQ(ierr);

          }
        }
      }

      ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      // bout_error("stopping");
    }
  }
    
  // Write J in binary for study
  ierr = PetscOptionsHasName(PETSC_NULL,"-J_write",&J_write);CHKERRQ(ierr);
  if (J_write){ /* write J into a binary file for viewing its data structure */
    PetscViewer    viewer;
    ierr = PetscPrintf(comm,"[%d] writing J in binary to data_petsc/J.dat...\n",rank);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(comm,"data_petsc/J.dat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = MatView(J,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  }
   
  // Create coloring context of J to be used during time stepping 
  ierr = MatGetColoring(J,MATCOLORING_SL,&iscoloring);CHKERRQ(ierr); 
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

  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);   // enable PETSc runtime options

  return(0);
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

PetscErrorCode PetscSolver::run()
{
  integer steps;
  BoutReal ftime;
  
  // Set when the next call to monitor is desired
  next_time = simtime + tstep;
  outputnext = false;

  PetscFunctionBegin;
  PetscFunctionReturn(TSStep(ts,&steps,&ftime));
}

/**************************************************************************
 * RHS function
 **************************************************************************/

PetscErrorCode PetscSolver::rhs(TS ts, BoutReal t, Vec udata, Vec dudata)
{
  TRACE("Running RHS: Petsc31Solver::rhs(%e)", t);

  int flag;
  BoutReal *udata_array, *dudata_array;

  PetscFunctionBegin;

  // Load state from PETSc
  VecGetArray(udata, &udata_array);
  load_vars(udata_array);
  VecRestoreArray(udata, &udata_array);

  // Call RHS function
  flag = run_rhs(t);

  // Save derivatives to PETSc
  VecGetArray(dudata, &dudata_array);
  save_derivs(dudata_array);
  VecRestoreArray(dudata, &dudata_array);

  simtime = t; // Update the simulation time

  // Decide whether to call the monitor
  if(t >= next_time) {
  //if(outputnext) {
    // NOTE: Not using if(t >= next_time) to avoid floating-point comparisons
    
    iteration++; // Increment the 'iteration' number. keeps track of outputs
    
    // Call the monitor function
    /*
    if(call_monitors(simtime, iteration, nout)) {
      // User signalled to quit
      
      PetscFunctionReturn(1);
    }
    */
   
    // Reset iteration and wall-time count
    rhs_ncalls = 0;

    outputnext = false;
    next_time = simtime + tstep; // Set the next output time
  }

  PetscFunctionReturn(0);
}

/**************************************************************************
 * PRIVATE FUNCTIONS
 **************************************************************************/


/**************************************************************************
 * Static functions which can be used for PETSc callbacks
 **************************************************************************/
#undef __FUNCT__  
#define __FUNCT__ "Petsc31Solver::solver_f"
PetscErrorCode solver_f(TS ts, BoutReal t, Vec globalin, Vec globalout, void *f_data)
{
  PetscSolver *s;
  
  PetscFunctionBegin;
  s = (PetscSolver*) f_data;
  PetscFunctionReturn(s->rhs(ts, t, globalin, globalout));
}

#undef __FUNCT__  
#define __FUNCT__ "Petsc31Solver::PreUpdate"
PetscErrorCode PreStep(TS ts) 
{
  PetscSolver *s;
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
#define __FUNCT__ "Petsc31Solver::PostUpdate"
PetscErrorCode PostStep(TS ts) 
{
  PetscFunctionReturn(0);
}

#undef __FUNCT__  

#endif
