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

#include "petsc.hxx"

#ifdef BOUT_HAS_PETSC_DEV

#include <globals.hxx>

#include <stdlib.h>

#include <interpolation.hxx> // Cell interpolation


extern PetscErrorCode solver_f(TS ts, BoutReal t, Vec globalin, Vec globalout, void *f_data);
extern PetscErrorCode solver_rhsjacobian(TS ts,BoutReal t,Vec globalin,Mat *J,Mat *Jpre,MatStructure *str,void *f_data);
extern PetscErrorCode solver_if(TS,BoutReal,Vec,Vec,Vec,void*);  
extern PetscErrorCode solver_ijacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);

PetscSolver::PetscSolver()
{
  has_constraints = false; // No constraints
  this->J = 0;
  this->matfdcoloring = 0;
}

PetscSolver::~PetscSolver()
{

  if(initialised) {
    // Free CVODE memory
    
    VecDestroy(&u);
    if (J) {MatDestroy(&J);}
    if (matfdcoloring){MatFDColoringDestroy(&matfdcoloring);}
    TSDestroy(&ts);

    initialised = false;
  }
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int PetscSolver::init(rhsfunc f, int argc, char **argv, bool restarting, int NOUT, BoutReal TIMESTEP)
{
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

#ifdef CHECK
  int msg_point = msg_stack.push("Initialising PETSc solver");
#endif

  /// Call the generic initialisation first
  Solver::init(f, argc, argv, restarting, NOUT, TIMESTEP);

  output.write("Initialising PETSc solver\n");

  PetscInt n2d = n2Dvars();       // Number of 2D variables
  PetscInt n3d = n3Dvars();       // Number of 3D variables
  PetscInt local_N = getLocalN(); // Number of evolving variables on this processor

  /********** Get total problem size **********/
  if(MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    output.write("\tERROR: MPI_Allreduce failed!\n");
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
  if(save_vars(udata)) {
    bout_error("\tError: Initial variable value not set\n");
    return(1);
  }
  ierr = VecRestoreArray(u,&udata);CHKERRQ(ierr);

  PetscBool       J_load;
  MatStructure    J_structure; 
  PetscMPIInt     rank;
  char            load_file[PETSC_MAX_PATH_LEN];  /* jacobian input file name */
  PetscBool       J_write=PETSC_FALSE,J_slowfd=PETSC_FALSE;
  ISColoring      iscoloring;
  
  // Create timestepper 
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = TSCreate(BoutComm::get(),&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSGL);CHKERRQ(ierr);
  //ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  //ierr = TSSetApplicationContext(ts, this);CHKERRQ(ierr); // do we need this?

  // Set user provided RHSFunction
  ierr = TSSetRHSFunction(ts,solver_f,this);CHKERRQ(ierr); // needed for ts_type=sundials
  ierr = TSSetIFunction(ts,solver_if,this);CHKERRQ(ierr);
  
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
  PetscBool       pcnone=PETSC_FALSE,sundialstype;

  // set default ksp and pc 
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr); 
  //pcnone = PETSC_TRUE;

#if defined(USE_SUNDIALS)
  // this block causes '-ts_type gl -snes_mf_operator' crash - further work on this! 
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  ierr = TSSetUp(ts);CHKERRQ(ierr); // for getting pctype below
  ierr = PetscTypeCompare((PetscObject)pc,PCNONE,&pcnone);CHKERRQ(ierr);

  ierr = PetscTypeCompare((PetscObject)ts,TSSUNDIALS,&sundialstype);CHKERRQ(ierr);
  if (sundialstype){ // mv sundials procedural calls here

    output.write("\tSundialstype %d, PCNONE %d\n",sundialstype,pcnone);
    if (pcnone) return(0); 
  }
#endif


  //Create Jacobian matrix to be used by preconditioner
  output.write("\tGet Jacobian matrix .... tstart %g, J localsize %d\n",simtime,local_N);
  ierr = PetscOptionsGetString(PETSC_NULL,"-J_load",load_file,PETSC_MAX_PATH_LEN-1,&J_load);CHKERRQ(ierr);
  if(J_load){
    PetscViewer     fd;   
    if (!rank){
      ierr = PetscPrintf(PETSC_COMM_SELF,"load Jmat ...\n");CHKERRQ(ierr);
    }
    ierr = PetscViewerBinaryOpen(comm,load_file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = MatLoad(J, fd);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    
  } else { // create Jacobian matrix by slow fd

    PetscInt MXSUB = mesh->xend - mesh->xstart + 1;
    PetscInt MYSUB = mesh->yend - mesh->ystart + 1;

    PetscInt nx = mesh->xend;//MXSUB;
    PetscInt ny = mesh->yend;//MYSUB;

    /* number of z points (need to subtract one because of historical reasons that MZ has an extra point) */
    PetscInt nz  = mesh->ngz - 1;

    /* number of degrees (variables) at each grid point */
    if(n2Dvars() != 0) {
      bout_error("PETSc solver can't handle 2D variables yet. Sorry\n");
    }
    
    PetscInt dof = n3Dvars();

    /* Stencil width. Hardcoded to 2 until there is a public method to get mesh->MXG */
    PetscInt n = local_N; //mesh->xend*mesh->yend*nz*dof; //<- that doesn't seem to work. Why is n3Dvars()*nz?
    // PetscInt n = MXSUB*MYSUB*nz*dof;
    PetscInt sw = 2;
    PetscInt dim = 3;
    PetscInt cols = sw*2*3+1;
    PetscInt prealloc; // = cols*dof;
    PetscInt preallocblock = cols*dof*dof; //prealloc*dof;

    ierr = MatCreate(comm,&J);CHKERRQ(ierr);
    ierr = MatSetType(J, MATBAIJ);CHKERRQ(ierr);
    cout << "n: " << n << "\t\t local_N: " << local_N << endl;
    ierr = MatSetSizes(J,local_N, local_N, PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);CHKERRQ(ierr);

    // Get nonzero pattern of J - color_none !!!
    prealloc = cols*dof*dof;
    ierr = MatSeqAIJSetPreallocation(J,prealloc,PETSC_NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(J,prealloc,PETSC_NULL,prealloc,PETSC_NULL);CHKERRQ(ierr);
    
    prealloc = cols;
    ierr = MatSeqBAIJSetPreallocation(J,dof,prealloc,PETSC_NULL);CHKERRQ(ierr);   
    ierr = MatMPIBAIJSetPreallocation(J,dof,prealloc,PETSC_NULL,prealloc,PETSC_NULL);CHKERRQ(ierr);

    ierr = TSSetIJacobian(ts,J,J,solver_ijacobian,this);CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts);CHKERRQ(ierr);   // enable PETSc runtime options
    ierr = TSSetUp(ts);CHKERRQ(ierr);            // for getting pctype below
    ierr = PetscTypeCompare((PetscObject)pc,PCNONE,&pcnone);CHKERRQ(ierr);
    if (pcnone){
      printf("pcnone %d, use dummy Jacobian... \n",pcnone);
      return(0);
    }
   
    ierr = PetscOptionsHasName(PETSC_NULL,"-J_slowfd",&J_slowfd);CHKERRQ(ierr);
    if (J_slowfd){ // create Jacobian matrix by slow fd
      ierr = TSSetRHSJacobian(ts,J,J,TSDefaultComputeJacobian,this);CHKERRQ(ierr);

      // Create coloring context of J to be used during time stepping - solver_f seems not linked to ts!!!
      if (!rank){ierr = PetscPrintf(PETSC_COMM_SELF,"compute Jmat by slow fd...\n");CHKERRQ(ierr);}
      ierr = TSComputeRHSJacobian(ts,simtime,u,&J,&J,&J_structure);CHKERRQ(ierr); // this J is correct!
#if defined(COLOR)          
      ierr = MatGetColoring(J,MATCOLORINGSL,&iscoloring);CHKERRQ(ierr); 
      ierr = MatFDColoringCreate(J,iscoloring,&matfdcoloring);CHKERRQ(ierr);
      ierr = MatFDColoringSetFromOptions(matfdcoloring);CHKERRQ(ierr);
      ierr = ISColoringDestroy(&iscoloring);CHKERRQ(ierr);

      ierr = MatFDColoringSetFunction(matfdcoloring,(PetscErrorCode (*)(void))solver_f,this);CHKERRQ(ierr);
      ierr = TSSetRHSJacobian(ts,J,J,TSDefaultComputeJacobianColor,matfdcoloring);CHKERRQ(ierr);  
  
      // Test TSDefaultComputeJacobianColor()
      ierr = PetscPrintf(comm,"\n[%d] Test TSDefaultComputeJacobianColor() ...\n",rank);
      ierr = TSComputeRHSJacobian(ts,simtime,u,&J,&J,&J_structure);CHKERRQ(ierr); // this J=0 !!!
      //ierr = TSDefaultComputeJacobianColor(ts,simtime,u,&J,&J,&J_structure, (MatFDColoring)matfdcoloring);CHKERRQ(ierr); // this J=0 !!!
      
      ierr = MatView(J,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
    } else { // get sparse pattern of the Jacobian
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
      // CHECK PETSC_COPY_VALUES
      ierr = ISLocalToGlobalMappingCreate(comm, n, ltog_array, PETSC_COPY_VALUES, &ltog);CHKERRQ(ierr);
      ierr = ISLocalToGlobalMappingBlock(ltog, dof, &ltogb);CHKERRQ(ierr);

      ierr = MatSetBlockSize(J, dof);CHKERRQ(ierr);
      // CHECK WITH HONG ABOUT THE BELOW CALL
      ierr = MatSetLocalToGlobalMapping(J, ltog, ltog);CHKERRQ(ierr);
      ierr = MatSetLocalToGlobalMappingBlock(J, ltogb, ltogb);CHKERRQ(ierr);
      ierr = MatSetStencil(J, dim, dims, starts, dof);CHKERRQ(ierr);

      bool xperiodic = false;
      // Need to figure out how to tell if y is periodic
      bool yperiodic = true;
      
      for(k=0;k<nz;k++) {
        cout << "----- " << k << " -----" << endl;
        for(j=mesh->ystart; j <= mesh->yend; j++) {
          // cout << "j " << mesh->YGLOBAL(j) << ": ";
          gj = mesh->YGLOBAL(j);
          
          // X range depends on whether there are X boundaries
          int xmin = mesh->xstart;
          if(mesh->firstX())
            xmin = 0; // This processor includes a boundary region
          int xmax = mesh->xend;
          if(mesh->lastX())
            xmax = mesh->ngx-1;
            
          for(i=xmin; i <= xmax; i++) {
            gi = mesh->XGLOBAL(i);
            
            // Check if X and Y are periodic
            yperiodic = mesh->surfaceClosed(i);
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
          // cout << endl;
        }
        // cout << endl;
      }

      ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      // bout_error("stopping");
      ierr = ISLocalToGlobalMappingDestroy(&ltog);CHKERRQ(ierr);
      ierr = ISLocalToGlobalMappingDestroy(&ltogb);CHKERRQ(ierr);
      ierr = TSSetRHSJacobian(ts,J,J,TSDefaultComputeJacobian,this);CHKERRQ(ierr); // remove!!!
    }
  }
#if defined(COLOR)    
  // Create coloring context of J to be used during time stepping 
  ierr = MatGetColoring(J,MATCOLORINGSL,&iscoloring);CHKERRQ(ierr); 
  ierr = MatFDColoringCreate(J,iscoloring,&matfdcoloring);CHKERRQ(ierr);
  ierr = ISColoringDestroy(&iscoloring);CHKERRQ(ierr);
  ierr = MatFDColoringSetFunction(matfdcoloring,(PetscErrorCode (*)(void))solver_f,this);CHKERRQ(ierr);
  ierr = MatFDColoringSetFromOptions(matfdcoloring);CHKERRQ(ierr);
  ierr = TSSetRHSJacobian(ts,J,J,TSDefaultComputeJacobianColor,matfdcoloring);CHKERRQ(ierr);
#endif //OLD

  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);   // enable PETSc runtime options

  // Test TSComputeRHSJacobian()
  // Write J in binary for study - see ~petsc/src/mat/examples/tests/ex124.c
  ierr = PetscOptionsHasName(PETSC_NULL,"-J_write",&J_write);CHKERRQ(ierr);
  if (J_write){ 
    PetscViewer    viewer;
    ierr = PetscPrintf(comm,"\n[%d] Test TSComputeRHSJacobian() ...\n",rank);
    ierr = TSComputeRHSJacobian(ts,simtime,u,&J,&J,&J_structure);CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(comm, "[%d] TSComputeRHSJacobian is done\n",rank);
    ierr = PetscSynchronizedFlush(comm);CHKERRQ(ierr);
    if (J_slowfd){
      ierr = PetscPrintf(comm,"[%d] writing J in binary to data_petsc/Jrhs_dense.dat...\n",rank);CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(comm,"data_petsc/Jrhs_dense.dat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(comm,"[%d] writing J in binary to data_petsc/Jrhs_sparse.dat...\n",rank);CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(comm,"data_petsc/Jrhs_sparse.dat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    }
    ierr = MatView(J,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return(0);
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

PetscErrorCode PetscSolver::run(MonitorFunc mon)
{
  integer steps;
  BoutReal ftime;
  
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

PetscErrorCode PetscSolver::rhs(TS ts, BoutReal t, Vec udata, Vec dudata)
{
  int flag;
  BoutReal *udata_array, *dudata_array;

  PetscFunctionBegin;
#ifdef CHECK
  int msg_point = msg_stack.push("Running RHS: PetscSolver::rhs(%e)", t);
#endif

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
    if(monitor(simtime, iteration, nout)) {
      // User signalled to quit
      
      // Write restart to a different file
      char restartname[512];
      sprintf(restartname, "data/BOUT.final.%d.pdb", MYPE);
      restart.write(restartname);
      
      output.write("Monitor signalled to quit. Returning\n");
      
      PetscFunctionReturn(1);
    }
    */
   
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
void PetscSolver::loop_vars_op(int jx, int jy, BoutReal *udata, int &p, SOLVER_VAR_OP op)
{
  BoutReal **d2d, ***d3d;
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
void PetscSolver::loop_vars(BoutReal *udata, SOLVER_VAR_OP op)
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

void PetscSolver::load_vars(BoutReal *udata)
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
int PetscSolver::save_vars(BoutReal *udata)
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

void PetscSolver::save_derivs(BoutReal *dudata)
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
 * Static functions which can be used for PETSc callbacks
 **************************************************************************/
#undef __FUNCT__  
#define __FUNCT__ "solver_f"
PetscErrorCode solver_f(TS ts, BoutReal t, Vec globalin, Vec globalout, void *f_data)
{
  PetscSolver *s;
  
  PetscFunctionBegin;
  //printf("solver_f(), t %g\n",t);
  s = (PetscSolver*) f_data;
  PetscFunctionReturn(s->rhs(ts, t, globalin, globalout));
}

/*
  FormIFunction = Udot - RHSFunction
*/
#undef __FUNCT__  
#define __FUNCT__ "solver_if"
PetscErrorCode solver_if(TS ts, BoutReal t, Vec globalin,Vec globalindot, Vec globalout, void *f_data)
{
  PetscErrorCode ierr;
  PetscReal      unorm,fnorm;
  
  PetscFunctionBegin;
  ierr = solver_f(ts,t, globalin,globalout, (void *)f_data);CHKERRQ(ierr);
  //ierr = VecNorm(globalin,NORM_INFINITY,&unorm);
  //ierr = VecNorm(globalout,NORM_INFINITY,&fnorm);
  //printf("      solver_if(), t %g, unorm %g, globalout: %g, ",t,unorm,fnorm);

  ierr = VecAYPX(globalout,-1.0,globalindot);CHKERRQ(ierr); // globalout = globalindot + (-1)globalout
  //ierr = VecNorm(globalout,NORM_INFINITY,&fnorm);
  //printf(" udot-rhs: %g\n",fnorm);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "solver_rhsjacobian"
PetscErrorCode solver_rhsjacobian(TS ts,BoutReal t,Vec globalin,Mat *J,Mat *Jpre,MatStructure *str,void *f_data)                          
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  //printf("       solver_jacobian ... a dummy function\n");
  ierr = MatZeroEntries(*Jpre);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(*Jpre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*Jpre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (*J != *Jpre){
    ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  FormIJacobian() - Compute IJacobian = dF/dU + a dF/dUdot  - a dummy matrix used for pc=none
*/
#undef __FUNCT__
#define __FUNCT__ "solver_ijacobian"
PetscErrorCode solver_ijacobian(TS ts,BoutReal t,Vec globalin,Vec globalindot,PetscReal a,Mat *J,Mat *Jpre,MatStructure *str,void *f_data)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = solver_rhsjacobian(ts,t,globalin,J,Jpre,str,(void *)f_data);CHKERRQ(ierr);
  //*Jpre + a
  PetscFunctionReturn(0);
}

//-----------------------------------------

#undef __FUNCT__  
#define __FUNCT__ "PetscSolver::PreUpdate"
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
#define __FUNCT__ "PetscSolver::PostUpdate"
PetscErrorCode PostStep(TS ts) 
{
  PetscFunctionReturn(0);
}

#undef __FUNCT__  

#endif
