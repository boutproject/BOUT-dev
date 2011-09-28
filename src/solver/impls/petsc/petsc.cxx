/**************************************************************************
 * Interface to PETSc solver
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

#include <private/tsimpl.h>

#include <globals.hxx>

#include <stdlib.h>

#include <interpolation.hxx> // Cell interpolation


extern PetscErrorCode solver_f(TS ts, BoutReal t, Vec globalin, Vec globalout, void *f_data);
extern PetscErrorCode solver_rhsjacobian(TS ts,BoutReal t,Vec globalin,Mat *J,Mat *Jpre,MatStructure *str,void *f_data);
extern PetscErrorCode solver_if(TS,BoutReal,Vec,Vec,Vec,void*);
extern PetscErrorCode solver_ijacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode solver_ijacobianfd(TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);

PetscSolver::PetscSolver()
{
  has_constraints = false; // No constraints
  this->J = 0;
  this->Jmf = 0;
  this->matfdcoloring = 0;
  this->interpolate = PETSC_TRUE;
  PetscLogEventRegister("loop_vars",PETSC_VIEWER_CLASSID,&loop_event);
  PetscLogEventRegister("solver_f",PETSC_VIEWER_CLASSID,&solver_event);
  PetscLogEventRegister("PetscSolver::init",PETSC_VIEWER_CLASSID,&init_event);
}

PetscSolver::~PetscSolver()
{

  if(initialised) {
    // Free CVODE memory

    VecDestroy(&u);
    if (J) {MatDestroy(&J);}
    if (Jmf) {MatDestroy(&Jmf);}
    if (matfdcoloring) {MatFDColoringDestroy(&matfdcoloring);}
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
  PetscMPIInt     rank;

#ifdef CHECK
  int msg_point = msg_stack.push("Initialising PETSc solver");
#endif

  PetscLogEventBegin(init_event,0,0,0,0);
  /// Call the generic initialisation first
  Solver::init(f, argc, argv, restarting, NOUT, TIMESTEP);

  output.write("Initialising PETSc solver\n");
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);

  // Save NOUT and TIMESTEP for use later
  nout = NOUT;
  tstep = TIMESTEP;

  // Set the rhs solver function
  // func = f;

  PetscInt n2d = n2Dvars();       // Number of 2D variables
  PetscInt n3d = n3Dvars();       // Number of 3D variables
  PetscInt local_N = getLocalN(); // Number of evolving variables on this processor

  /********** Get total problem size **********/
  if(MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    output.write("\tERROR: MPI_Allreduce failed!\n");
    return 1;
  }

  ierr = PetscPrintf(comm,"\t3d fields = %d, 2d fields = %d neq=%d\n",n3d, n2d, neq);CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(comm, "\t[%d] local_N %d\n",rank,local_N);CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(comm);CHKERRQ(ierr);

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
  PetscReal norm;
  ierr = VecNorm(u,NORM_1,&norm);CHKERRQ(ierr);
  if (!rank){ierr = PetscPrintf(PETSC_COMM_SELF,"initial |u| = %g\n",norm);}
  //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  PetscBool       J_load;
  MatStructure    J_structure;
  char            load_file[PETSC_MAX_PATH_LEN];  /* jacobian input file name */
  PetscBool       J_write=PETSC_FALSE,J_slowfd=PETSC_FALSE;
  ISColoring      iscoloring;

  // Create timestepper
  ierr = TSCreate(BoutComm::get(),&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  ierr = TSSetApplicationContext(ts, this);CHKERRQ(ierr);

  // Set user provided RHSFunction
  // Need to duplicate the solution vector for the residual
  Vec rhs_vec;
  ierr = VecDuplicate(u,&rhs_vec);
  ierr = TSSetRHSFunction(ts,rhs_vec,solver_f,this);CHKERRQ(ierr);
  ierr = VecDestroy(&rhs_vec);

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
  //printf("abstol %g, reltol %g\n",abstol,reltol); why reltol=1.e-7?

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
  next_output = simtime;

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
  const TSType    tstype;
  PetscBool       pcnone=PETSC_TRUE;

  ierr = TSSetExactFinalTime(ts,PETSC_TRUE);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(PETSC_NULL,"-interpolate",&interpolate,PETSC_NULL);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,PetscMonitor,this,PETSC_NULL);CHKERRQ(ierr);

  // Default to matrix-free
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESSetTolerances(snes,abstol,reltol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  ierr = MatCreateSNESMF(snes,&Jmf);CHKERRQ(ierr);
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,Jmf,Jmf,MatMFFDComputeJacobian,this);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

  // Default to no preconditioner
  ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);   // enable PETSc runtime options

  ierr = PCGetType(pc,&pctype);CHKERRQ(ierr);
  ierr = TSGetType(ts,&tstype);CHKERRQ(ierr);
  output.write("\tTS type %s, PC type %s\n",tstype,pctype);

  ierr = PetscTypeCompare((PetscObject)pc,PCNONE,&pcnone);CHKERRQ(ierr);
  if (pcnone) return(0);

  // Create Jacobian matrix to be used by preconditioner
  ierr = PetscPrintf(PETSC_COMM_WORLD," Get Jacobian matrix at simtime %g\n",simtime);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-J_load",load_file,PETSC_MAX_PATH_LEN-1,&J_load);CHKERRQ(ierr);
  if(J_load) {
    PetscViewer     fd;
    ierr = PetscPrintf(PETSC_COMM_WORLD," Load Jmat ...local_N %d, neq %d\n",local_N,neq);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(comm,load_file,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
    //ierr = MatSetType(J, MATBAIJ);CHKERRQ(ierr);
    ierr = MatSetSizes(J,local_N,local_N,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);CHKERRQ(ierr);
    ierr = MatLoad(J, fd);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  } else { // create Jacobian matrix

    PetscInt MXSUB = mesh->xend - mesh->xstart + 1;
    PetscInt MYSUB = mesh->yend - mesh->ystart + 1;

    PetscInt nx = mesh->xend;//MXSUB;
    PetscInt ny = mesh->yend;//MYSUB;

    /* number of z points (need to subtract one because of historical reasons that MZ has an extra point) */
    PetscInt nz  = mesh->ngz - 1;

    /* number of degrees (variables) at each grid point */
    if(n2Dvars() != 0) bout_error("PETSc solver can't handle 2D variables yet. Sorry\n");

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
    ierr = MatSetSizes(J,local_N, local_N, neq,neq);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);CHKERRQ(ierr);

    // Get nonzero pattern of J - color_none !!!
    prealloc = cols*dof*dof;
    ierr = MatSeqAIJSetPreallocation(J,prealloc,PETSC_NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(J,prealloc,PETSC_NULL,prealloc,PETSC_NULL);CHKERRQ(ierr);

    prealloc = cols; // why nonzeros=295900, allocated nonzeros=2816000/12800000 (*dof*dof), number of mallocs used during MatSetValues calls =256?
    ierr = MatSeqBAIJSetPreallocation(J,dof,prealloc,PETSC_NULL);CHKERRQ(ierr);
    ierr = MatMPIBAIJSetPreallocation(J,dof,prealloc,PETSC_NULL,prealloc,PETSC_NULL);CHKERRQ(ierr);

    ierr = PetscOptionsHasName(PETSC_NULL,"-J_slowfd",&J_slowfd);CHKERRQ(ierr);
    if (J_slowfd) { // create Jacobian matrix by slow fd
      MatStructure flg;

      ierr = TSSetRHSJacobian(ts,J,J,TSDefaultComputeJacobian,PETSC_NULL);CHKERRQ(ierr); //cannot set TSSetRHSJacobian(ts,Jmf,J,...)?
      ierr = PetscPrintf(PETSC_COMM_WORLD,"SNESComputeJacobian J by slow fd...\n");CHKERRQ(ierr);

      ierr = TSComputeRHSJacobian(ts,simtime,u,&J,&J,&flg);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"compute J by slow fd is done.\n");CHKERRQ(ierr);
      //ierr = MatView(J,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    } else { // get sparse pattern of the Jacobian
      ierr = PetscPrintf(PETSC_COMM_WORLD,"get sparse pattern of the Jacobian...\n");CHKERRQ(ierr);

      ISLocalToGlobalMapping ltog, ltogb;
      PetscInt i, j, k, d, s;
      PetscInt gi, gj;

      MatStencil stencil[cols];

      PetscScalar one[preallocblock];
      for (i = 0; i < preallocblock; i++) one[i] = 1.0;

      // Change this block for parallel. Currently this is just the identity
      // map since advect1d has no branch cuts (and we are only testing
      // single processor now)
      PetscInt ltog_array[n];
      for (i = 0; i < n; i++) ltog_array[i] = i;

      // Also change this for parallel. This defines the 'global stencil'
      // where starts are the starting index in each dimension and dims
      // are the size
      PetscInt starts[3], dims[3];
      starts[0] = starts[1] = starts[2] = 0;

      // Only for advect1d, need to figure this out
      nx = 5;
      ny = 128;

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

      printf(" dof %d,dim %d: %d %d %d\n",dof,dim,dims[0],dims[1],dims[2]);
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
            //printf("stencil: %d, %d, %d; -- %d %d %d %d\n",gi,gj,k,stencil[0].i,stencil[0].j,stencil[0].k,stencil[0].c);

          }
          // cout << endl;
        }
        // cout << endl;
      }

      ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

      //J(2585:n-1,:) has no entry !!!
      // ierr = PetscPrintf(PETSC_COMM_WORLD,"Sparse pattern:\n");
      // ierr = MatView(J,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      // bout_error("stopping");
      ierr = ISLocalToGlobalMappingDestroy(&ltog);CHKERRQ(ierr);
      ierr = ISLocalToGlobalMappingDestroy(&ltogb);CHKERRQ(ierr);
    }
  }

  // Create coloring context of J to be used during time stepping
  ierr = PetscPrintf(PETSC_COMM_WORLD," Create coloring ...\n");
  ierr = MatGetColoring(J,MATCOLORINGSL,&iscoloring);CHKERRQ(ierr);
  ierr = MatFDColoringCreate(J,iscoloring,&matfdcoloring);CHKERRQ(ierr);
  ierr = MatFDColoringSetFromOptions(matfdcoloring);CHKERRQ(ierr);
  ierr = ISColoringDestroy(&iscoloring);CHKERRQ(ierr);
  ierr = MatFDColoringSetFunction(matfdcoloring,(PetscErrorCode (*)(void))solver_f,this);CHKERRQ(ierr);
  ierr = TSSetRHSJacobian(ts,J,J,TSDefaultComputeJacobianColor,matfdcoloring);CHKERRQ(ierr);

  // Write J in binary for study - see ~petsc/src/mat/examples/tests/ex124.c
  ierr = PetscOptionsHasName(PETSC_NULL,"-J_write",&J_write);CHKERRQ(ierr);
  if (J_write){
    PetscViewer    viewer;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n[%d] Test TSComputeRHSJacobian() ...\n",rank);
    ierr = TSComputeRHSJacobian(ts,simtime,u,&J,&J,&J_structure);CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(comm, "[%d] TSComputeRHSJacobian is done\n",rank);
    ierr = PetscSynchronizedFlush(comm);CHKERRQ(ierr);
    if (J_slowfd) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"[%d] writing J in binary to data/Jrhs_dense.dat...\n",rank);CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(comm,"data/Jrhs_dense.dat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"[%d] writing J in binary to data/Jrhs_sparse.dat...\n",rank);CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(comm,"data/Jrhs_sparse.dat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    }
    ierr = MatView(J,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  PetscLogEventEnd(init_event,0,0,0,0);
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
  // next_output = simtime + tstep;
  monitor = mon; // Store the monitor function pointer

  PetscFunctionBegin;
  PetscFunctionReturn(TSSolve(ts,u,&ftime));
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

  PetscLogEventBegin(loop_event,0,0,0,0);
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
  PetscLogEventEnd(loop_event,0,0,0,0);
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
  s = (PetscSolver*) f_data;
  PetscLogEventBegin(s->solver_event,0,0,0,0);
  s->rhs(ts, t, globalin, globalout);
  PetscLogEventEnd(s->solver_event,0,0,0,0);
  PetscFunctionReturn(0);
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
  //ierr = MatZeroEntries(*Jpre);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(*Jpre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*Jpre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (*J != *Jpre){
    ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  solver_ijacobian - Compute IJacobian = dF/dU + a dF/dUdot  - a dummy matrix used for pc=none
*/
#undef __FUNCT__
#define __FUNCT__ "solver_ijacobian"
PetscErrorCode solver_ijacobian(TS ts,BoutReal t,Vec globalin,Vec globalindot,PetscReal a,Mat *J,Mat *Jpre,MatStructure *str,void *f_data)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  //printf("     solver_ijacobian...\n");
  ierr = solver_rhsjacobian(ts,t,globalin,J,Jpre,str,(void *)f_data);CHKERRQ(ierr);
  //*Jpre + a
  PetscFunctionReturn(0);
}

/*
  solver_ijacobianfd - Compute IJacobian = dF/dU + a dF/dUdot  using finite deference - not implemented yet
*/
#undef __FUNCT__
#define __FUNCT__ "solver_ijacobianfd"
PetscErrorCode solver_ijacobianfd(TS ts,BoutReal t,Vec globalin,Vec globalindot,PetscReal a,Mat *J,Mat *Jpre,MatStructure *str,void *f_data)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = solver_rhsjacobian(ts,t,globalin,J,Jpre,str,(void *)f_data);CHKERRQ(ierr);
  //*Jpre + a
  PetscFunctionReturn(0);
}
//-----------------------------------------

#undef __FUNCT__
#define __FUNCT__ "PetscMonitor"
PetscErrorCode PetscMonitor(TS ts,PetscInt step,PetscReal t,Vec X,void *ctx)
{
  PetscErrorCode ierr;
  PetscSolver *s = (PetscSolver *)ctx;
  PetscReal tfinal, dt;
  Vec interpolatedX;
  const PetscScalar *x;
  static int i = 0;

  PetscFunctionBegin;
  ierr = TSGetTimeStep(ts, &dt);CHKERRQ(ierr);
  ierr = TSGetDuration(ts, PETSC_NULL, &tfinal);CHKERRQ(ierr);

  /* Duplicate the solution vector X into a work vector */
  ierr = VecDuplicate(X,&interpolatedX);CHKERRQ(ierr);
  while (s->next_output <= t && s->next_output <= tfinal) {
    if (s->interpolate) ierr = TSInterpolate(ts,s->next_output,interpolatedX);CHKERRQ(ierr);

    /* Place the interpolated values into the global variables */
    ierr = VecGetArrayRead(interpolatedX,&x);CHKERRQ(ierr);
    s->load_vars((BoutReal *)x);
    ierr = VecRestoreArrayRead(interpolatedX,&x);CHKERRQ(ierr);

    if (s->monitor(simtime,i++,s->nout)) {
      s->restart.write("%s/BOUT.final.%d.%s", s->restartdir.c_str(), s->MYPE, s->restartext.c_str());

      output.write("Monitor signalled to quit. Returning\n");
    }

    s->next_output += s->tstep;
    simtime = s->next_output;
  }

  /* Done with vector, so destroy it */
  ierr = VecDestroy(&interpolatedX);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#endif
