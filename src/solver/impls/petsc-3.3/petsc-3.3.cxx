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

#ifdef BOUT_HAS_PETSC_3_3

#include "petsc-3.3.hxx"

#include <petsc-private/tsimpl.h>
#include <petsc.h>

#include <globals.hxx>

#include <stdlib.h>

#include <interpolation.hxx> // Cell interpolation
#include <msg_stack.hxx>
#include <output.hxx>
#include <boutcomm.hxx>

static char help[] = "BOUT++: Uses finite difference methods to solve plasma fluid problems in curvilinear coordinates";

extern PetscErrorCode solver_f(TS ts, BoutReal t, Vec globalin, Vec globalout, void *f_data);
extern PetscErrorCode solver_rhsjacobian(TS ts,BoutReal t,Vec globalin,Mat *J,Mat *Jpre,MatStructure *str,void *f_data);
extern PetscErrorCode solver_if(TS,BoutReal,Vec,Vec,Vec,void*);

extern PetscErrorCode solver_ijacobianfd(TS,PetscReal,Vec,Vec,PetscReal,Mat*,Mat*,MatStructure*,void*);

/// KSP preconditioner PCShell routines for physics preconditioners
extern PetscErrorCode PhysicsPCApply(PC,Vec x,Vec y);
extern PetscErrorCode PhysicsJacobianApply(Mat J, Vec x, Vec y);
extern PetscErrorCode PhysicsSNESApply(SNES,Vec);

PetscSolver::PetscSolver(Options *opts) : Solver(opts) {
  has_constraints = false; // No constraints
  J = 0;
  Jmf = 0;
  matfdcoloring = 0;
  interpolate = PETSC_TRUE;
  initialised = false;
  bout_snes_time = .0;
  
  jacfunc = NULL;

  output_flag = PETSC_FALSE;
}

PetscSolver::~PetscSolver() {
  if(initialised) {
    // Free memory

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

int PetscSolver::init(int NOUT, BoutReal TIMESTEP) {
  TRACE("Initialising PETSc-3.3 solver");

  PetscErrorCode  ierr;
  int             neq;
  int             mudq, mldq, mukeep, mlkeep;
  bool            use_precon, use_jacobian;
  int             precon_dimens;
  BoutReal        precon_tol;
  MPI_Comm        comm = PETSC_COMM_WORLD;
  PetscMPIInt     rank;

  PetscFunctionBegin;
  PetscLogEventRegister("PetscSolver::init",PETSC_VIEWER_CLASSID,&init_event);
  PetscLogEventRegister("loop_vars",PETSC_VIEWER_CLASSID,&loop_event);
  PetscLogEventRegister("solver_f",PETSC_VIEWER_CLASSID,&solver_event);

  /// Call the generic initialisation first
  Solver::init(NOUT, TIMESTEP);

  ierr = PetscLogEventBegin(init_event,0,0,0,0);CHKERRQ(ierr);
  output.write("Initialising PETSc-3.3 solver\n");
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
    output_error.write("\tERROR: MPI_Allreduce failed!\n");
    ierr = PetscLogEventEnd(init_event,0,0,0,0);CHKERRQ(ierr);
    PetscFunctionReturn(1);
  }

  ierr = VecCreate(BoutComm::get(), &u);CHKERRQ(ierr);
  ierr = VecSetSizes(u, local_N, PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);

  ////////// SAVE INITIAL STATE TO PETSc VECTOR ///////////
  // Set pointer to data array in vector u.
  BoutReal *udata;

  ierr = VecGetArray(u,&udata);CHKERRQ(ierr);
  save_vars(udata);
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
#ifdef PETSC_HAS_SUNDIALS
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
#endif
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
  OPTION(options, use_jacobian, false);
  OPTION(options, precon_dimens, 50);
  OPTION(options, precon_tol, 1.0e-4);
  OPTION(options, diagnose,     false);

  BoutReal abstol, reltol;
  options->get("ATOL", abstol, 1.0e-12);
  options->get("RTOL", reltol, 1.0e-5);

#ifdef PETSC_HAS_SUNDIALS
  // Set Sundials tolerances
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
#endif

  // Initial time and timestep. By default just use TIMESTEP
  BoutReal start_timestep;
  OPTION(options, start_timestep, TIMESTEP);
  ierr = TSSetInitialTimeStep(ts,simtime,start_timestep);CHKERRQ(ierr);
  next_output = simtime;

  // Maximum number of steps
  int mxstep;
  OPTION(options, mxstep, 500); // Number of steps between outputs
  mxstep *= NOUT; // Total number of steps
  PetscReal tfinal = simtime + NOUT*TIMESTEP; // Final output time'=
  output.write("\tSet mxstep %d, tfinal %g, simtime %g\n",mxstep,tfinal,simtime);
  ierr = TSSetDuration(ts,mxstep,tfinal);CHKERRQ(ierr);

  // Set the current solution
  ierr = TSSetSolution(ts,u);CHKERRQ(ierr);

  // Create RHSJacobian J
  SNES            snes, psnes;
  KSP             ksp, nksp;
  PC              pc, npc;
  const PCType    pctype;
  const TSType    tstype;
  PetscBool       pcnone=PETSC_TRUE;

  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,PETSC_TRUE);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(PETSC_NULL,"-interpolate",&interpolate,PETSC_NULL);CHKERRQ(ierr);

  // Check for -output_name to see if user specified a "performance"
  // run, if they didn't then use the standard monitor function. TODO:
  // use PetscFList
  ierr = PetscOptionsGetString(PETSC_NULL,"-output_name",this->output_name, sizeof this->output_name,&output_flag);CHKERRQ(ierr);

  // If the output_name is not specified then use the standard monitor function
  if(output_flag) {
    ierr = SNESMonitorSet(snes,PetscSNESMonitor,this,PETSC_NULL);CHKERRQ(ierr);
  } else {
    ierr = TSMonitorSet(ts,PetscMonitor,this,PETSC_NULL);CHKERRQ(ierr);
  }

  ierr = SNESSetTolerances(snes,abstol,reltol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  // Matrix free Jacobian

  if(use_jacobian && (jacfunc != NULL)) {
    // Use a user-supplied Jacobian function
    ierr = MatCreateShell(comm, local_N, local_N, neq, neq, this, &Jmf);
    ierr = MatShellSetOperation(Jmf, MATOP_MULT, (void (*)(void)) PhysicsJacobianApply); CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts, Jmf, Jmf, solver_ijacobian, this); CHKERRQ(ierr);
  }else {
    // Use finite difference approximation
    ierr = MatCreateSNESMF(snes,&Jmf);CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes,Jmf,Jmf,MatMFFDComputeJacobian,this);CHKERRQ(ierr);
  }

  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

  if(use_precon && have_user_precon()) {

    ierr = SNESGetPC(snes,&psnes);CHKERRQ(ierr);
    ierr = SNESGetKSP(psnes,&nksp);CHKERRQ(ierr);
    ierr = KSPGetPC(nksp,&npc);CHKERRQ(ierr);
    ierr = SNESSetType(psnes,SNESSHELL);CHKERRQ(ierr);
    ierr = SNESShellSetSolve(psnes,PhysicsSNESApply);CHKERRQ(ierr);

    // Use a user-supplied preconditioner

    // Tell PETSc we're using a "shell" preconditioner
    ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);

    // Set routine for applying preconditioner
    ierr = PCShellSetApply(pc,PhysicsPCApply);CHKERRQ(ierr);

    // Set context to this solver object
    ierr = PCShellSetContext(pc,this);CHKERRQ(ierr);

    // Set name of preconditioner
    ierr = PCShellSetName(pc,"PhysicsPreconditioner");CHKERRQ(ierr);

    // Need a callback for IJacobian to get shift 'alpha'
    ierr = TSSetIJacobian(ts, Jmf, Jmf, solver_ijacobian, this);

    // Use right preconditioner
    ierr = KSPSetPCSide(ksp, PC_RIGHT);CHKERRQ(ierr);

  }else {
    // Default to no preconditioner
    ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
  }
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);   // enable PETSc runtime options

  ierr = PCGetType(pc,&pctype);CHKERRQ(ierr);
  ierr = TSGetType(ts,&tstype);CHKERRQ(ierr);
  output.write("\tTS type %s, PC type %s\n",tstype,pctype);

  ierr = PetscObjectTypeCompare((PetscObject)pc,PCNONE,&pcnone);CHKERRQ(ierr);
  if (pcnone) {
    ierr = PetscLogEventEnd(init_event,0,0,0,0);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = PetscObjectTypeCompare((PetscObject)pc,PCSHELL,&pcnone);CHKERRQ(ierr);
  if (pcnone) {
    ierr = PetscLogEventEnd(init_event,0,0,0,0);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

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

    //PetscInt MXSUB = mesh->xend - mesh->xstart + 1;
    //PetscInt MYSUB = mesh->yend - mesh->ystart + 1;

    PetscInt nx = mesh->xend;//MXSUB;
    PetscInt ny = mesh->yend;//MYSUB;

    /* number of z points */
    PetscInt nz  = mesh->LocalNz;

    /* number of degrees (variables) at each grid point */
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

      ierr = SNESSetJacobian(snes,J,J,SNESDefaultComputeJacobian,PETSC_NULL);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"SNESComputeJacobian J by slow fd...\n");CHKERRQ(ierr);

      ierr = TSComputeRHSJacobian(ts,simtime,u,&J,&J,&flg);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"compute J by slow fd is done.\n");CHKERRQ(ierr);
      //ierr = MatView(J,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    } else { // get sparse pattern of the Jacobian
      throw BoutException("Path followed in PETSc solver not yet implemented in general "
                          "-- experimental hard coded values here. Sorry\n");

      ierr = PetscPrintf(PETSC_COMM_WORLD,"get sparse pattern of the Jacobian...\n");CHKERRQ(ierr);

      if(n2Dvars() != 0) throw BoutException("PETSc solver can't handle 2D variables yet. Sorry\n");

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

      output.write(" dof %d,dim %d: %d %d %d\n",dof,dim,dims[0],dims[1],dims[2]);
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
  ierr = SNESSetJacobian(snes,J,J,SNESDefaultComputeJacobianColor,matfdcoloring);CHKERRQ(ierr);

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

  ierr = PetscLogEventEnd(init_event,0,0,0,0);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

PetscErrorCode PetscSolver::run() {
  PetscErrorCode ierr;
  //integer steps;
  BoutReal ftime;
  FILE *fp = NULL;

  // Set when the next call to monitor is desired
  // next_output = simtime + tstep;

  PetscFunctionBegin;

  if(this->output_flag) {
    prev_linear_its = 0;
    bout_snes_time = MPI_Wtime();
  }

  ierr = TSSolve(ts,u,&ftime);CHKERRQ(ierr);

  // Gawd, everything is a hack
  if(this->output_flag) {
    ierr = PetscFOpen(PETSC_COMM_WORLD, this->output_name, "w", &fp);CHKERRQ(ierr);
    ierr = PetscFPrintf(PETSC_COMM_WORLD, fp, "SNES Iteration, KSP Iterations, Wall Time, Norm\n");CHKERRQ(ierr);
    for(int i =0;i < snes_list.size();i++) {
      ierr = PetscFPrintf(PETSC_COMM_WORLD, fp, "%i, %i, %e, %e\n", snes_list[i].it, snes_list[i].linear_its, snes_list[i].time, snes_list[i].norm);CHKERRQ(ierr);
    }
    ierr = PetscFClose(PETSC_COMM_WORLD, fp);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/**************************************************************************
 * RHS function
 **************************************************************************/

PetscErrorCode PetscSolver::rhs(TS ts, BoutReal t, Vec udata, Vec dudata) {
  TRACE("Running RHS: Petsc33Solver::rhs(%e)", t);

  BoutReal *udata_array, *dudata_array;

  PetscFunctionBegin;

  // Load state from PETSc
  VecGetArray(udata, &udata_array);
  load_vars(udata_array);
  VecRestoreArray(udata, &udata_array);

  // Call RHS function
  int flag = run_rhs(t);

  // Save derivatives to PETSc
  VecGetArray(dudata, &dudata_array);
  save_derivs(dudata_array);
  VecRestoreArray(dudata, &dudata_array);

  PetscFunctionReturn(0);
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

PetscErrorCode PetscSolver::pre(PC pc, Vec x, Vec y) {
  TRACE("Petsc33Solver::pre()");

  BoutReal *data;

  if(diagnose)
    output << "Preconditioning" << endl;

  // Load state
  VecGetArray(state, &data);
  load_vars(data);
  VecRestoreArray(state, &data);

  // Load vector to be inverted into F_vars
  VecGetArray(x, &data);
  load_derivs(data);
  VecRestoreArray(x, &data);

  // Call the preconditioner
  run_precon(ts_time, 1./shift, 0.0);

  // Save the solution from time derivatives
  VecGetArray(y, &data);
  save_derivs(data);
  VecRestoreArray(y, &data);

  // Petsc's definition of Jacobian differs by a factor from Sundials'
  PetscErrorCode ierr = VecScale(y, shift); CHKERRQ(ierr);

   return 0;
 }

/**************************************************************************
 * User-supplied Jacobian function J(state) * x = y
 **************************************************************************/

PetscErrorCode PetscSolver::jac(Vec x, Vec y) {
  TRACE("Petsc33Solver::jac()");

  BoutReal *data;

  if(diagnose)
    output << "Jacobian evaluation\n";

  // Load state
  VecGetArray(state, &data);
  load_vars(data);
  VecRestoreArray(state, &data);

  // Load vector to be operated on into F_vars
  VecGetArray(x, &data);
  load_derivs(data);
  VecRestoreArray(x, &data);

  // Call the Jacobian function
  (*jacfunc)(ts_time);

  // Save the solution from time derivatives
  VecGetArray(y, &data);
  save_derivs(data);
  VecRestoreArray(y, &data);

  // y = a * x - y
  int ierr = VecAXPBY(y, shift, -1.0, x);

  return 0;
}

/**************************************************************************
 * Static functions which can be used for PETSc callbacks
 **************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "solver_f"
PetscErrorCode solver_f(TS ts, BoutReal t, Vec globalin, Vec globalout, void *f_data) {
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
PetscErrorCode solver_if(TS ts, BoutReal t, Vec globalin,Vec globalindot, Vec globalout, void *f_data) {
  PetscErrorCode ierr;
  //PetscReal      unorm,fnorm;

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
PetscErrorCode solver_rhsjacobian(TS ts,BoutReal t,Vec globalin,Mat *J,Mat *Jpre,MatStructure *str,void *f_data) {
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
PetscErrorCode solver_ijacobian(TS ts,BoutReal t,Vec globalin,Vec globalindot,PetscReal a,Mat *J,Mat *Jpre,MatStructure *str,void *f_data) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = solver_rhsjacobian(ts,t,globalin,J,Jpre,str,(void *)f_data);CHKERRQ(ierr);

  ////// Save data for preconditioner
  PetscSolver *solver = (PetscSolver*) f_data;

  if(solver->diagnose)
    output << "Saving state, t = " << t << ", a = " << a << endl;

  solver->shift = a; // Save the shift 'a'
  solver->state = globalin;  // Save system state
  solver->ts_time = t;

  PetscFunctionReturn(0);
}

/*
  solver_ijacobianfd - Compute IJacobian = dF/dU + a dF/dUdot  using finite deference - not implemented yet
*/
#undef __FUNCT__
#define __FUNCT__ "solver_ijacobianfd"
PetscErrorCode solver_ijacobianfd(TS ts,BoutReal t,Vec globalin,Vec globalindot,PetscReal a,Mat *J,Mat *Jpre,MatStructure *str,void *f_data) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = solver_rhsjacobian(ts,t,globalin,J,Jpre,str,(void *)f_data);CHKERRQ(ierr);
  //*Jpre + a
  PetscFunctionReturn(0);
}
//-----------------------------------------

#undef __FUNCT__
#define __FUNCT__ "PhysicsSNESApply"
PetscErrorCode PhysicsSNESApply(SNES snes, Vec x) {
  PetscErrorCode ierr;
  Vec F,Fout;
  PetscReal fnorm = 0., foutnorm = 0., dot=0.;
  KSP ksp;
  PC pc;
  Mat A,B;
  MatStructure diff = DIFFERENT_NONZERO_PATTERN;

  PetscFunctionBegin;
  ierr = SNESGetJacobian(snes, &A, &B, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);
  ierr = SNESComputeJacobian(snes, x, &A, &B, &diff);CHKERRQ(ierr);
  ierr = SNESGetKSP(snes, &ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
  ierr = SNESGetFunction(snes,&F,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = SNESComputeFunction(snes, x, F);CHKERRQ(ierr);
  ierr = SNESGetSolutionUpdate(snes, &Fout);CHKERRQ(ierr);

  ierr = PCApply(pc,F,Fout);CHKERRQ(ierr);
  ierr = VecNorm(Fout, NORM_2, &foutnorm);CHKERRQ(ierr);
  ierr = VecAXPY(x, -1., Fout);CHKERRQ(ierr);
  ierr = SNESComputeFunction(snes, x, F);CHKERRQ(ierr);
  ierr = VecNorm(F,NORM_2,&fnorm);CHKERRQ(ierr);
  ierr = VecDot(F,Fout,&dot);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, " (Debug) function norm: %g, P(f) norm %g, F \\cdot Fout %g  ", fnorm, foutnorm, dot);CHKERRQ(ierr);
  ierr = SNESSetFunctionNorm(snes, fnorm);CHKERRQ(ierr);
  ierr = SNESMonitor(snes,0,fnorm);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PhysicsPCApply"
PetscErrorCode PhysicsPCApply(PC pc,Vec x,Vec y) {
  int ierr;

  // Get the context
  PetscSolver *s;
  ierr = PCShellGetContext(pc,(void**)&s);CHKERRQ(ierr);

  PetscFunctionReturn(s->pre(pc, x, y));
}

#undef __FUNCT__
#define __FUNCT__ "PhysicsJacobianApply"
PetscErrorCode PhysicsJacobianApply(Mat J, Vec x, Vec y) {
  // Get the context
  PetscSolver *s;
  int ierr = MatShellGetContext(J, (void**)&s); CHKERRQ(ierr);
  PetscFunctionReturn(s->jac(x, y));
}

#undef __FUNCT__
#define __FUNCT__ "PetscMonitor"
PetscErrorCode PetscMonitor(TS ts,PetscInt step,PetscReal t,Vec X,void *ctx) {
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

    if (s->call_monitors(simtime,i++,s->nout)) {
      PetscFunctionReturn(1);
    }

    // Reset counters
    s->rhs_ncalls = 0;

    s->next_output += s->tstep;
    simtime = s->next_output;
  }

  /* Done with vector, so destroy it */
  ierr = VecDestroy(&interpolatedX);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscSNESMonitor"
PetscErrorCode PetscSNESMonitor(SNES snes, PetscInt its, PetscReal norm, void *ctx) {
  PetscErrorCode ierr;
  PetscInt linear_its=0;
  BoutReal tmp = .0;
  snes_info row;
  PetscSolver *s = (PetscSolver*)ctx;

  PetscFunctionBegin;

  if(!its) s->prev_linear_its = 0;
  ierr = SNESGetLinearSolveIterations(snes, &linear_its);CHKERRQ(ierr);
  tmp = MPI_Wtime();

  row.it = its;
  s->prev_linear_its = row.linear_its = linear_its-s->prev_linear_its;
  row.time = tmp-s->bout_snes_time;
  row.norm = norm;

  s->snes_list.push_back(row);

  PetscFunctionReturn(0);
}

#endif // BOUT_HAS_PETSC_3_3
