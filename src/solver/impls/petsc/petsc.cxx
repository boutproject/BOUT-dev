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

#include "bout/build_defines.hxx"

#include "petsc.hxx"

#if BOUT_HAS_PETSC

#include <petsc.h>

#include <bout/boutcomm.hxx>

#include <cstdlib>

#include <bout/interpolation.hxx> // Cell interpolation
#include <bout/msg_stack.hxx>
#include <bout/output.hxx>

extern PetscErrorCode solver_f(TS ts, BoutReal t, Vec globalin, Vec globalout,
                               void* f_data);

#if PETSC_VERSION_GE(3, 5, 0)
extern PetscErrorCode solver_rhsjacobian(TS ts, BoutReal t, Vec globalin, Mat J, Mat Jpre,
                                         void* f_data);
extern PetscErrorCode solver_ijacobianfd(TS, PetscReal, Vec, Vec, PetscReal, Mat, Mat,
                                         void*);
#else
extern PetscErrorCode solver_rhsjacobian(TS ts, BoutReal t, Vec globalin, Mat* J,
                                         Mat* Jpre, MatStructure* str, void* f_data);
extern PetscErrorCode solver_ijacobianfd(TS, PetscReal, Vec, Vec, PetscReal, Mat*, Mat*,
                                         MatStructure*, void*);
#endif

extern PetscErrorCode solver_if(TS, BoutReal, Vec, Vec, Vec, void*);

/// KSP preconditioner PCShell routines for physics preconditioners
extern PetscErrorCode PhysicsPCApply(PC, Vec x, Vec y);
extern PetscErrorCode PhysicsJacobianApply(Mat J, Vec x, Vec y);
extern PetscErrorCode PhysicsSNESApply(SNES, Vec);

PetscSolver::PetscSolver(Options* opts)
    : Solver(opts),
      diagnose(
          (*options)["diagnose"].doc("Enable some diagnostic output").withDefault(false)),
      adaptive(
          (*options)["adaptive"].doc("Use adaptive timestepping").withDefault(false)),
      use_precon((*options)["use_precon"]
                     .doc("Use user-supplied preconditioning function")
                     .withDefault(false)),
      use_jacobian((*options)["use_jacobian"]
                       .doc("Use user-supplied Jacobian function")
                       .withDefault(false)),
      abstol((*options)["atol"].doc("Absolute tolerance").withDefault(1.0e-12)),
      reltol((*options)["rtol"].doc("Relative tolerance").withDefault(1.0e-5)),
      adams_moulton((*options)["adams_moulton"]
                        .doc("Use Adams-Moulton implicit multistep method instead of BDF "
                             "(requires PETSc to have been built with SUNDIALS)")
                        .withDefault(false)),
      start_timestep((*options)["start_timestep"]
                         .doc("Initial internal timestep (defaults to output timestep)")
                         .withDefault(getOutputTimestep())),
      mxstep(
          (*options)["mxstep"].doc("Number of steps between outputs").withDefault(500)) {}

PetscSolver::~PetscSolver() {
  VecDestroy(&u);
  MatDestroy(&J);
  MatDestroy(&Jmf);
  MatFDColoringDestroy(&matfdcoloring);
  TSDestroy(&ts);
}

/**************************************************************************
 * Initialise
 **************************************************************************/

int PetscSolver::init() {
  PetscErrorCode ierr;
  int neq;
  MPI_Comm comm = PETSC_COMM_WORLD;
  PetscMPIInt rank;

  TRACE("Initialising PETSc-dev solver");

  PetscFunctionBegin;
  PetscLogEventRegister("PetscSolver::init", PETSC_VIEWER_CLASSID, &init_event);
  PetscLogEventRegister("loop_vars", PETSC_VIEWER_CLASSID, &loop_event);
  PetscLogEventRegister("solver_f", PETSC_VIEWER_CLASSID, &solver_event);

  Solver::init();

  ierr = PetscLogEventBegin(init_event, 0, 0, 0, 0);
  CHKERRQ(ierr);
  output.write("Initialising PETSc-dev solver\n");
  ierr = bout::globals::mpi->MPI_Comm_rank(comm, &rank);
  CHKERRQ(ierr);

  PetscInt local_N = getLocalN(); // Number of evolving variables on this processor

  /********** Get total problem size **********/
  if (bout::globals::mpi->MPI_Allreduce(&local_N, &neq, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    output_error.write("\tERROR: MPI_Allreduce failed!\n");
    ierr = PetscLogEventEnd(init_event, 0, 0, 0, 0);
    CHKERRQ(ierr);
    PetscFunctionReturn(1);
  }

  ierr = VecCreate(BoutComm::get(), &u);
  CHKERRQ(ierr);
  ierr = VecSetSizes(u, local_N, PETSC_DECIDE);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);
  CHKERRQ(ierr);

  ////////// SAVE INITIAL STATE TO PETSc VECTOR ///////////
  // Set pointer to data array in vector u.
  BoutReal* udata;

  ierr = VecGetArray(u, &udata);
  CHKERRQ(ierr);
  save_vars(udata);
  ierr = VecRestoreArray(u, &udata);
  CHKERRQ(ierr);

  PetscReal norm;
  ierr = VecNorm(u, NORM_1, &norm);
  CHKERRQ(ierr);
  output_info << "initial |u| = " << norm << "\n";

  PetscBool J_load;
  char load_file[PETSC_MAX_PATH_LEN]; /* jacobian input file name */
  PetscBool J_write = PETSC_FALSE, J_slowfd = PETSC_FALSE;

  // Create timestepper
  ierr = TSCreate(BoutComm::get(), &ts);
  CHKERRQ(ierr);
  ierr = TSSetProblemType(ts, TS_NONLINEAR);
  CHKERRQ(ierr);
#if PETSC_HAS_SUNDIALS
  ierr = TSSetType(ts, TSSUNDIALS);
  CHKERRQ(ierr);
#else
  ierr = TSSetType(ts, TSRK);
  CHKERRQ(ierr);
#endif
  ierr = TSSetApplicationContext(ts, this);
  CHKERRQ(ierr);

  // Set user provided RHSFunction
  // Need to duplicate the solution vector for the residual
  Vec rhs_vec;
  ierr = VecDuplicate(u, &rhs_vec);
  CHKERRQ(ierr);
  ierr = TSSetRHSFunction(ts, rhs_vec, solver_f, this);
  CHKERRQ(ierr);
  ierr = VecDestroy(&rhs_vec);
  CHKERRQ(ierr);

  // Set up adaptive time-stepping
  TSAdapt adapt;
  ierr = TSGetAdapt(ts, &adapt);
  CHKERRQ(ierr);
  if (adaptive) {
    ierr = TSAdaptSetType(adapt, TSADAPTBASIC);
    CHKERRQ(ierr);
  } else {
    ierr = TSAdaptSetType(adapt, TSADAPTNONE);
    CHKERRQ(ierr);
  }

  // Set default absolute/relative tolerances
  ierr = TSSetTolerances(ts, abstol, nullptr, reltol, nullptr);
  CHKERRQ(ierr);

#if PETSC_VERSION_LT(3, 5, 0)
  ierr = TSRKSetTolerance(ts, reltol);
  CHKERRQ(ierr);
#endif

#if PETSC_HAS_SUNDIALS
  // Set Sundials tolerances
  ierr = TSSundialsSetTolerance(ts, abstol, reltol);
  CHKERRQ(ierr);

  // Select Sundials Adams-Moulton or BDF method
  if (adams_moulton) {
    output.write("\tUsing Adams-Moulton implicit multistep method\n");
    ierr = TSSundialsSetType(ts, SUNDIALS_ADAMS);
    CHKERRQ(ierr);
  } else {
    output.write("\tUsing BDF method\n");
    ierr = TSSundialsSetType(ts, SUNDIALS_BDF);
    CHKERRQ(ierr);
  }
#endif

  // Initial time and timestep
  ierr = TSSetTime(ts, simtime);
  CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts, start_timestep);
  CHKERRQ(ierr);
  next_output = simtime;

  // Total number of steps
  PetscInt total_steps = mxstep * getNumberOutputSteps();
  // Final output time
  PetscReal tfinal = simtime + (getNumberOutputSteps() * getOutputTimestep());
  output.write("\tSet total_steps {:d}, tfinal {:g}, simtime {:g}\n", total_steps, tfinal,
               simtime);

#if PETSC_VERSION_GE(3, 8, 0)
  ierr = TSSetMaxSteps(ts, total_steps);
  CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts, tfinal);
  CHKERRQ(ierr);
#else
  ierr = TSSetDuration(ts, total_steps, tfinal);
  CHKERRQ(ierr);
#endif

  // Set the current solution
  ierr = TSSetSolution(ts, u);
  CHKERRQ(ierr);

  // Create RHSJacobian J
  SNES snes, psnes;
  KSP ksp, nksp;
  PC pc, npc;
  PCType pctype;
  TSType tstype;
  PetscBool pcnone = PETSC_TRUE;

  ierr = TSGetSNES(ts, &snes);
  CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_INTERPOLATE);
  CHKERRQ(ierr);

#if PETSC_VERSION_GE(3, 7, 0)
  ierr = PetscOptionsGetBool(nullptr, nullptr, "-interpolate", &interpolate, nullptr);
  CHKERRQ(ierr);
#else
  ierr = PetscOptionsGetBool(nullptr, "-interpolate", &interpolate, nullptr);
  CHKERRQ(ierr);
#endif

  // Check for -output_name to see if user specified a "performance"
  // run, if they didn't then use the standard monitor function. TODO:
  // use PetscFList
#if PETSC_VERSION_GE(3, 7, 0)
  ierr = PetscOptionsGetString(nullptr, nullptr, "-output_name", this->output_name,
                               sizeof this->output_name, &output_flag);
  CHKERRQ(ierr);
#else
  ierr = PetscOptionsGetString(nullptr, "-output_name", this->output_name,
                               sizeof this->output_name, &output_flag);
  CHKERRQ(ierr);
#endif

  // If the output_name is not specified then use the standard monitor function
  if (output_flag != 0U) {
    ierr = SNESMonitorSet(snes, PetscSNESMonitor, this, nullptr);
    CHKERRQ(ierr);
  } else {
    ierr = TSMonitorSet(ts, PetscMonitor, this, nullptr);
    CHKERRQ(ierr);
  }

  ierr = SNESSetTolerances(snes, abstol, reltol, PETSC_DEFAULT, PETSC_DEFAULT,
                           PETSC_DEFAULT);
  CHKERRQ(ierr);

  // Matrix free Jacobian

  if (use_jacobian and hasJacobian()) {
    // Use a user-supplied Jacobian function
    ierr = MatCreateShell(comm, local_N, local_N, neq, neq, this, &Jmf);
    CHKERRQ(ierr);
    ierr = MatShellSetOperation(Jmf, MATOP_MULT,
                                reinterpret_cast<void (*)()>(PhysicsJacobianApply));
    CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts, Jmf, Jmf, solver_ijacobian, this);
    CHKERRQ(ierr);
  } else {
    // Use finite difference approximation
    ierr = MatCreateSNESMF(snes, &Jmf);
    CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, Jmf, Jmf, MatMFFDComputeJacobian, this);
    CHKERRQ(ierr);
  }

  ierr = SNESGetKSP(snes, &ksp);
  CHKERRQ(ierr);

  ierr = KSPSetTolerances(ksp, reltol, abstol, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ(ierr);

  ierr = KSPGetPC(ksp, &pc);
  CHKERRQ(ierr);

  if (use_precon and hasPreconditioner()) {

#if PETSC_VERSION_GE(3, 5, 0)
    ierr = SNESGetNPC(snes, &psnes);
    CHKERRQ(ierr);
#else
    ierr = SNESGetPC(snes, &psnes);
    CHKERRQ(ierr);
#endif
    ierr = SNESGetKSP(psnes, &nksp);
    CHKERRQ(ierr);
    ierr = KSPGetPC(nksp, &npc);
    CHKERRQ(ierr);
    ierr = SNESSetType(psnes, SNESSHELL);
    CHKERRQ(ierr);
    ierr = SNESShellSetSolve(psnes, PhysicsSNESApply);
    CHKERRQ(ierr);

    // Use a user-supplied preconditioner

    // Tell PETSc we're using a "shell" preconditioner
    ierr = PCSetType(pc, PCSHELL);
    CHKERRQ(ierr);

    // Set routine for applying preconditioner
    ierr = PCShellSetApply(pc, PhysicsPCApply);
    CHKERRQ(ierr);

    // Set context to this solver object
    ierr = PCShellSetContext(pc, this);
    CHKERRQ(ierr);

    // Set name of preconditioner
    ierr = PCShellSetName(pc, "PhysicsPreconditioner");
    CHKERRQ(ierr);

    // Need a callback for IJacobian to get shift 'alpha'
    ierr = TSSetIJacobian(ts, Jmf, Jmf, solver_ijacobian, this);
    CHKERRQ(ierr);

    // Use right preconditioner
    ierr = KSPSetPCSide(ksp, PC_RIGHT);
    CHKERRQ(ierr);

  } else {
    // Default to no preconditioner
    ierr = PCSetType(pc, PCNONE);
    CHKERRQ(ierr);
  }
  ierr = TSSetFromOptions(ts);
  CHKERRQ(ierr); // enable PETSc runtime options

  ierr = PCGetType(pc, &pctype);
  CHKERRQ(ierr);
  ierr = TSGetType(ts, &tstype);
  CHKERRQ(ierr);
  output.write("\tTS type {:s}, PC type {:s}\n", tstype, pctype);

  ierr = PetscObjectTypeCompare(reinterpret_cast<PetscObject>(pc), PCNONE, &pcnone);
  CHKERRQ(ierr);
  if (pcnone) {
    ierr = PetscLogEventEnd(init_event, 0, 0, 0, 0);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = PetscObjectTypeCompare(reinterpret_cast<PetscObject>(pc), PCSHELL, &pcnone);
  CHKERRQ(ierr);
  if (pcnone) {
    ierr = PetscLogEventEnd(init_event, 0, 0, 0, 0);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  // Create Jacobian matrix to be used by preconditioner
  output_info << " Get Jacobian matrix at simtime " << simtime << "\n";
#if PETSC_VERSION_GE(3, 7, 0)
  ierr = PetscOptionsGetString(nullptr, nullptr, "-J_load", load_file,
                               PETSC_MAX_PATH_LEN - 1, &J_load);
  CHKERRQ(ierr);
#else
  ierr = PetscOptionsGetString(nullptr, "-J_load", load_file, PETSC_MAX_PATH_LEN - 1,
                               &J_load);
  CHKERRQ(ierr);
#endif
  if (J_load) {
    PetscViewer fd;
    output_info << " Load Jmat ...local_N " << local_N << " neq " << neq << "\n";
    ierr = PetscViewerBinaryOpen(comm, load_file, FILE_MODE_READ, &fd);
    CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &J);
    CHKERRQ(ierr);
    //ierr = MatSetType(J, MATBAIJ);CHKERRQ(ierr);
    ierr = MatSetSizes(J, local_N, local_N, PETSC_DECIDE, PETSC_DECIDE);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);
    CHKERRQ(ierr);
    ierr = MatLoad(J, fd);
    CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);
    CHKERRQ(ierr);
  } else { // create Jacobian matrix

    /* number of degrees (variables) at each grid point */
    PetscInt dof = n3Dvars();

    // Maximum allowable size of stencil in x is the number of guard cells
    PetscInt stencil_width_estimate = options->operator[]("stencil_width_estimate")
                                          .withDefault(bout::globals::mesh->xstart);
    // This is the stencil in each direction (*2) along each dimension
    // (*3), plus the point itself. Not sure if this is correct
    // though, on several levels:
    //   1. Ignores corner points used in e.g. brackets
    //   2. Could have different stencil widths in each dimension
    //   3. FFTs couple every single point together
    PetscInt cols = stencil_width_estimate * 2 * 3 + 1;
    PetscInt prealloc; // = cols*dof;

    ierr = MatCreate(comm, &J);
    CHKERRQ(ierr);
    ierr = MatSetType(J, MATBAIJ);
    CHKERRQ(ierr);
    ierr = MatSetSizes(J, local_N, local_N, neq, neq);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);
    CHKERRQ(ierr);

    // Get nonzero pattern of J - color_none !!!
    prealloc = cols * dof * dof;
    ierr = MatSeqAIJSetPreallocation(J, prealloc, nullptr);
    CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(J, prealloc, nullptr, prealloc, nullptr);
    CHKERRQ(ierr);

    // why nonzeros=295900, allocated nonzeros=2816000/12800000 (*dof*dof), number of mallocs used during MatSetValues calls =256?
    prealloc = cols;
    ierr = MatSeqBAIJSetPreallocation(J, dof, prealloc, nullptr);
    CHKERRQ(ierr);
    ierr = MatMPIBAIJSetPreallocation(J, dof, prealloc, nullptr, prealloc, nullptr);
    CHKERRQ(ierr);
#if PETSC_VERSION_GE(3, 7, 0)
    ierr = PetscOptionsHasName(nullptr, nullptr, "-J_slowfd", &J_slowfd);
    CHKERRQ(ierr);
#else
    ierr = PetscOptionsHasName(nullptr, "-J_slowfd", &J_slowfd);
    CHKERRQ(ierr);
#endif
    if (J_slowfd != 0U) { // create Jacobian matrix by slow fd
      ierr = SNESSetJacobian(snes, J, J, SNESComputeJacobianDefault, nullptr);
      CHKERRQ(ierr);
      output_info << "SNESComputeJacobian J by slow fd...\n";

#if PETSC_VERSION_GE(3, 5, 0)
      ierr = TSComputeRHSJacobian(ts, simtime, u, J, J);
      CHKERRQ(ierr);
#else
      MatStructure flg;
      ierr = TSComputeRHSJacobian(ts, simtime, u, &J, &J, &flg);
      CHKERRQ(ierr);
#endif
      output_info << "compute J by slow fd is done.\n";
    } else { // get sparse pattern of the Jacobian
      throw BoutException("Sorry, unimplemented way of setting PETSc preconditioner. "
                          "Either set a preconditioner function with 'setPrecon' "
                          "in your code, or load a matrix with '-Jload' on the "
                          "command line, or calculate with finite differences "
                          "with '-J_slowfd' (on command line)");
    }
  }

  // Create coloring context of J to be used during time stepping

  output_info << " Create coloring ...\n";

  ISColoring iscoloring;
#if PETSC_VERSION_GE(3, 5, 0)
  MatColoring coloring;
  MatColoringCreate(Jmf, &coloring);
  MatColoringSetType(coloring, MATCOLORINGSL);
  MatColoringSetFromOptions(coloring);
  // Calculate index sets
  MatColoringApply(coloring, &iscoloring);
  MatColoringDestroy(&coloring);
#else
  ierr = MatGetColoring(J, MATCOLORINGSL, &iscoloring);
  CHKERRQ(ierr);
#endif
  ierr = MatFDColoringCreate(J, iscoloring, &matfdcoloring);
  CHKERRQ(ierr);

  ierr = MatFDColoringSetFromOptions(matfdcoloring);
  CHKERRQ(ierr);
  ierr = ISColoringDestroy(&iscoloring);
  CHKERRQ(ierr);
  throw BoutException("Coloring is not working");
  ierr = SNESSetJacobian(snes, J, J, SNESComputeJacobianDefaultColor, matfdcoloring);
  CHKERRQ(ierr);

  // Write J in binary for study - see ~petsc/src/mat/examples/tests/ex124.c
#if PETSC_VERSION_GE(3, 7, 0)
  ierr = PetscOptionsHasName(nullptr, nullptr, "-J_write", &J_write);
  CHKERRQ(ierr);
#else
  ierr = PetscOptionsHasName(nullptr, "-J_write", &J_write);
  CHKERRQ(ierr);
#endif
  if (J_write != 0U) {
    PetscViewer viewer = nullptr;
    output_info.write("\n[{:d}] Test TSComputeRHSJacobian() ...\n", rank);
#if PETSC_VERSION_GE(3, 5, 0)
    ierr = TSComputeRHSJacobian(ts, simtime, u, J, J);
    CHKERRQ(ierr);
#else
    MatStructure J_structure;
    ierr = TSComputeRHSJacobian(ts, simtime, u, &J, &J, &J_structure);
    CHKERRQ(ierr);
#endif

    output.write("[{:d}] TSComputeRHSJacobian is done\n", rank);

    if (J_slowfd != 0U) {
      output_info.write("[{:d}] writing J in binary to data/Jrhs_dense.dat...\n", rank);
      ierr = PetscViewerBinaryOpen(comm, "data/Jrhs_dense.dat", FILE_MODE_WRITE, &viewer);
      CHKERRQ(ierr);
    } else {
      output_info.write("[{:d}] writing J in binary to data/Jrhs_sparse.dat...\n", rank);
      ierr =
          PetscViewerBinaryOpen(comm, "data/Jrhs_sparse.dat", FILE_MODE_WRITE, &viewer);
      CHKERRQ(ierr);
    }
    ierr = MatView(J, viewer);
    CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);
    CHKERRQ(ierr);
  }

  ierr = PetscLogEventEnd(init_event, 0, 0, 0, 0);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**************************************************************************
 * Run - Advance time
 **************************************************************************/

PetscErrorCode PetscSolver::run() {

  // Set when the next call to monitor is desired
  next_output = simtime + getOutputTimestep();

  PetscFunctionBegin;

  if (this->output_flag) {
    prev_linear_its = 0;
    bout_snes_time = bout::globals::mpi->MPI_Wtime();
  }

  CHKERRQ(TSSolve(ts, u));

  // Gawd, everything is a hack
  if ((this->output_flag != 0U) and (BoutComm::rank() == 0)) {
    Output petsc_info(output_name);
    // Don't write to stdout
    petsc_info.disable();
    petsc_info.write("SNES Iteration, KSP Iterations, Wall Time, Norm\n");
    for (const auto& info : snes_list) {
      petsc_info.write("{:d}, {:d}, {:e}, {:e}\n", info.it, info.linear_its, info.time,
                       info.norm);
    }
  }

  PetscFunctionReturn(0);
}

/**************************************************************************
 * RHS function
 **************************************************************************/

PetscErrorCode PetscSolver::rhs(TS UNUSED(ts), BoutReal t, Vec udata, Vec dudata) {
  TRACE("Running RHS: PetscSolver::rhs({:e})", t);

  const BoutReal* udata_array;
  BoutReal* dudata_array;

  PetscFunctionBegin;

  // Load state from PETSc
  VecGetArrayRead(udata, &udata_array);
  load_vars(const_cast<BoutReal*>(udata_array));
  VecRestoreArrayRead(udata, &udata_array);

  // Call RHS function
  run_rhs(t);

  // Save derivatives to PETSc
  VecGetArray(dudata, &dudata_array);
  save_derivs(dudata_array);
  VecRestoreArray(dudata, &dudata_array);

  PetscFunctionReturn(0);
}

/**************************************************************************
 * Preconditioner function
 **************************************************************************/

PetscErrorCode PetscSolver::pre(PC UNUSED(pc), Vec x, Vec y) {
  TRACE("PetscSolver::pre()");

  BoutReal* data;

  if (diagnose) {
    output << "Preconditioning" << endl;
  }

  // Load state
  VecGetArray(state, &data);
  load_vars(data);
  VecRestoreArray(state, &data);

  // Load vector to be inverted into F_vars
  VecGetArray(x, &data);
  load_derivs(data);
  VecRestoreArray(x, &data);

  // Call the preconditioner
  runPreconditioner(ts_time, 1. / shift, 0.0);

  // Save the solution from time derivatives
  VecGetArray(y, &data);
  save_derivs(data);
  VecRestoreArray(y, &data);

  // Petsc's definition of Jacobian differs by a factor from Sundials'
  PetscErrorCode ierr = VecScale(y, shift);
  CHKERRQ(ierr);

  return 0;
}

/**************************************************************************
 * User-supplied Jacobian function J(state) * x = y
 **************************************************************************/

PetscErrorCode PetscSolver::jac(Vec x, Vec y) {
  TRACE("PetscSolver::jac()");

  BoutReal* data;

  if (diagnose) {
    output << "Jacobian evaluation\n";
  }

  // Load state
  VecGetArray(state, &data);
  load_vars(data);
  VecRestoreArray(state, &data);

  // Load vector to be operated on into F_vars
  VecGetArray(x, &data);
  load_derivs(data);
  VecRestoreArray(x, &data);

  // Call the Jacobian function
  runJacobian(ts_time);

  // Save the solution from time derivatives
  VecGetArray(y, &data);
  save_derivs(data);
  VecRestoreArray(y, &data);

  // y = a * x - y
  PetscErrorCode ierr = VecAXPBY(y, shift, -1.0, x);
  CHKERRQ(ierr);

  return 0;
}

/**************************************************************************
 * Static functions which can be used for PETSc callbacks
 **************************************************************************/

PetscErrorCode solver_f(TS ts, BoutReal t, Vec globalin, Vec globalout, void* f_data) {
  PetscFunctionBegin;
  auto* s = static_cast<PetscSolver*>(f_data);
  PetscLogEventBegin(s->solver_event, 0, 0, 0, 0);
  s->rhs(ts, t, globalin, globalout);
  PetscLogEventEnd(s->solver_event, 0, 0, 0, 0);
  PetscFunctionReturn(0);
}

/*
  FormIFunction = Udot - RHSFunction
*/
PetscErrorCode solver_if(TS ts, BoutReal t, Vec globalin, Vec globalindot, Vec globalout,
                         void* f_data) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = solver_f(ts, t, globalin, globalout, f_data);
  CHKERRQ(ierr);

  ierr = VecAYPX(globalout, -1.0, globalindot);
  CHKERRQ(ierr); // globalout = globalindot + (-1)globalout
  PetscFunctionReturn(0);
}

#if PETSC_VERSION_GE(3, 5, 0)
PetscErrorCode solver_rhsjacobian(TS UNUSED(ts), BoutReal UNUSED(t), Vec UNUSED(globalin),
                                  Mat J, Mat Jpre, void* UNUSED(f_data)) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  if (J != Jpre) {
    ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
#else
PetscErrorCode solver_rhsjacobian([[maybe_unused]] TS ts, [[maybe_unused]] BoutReal t,
                                  [[maybe_unused]] Vec globalin, Mat* J, Mat* Jpre,
                                  [[maybe_unused]] MatStructure* str,
                                  [[maybe_unused]] void* f_data) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatAssemblyBegin(*Jpre, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*Jpre, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  if (*J != *Jpre) {
    ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
#endif

/*
  solver_ijacobian - Compute IJacobian = dF/dU + a dF/dUdot  - a dummy matrix used for pc=none
*/
#if PETSC_VERSION_GE(3, 5, 0)
PetscErrorCode solver_ijacobian(TS ts, BoutReal t, Vec globalin, Vec UNUSED(globalindot),
                                PetscReal a, Mat J, Mat Jpre, void* f_data) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = solver_rhsjacobian(ts, t, globalin, J, Jpre, f_data);
  CHKERRQ(ierr);

  ////// Save data for preconditioner
  auto* solver = static_cast<PetscSolver*>(f_data);

  if (solver->diagnose) {
    output << "Saving state, t = " << t << ", a = " << a << endl;
  }

  solver->shift = a;        // Save the shift 'a'
  solver->state = globalin; // Save system state
  solver->ts_time = t;

  PetscFunctionReturn(0);
}
#else
PetscErrorCode solver_ijacobian(TS ts, BoutReal t, Vec globalin,
                                [[maybe_unused]] Vec globalindot,
                                [[maybe_unused]] PetscReal a, Mat* J, Mat* Jpre,
                                MatStructure* str, void* f_data) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = solver_rhsjacobian(ts, t, globalin, J, Jpre, str, (void*)f_data);
  CHKERRQ(ierr);

  ////// Save data for preconditioner
  PetscSolver* solver = (PetscSolver*)f_data;

  if (solver->diagnose) {
    output << "Saving state, t = " << t << ", a = " << a << endl;
  }

  solver->shift = a;        // Save the shift 'a'
  solver->state = globalin; // Save system state
  solver->ts_time = t;

  PetscFunctionReturn(0);
}
#endif

/*
  solver_ijacobianfd - Compute IJacobian = dF/dU + a dF/dUdot  using finite deference - not implemented yet
*/
#if PETSC_VERSION_GE(3, 5, 0)
PetscErrorCode solver_ijacobianfd(TS ts, BoutReal t, Vec globalin,
                                  Vec UNUSED(globalindot), PetscReal UNUSED(a), Mat J,
                                  Mat Jpre, void* f_data) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = solver_rhsjacobian(ts, t, globalin, J, Jpre, f_data);
  CHKERRQ(ierr);
  //*Jpre + a
  PetscFunctionReturn(0);
}
#else
PetscErrorCode solver_ijacobianfd(TS ts, BoutReal t, Vec globalin, Vec globalindot,
                                  PetscReal a, Mat* J, Mat* Jpre, MatStructure* str,
                                  void* f_data) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = solver_rhsjacobian(ts, t, globalin, J, Jpre, str, (void*)f_data);
  CHKERRQ(ierr);
  //*Jpre + a
  PetscFunctionReturn(0);
}
#endif
//-----------------------------------------

PetscErrorCode PhysicsSNESApply(SNES snes, Vec x) {
  PetscErrorCode ierr;
  Vec F, Fout;
  PetscReal fnorm = 0., foutnorm = 0., dot = 0.;
  KSP ksp;
  PC pc;
  Mat A, B;

  PetscFunctionBegin;
  ierr = SNESGetJacobian(snes, &A, &B, nullptr, nullptr);
  CHKERRQ(ierr);
#if PETSC_VERSION_GE(3, 5, 0)
  ierr = SNESComputeJacobian(snes, x, A, B);
  CHKERRQ(ierr);
#else
  MatStructure diff = DIFFERENT_NONZERO_PATTERN;
  ierr = SNESComputeJacobian(snes, x, &A, &B, &diff);
  CHKERRQ(ierr);
#endif
  ierr = SNESGetKSP(snes, &ksp);
  CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc);
  CHKERRQ(ierr);
  ierr = SNESGetFunction(snes, &F, nullptr, nullptr);
  CHKERRQ(ierr);
  ierr = SNESComputeFunction(snes, x, F);
  CHKERRQ(ierr);
  ierr = SNESGetSolutionUpdate(snes, &Fout);
  CHKERRQ(ierr);

  ierr = PCApply(pc, F, Fout);
  CHKERRQ(ierr);
  ierr = VecNorm(Fout, NORM_2, &foutnorm);
  CHKERRQ(ierr);
  ierr = VecAXPY(x, -1., Fout);
  CHKERRQ(ierr);
  ierr = SNESComputeFunction(snes, x, F);
  CHKERRQ(ierr);
  ierr = VecNorm(F, NORM_2, &fnorm);
  CHKERRQ(ierr);
  ierr = VecDot(F, Fout, &dot);
  CHKERRQ(ierr);
  output_info << " (Debug) function norm: " << fnorm << ", P(f) norm " << foutnorm
              << ", F \\cdot Fout " << dot << "  ";
#if PETSC_VERSION_GE(3, 5, 0)
  Vec func;
  ierr = SNESGetFunction(snes, &func, nullptr, nullptr);
  CHKERRQ(ierr);
  ierr = VecNorm(func, NORM_2, &fnorm);
  CHKERRQ(ierr);
#else
  ierr = SNESSetFunctionNorm(snes, fnorm);
  CHKERRQ(ierr);
#endif
  ierr = SNESMonitor(snes, 0, fnorm);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PhysicsPCApply(PC pc, Vec x, Vec y) {
  int ierr;

  // Get the context
  PetscSolver* s;
  ierr = PCShellGetContext(pc, reinterpret_cast<void**>(&s));
  CHKERRQ(ierr);

  PetscFunctionReturn(s->pre(pc, x, y));
}

PetscErrorCode PhysicsJacobianApply(Mat J, Vec x, Vec y) {
  // Get the context
  PetscSolver* s;
  int ierr = MatShellGetContext(J, reinterpret_cast<void**>(&s));
  CHKERRQ(ierr);
  PetscFunctionReturn(s->jac(x, y));
}

PetscErrorCode PetscMonitor(TS ts, PetscInt UNUSED(step), PetscReal t, Vec X, void* ctx) {
  PetscErrorCode ierr;
  auto* s = static_cast<PetscSolver*>(ctx);
  PetscReal tfinal, dt;
  Vec interpolatedX;
  const PetscScalar* x;
  static int i = 0;

  PetscFunctionBegin;
  ierr = TSGetTimeStep(ts, &dt);
  CHKERRQ(ierr);

#if PETSC_VERSION_GE(3, 8, 0)
  ierr = TSGetMaxTime(ts, &tfinal);
  CHKERRQ(ierr);
#else
  ierr = TSGetDuration(ts, nullptr, &tfinal);
  CHKERRQ(ierr);
#endif

  /* Duplicate the solution vector X into a work vector */
  ierr = VecDuplicate(X, &interpolatedX);
  CHKERRQ(ierr);
  while (s->next_output <= t && s->next_output <= tfinal) {
    if (s->interpolate) {
      ierr = TSInterpolate(ts, s->next_output, interpolatedX);
      CHKERRQ(ierr);
    }

    /* Place the interpolated values into the global variables */
    ierr = VecGetArrayRead(interpolatedX, &x);
    CHKERRQ(ierr);
    s->load_vars(const_cast<BoutReal*>(x));
    ierr = VecRestoreArrayRead(interpolatedX, &x);
    CHKERRQ(ierr);

    if (s->call_monitors(simtime, i++, s->getNumberOutputSteps())) {
      PetscFunctionReturn(1);
    }

    s->next_output += s->getOutputTimestep();
    simtime = s->next_output;
  }

  /* Done with vector, so destroy it */
  ierr = VecDestroy(&interpolatedX);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PetscSNESMonitor(SNES snes, PetscInt its, PetscReal norm, void* ctx) {
  PetscErrorCode ierr;
  PetscInt linear_its = 0;
  BoutReal tmp = .0;
  snes_info row;
  auto* s = static_cast<PetscSolver*>(ctx);

  PetscFunctionBegin;

  if (!its) {
    s->prev_linear_its = 0;
  }
  ierr = SNESGetLinearSolveIterations(snes, &linear_its);
  CHKERRQ(ierr);
  tmp = bout::globals::mpi->MPI_Wtime();

  row.it = its;
  s->prev_linear_its = row.linear_its = linear_its - s->prev_linear_its;
  row.time = tmp - s->bout_snes_time;
  row.norm = norm;

  s->snes_list.push_back(row);

  PetscFunctionReturn(0);
}

#endif
