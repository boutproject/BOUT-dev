
#ifdef BOUT_HAS_PETSC

#include "snes.hxx"

#include <boutcomm.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include <output.hxx>

#include "petscsnes.h"

/*
 * PETSc callback function, which evaluates the nonlinear
 * function to be solved by SNES.
 *
 * This function assumes the context void pointer is a pointer
 * to an SNESSolver object.
 */
static PetscErrorCode FormFunction(SNES UNUSED(snes), Vec x, Vec f, void *ctx) {
  return static_cast<SNESSolver*>(ctx)->snes_function(x, f);
}

int SNESSolver::init(int nout, BoutReal tstep) {

  TRACE("Initialising SNES solver");
  
  /// Call the generic initialisation first
  if (Solver::init(nout, tstep))
    return 1;
  
  output << "\n\tSNES steady state solver\n";
  
  // Calculate number of variables
  nlocal = getLocalN();
  
  // Get total problem size
  int ntmp;
  if(MPI_Allreduce(&nlocal, &ntmp, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed!");
  }
  neq = ntmp;
  
  output.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n",
	       n3Dvars(), n2Dvars(), neq, nlocal);
  
  // Get options
  OPTION(options, mxstep, 500); // Maximum number of steps between outputs
  
  // Initialise PETSc components
  int ierr;
  
  // Vectors
  ierr = VecCreate(BoutComm::get(), &snes_x);CHKERRQ(ierr);
  ierr = VecSetSizes(snes_x, nlocal, PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(snes_x);CHKERRQ(ierr);
  
  VecDuplicate(snes_x,&snes_f);
  
  // Set initial guess at the solution from variables
  BoutReal *xdata;
  ierr = VecGetArray(snes_x,&xdata);CHKERRQ(ierr);
  save_vars(xdata);
  ierr = VecRestoreArray(snes_x,&xdata);CHKERRQ(ierr);
  
  // Nonlinear solver interface (SNES)
  SNESCreate(BoutComm::get(),&snes);
  
  // Set the callback function
  SNESSetFunction(snes,snes_f,FormFunction,this);
  
  // Set up the Jacobian
  //MatCreateSNESMF(snes,&Jmf);
  //SNESSetJacobian(snes,Jmf,Jmf,SNESComputeJacobianDefault,this);
  MatCreateAIJ(BoutComm::get(),
               nlocal,nlocal,  // Local sizes
               PETSC_DETERMINE, PETSC_DETERMINE, // Global sizes
               3,   // Number of nonzero entries in diagonal portion of local submatrix
               PETSC_NULL,
               0,   // Number of nonzeros per row in off-diagonal portion of local submatrix
               PETSC_NULL, 
               &Jmf);
#if PETSC_VERSION_GE(3,4,0)
  SNESSetJacobian(snes,Jmf,Jmf,SNESComputeJacobianDefault,this);
#else
  // Before 3.4
  SNESSetJacobian(snes,Jmf,Jmf,SNESDefaultComputeJacobian,this);
#endif
  MatSetOption(Jmf,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);

  // Set tolerances
  BoutReal atol, rtol; // Tolerances for SNES solver
  options->get("atol", atol, 1e-16);
  options->get("rtol", rtol, 1e-10);
  SNESSetTolerances(snes,atol,rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  
  // Get runtime options
  SNESSetFromOptions(snes);
  
  return 0;
}

int SNESSolver::run() {
  TRACE("SNESSolver::run()");
  
  /*
  output << "Computing Jacobian\n";
  MatStructure  flag;
  implicit_curtime = curtime;
  implicit_gamma = gamma;
  SNESComputeFunction(snes, snes_x, snes_f);
  SNESComputeJacobian(snes,snes_x,&Jmf,&Jmf,&flag);
  MatView(Jmf, 	PETSC_VIEWER_STDOUT_SELF);
  */
  SNESSolve(snes, nullptr, snes_x);

  // Find out if converged
  SNESConvergedReason reason;
  SNESGetConvergedReason(snes,&reason);
  if(reason < 0) {
    // Diverged
    throw BoutException("SNES failed to converge. Reason: {:d}\n", reason);
  }
  
  int its;
  SNESGetIterationNumber(snes,&its);
  
  //output << "Number of SNES iterations: " << its << endl;
  
  // Put the result into variables
  const BoutReal *xdata;
  int ierr;
  ierr = VecGetArrayRead(snes_x,&xdata);CHKERRQ(ierr);
  load_vars(const_cast<BoutReal*>(xdata));
  ierr = VecRestoreArrayRead(snes_x,&xdata);CHKERRQ(ierr);

  run_rhs(0.0); // Run RHS to calculate auxilliary variables
    
  /// Call the monitor function
  
  if(call_monitors(0.0, 1, 1)) {
    // User signalled to quit
  }
    
  return 0;
}

// f = rhs
PetscErrorCode SNESSolver::snes_function(Vec x, Vec f) {
  const BoutReal *xdata;
  BoutReal *fdata;
  int ierr;

  // Get data from PETSc into BOUT++ fields
  ierr = VecGetArrayRead(x,&xdata);CHKERRQ(ierr);
  load_vars(const_cast<BoutReal*>(xdata));
  ierr = VecRestoreArrayRead(x,&xdata);CHKERRQ(ierr);

  // Call RHS function
  run_rhs(0.0);
  
  // Copy derivatives back
  ierr = VecGetArray(f,&fdata);CHKERRQ(ierr);
  save_derivs(fdata);
  ierr = VecRestoreArray(f,&fdata);CHKERRQ(ierr);
  
  return 0;
}

#endif // BOUT_HAS_PETSC
