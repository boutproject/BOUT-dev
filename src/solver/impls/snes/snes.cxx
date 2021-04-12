
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
static PetscErrorCode FormFunction(SNES UNUSED(snes), Vec x, Vec f, void* ctx) {
  return static_cast<SNESSolver*>(ctx)->snes_function(x, f);
}

int SNESSolver::init(int nout, BoutReal tstep) {

  TRACE("Initialising SNES solver");

  /// Call the generic initialisation first
  if (Solver::init(nout, tstep) != 0) {
    return 1;
  }

  out_timestep = tstep; // Output timestep
  nsteps = nout;        // Save number of output steps

  output << "\n\tSNES steady state solver\n";

  // Calculate number of variables
  nlocal = getLocalN();

  // Get total problem size
  int ntmp;
  if (MPI_Allreduce(&nlocal, &ntmp, 1, MPI_INT, MPI_SUM, BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed!");
  }
  neq = ntmp;
  
  output.write("\t3d fields = %d, 2d fields = %d neq=%d, local_N=%d\n",
	       n3Dvars(), n2Dvars(), neq, nlocal);

  timestep =
      (*options)["timestep"].doc("Initial backward Euler timestep").withDefault(1.0);

  diagnose =
      (*options)["diagnose"].doc("Print additional diagnostics").withDefault(false);

  predictor = (*options)["predictor"].doc("Use linear predictor?").withDefault(true);

  // Initialise PETSc components
  int ierr;

  // Vectors
  ierr = VecCreate(BoutComm::get(), &snes_x);
  CHKERRQ(ierr); // NOLINT
  ierr = VecSetSizes(snes_x, nlocal, PETSC_DECIDE);
  CHKERRQ(ierr); // NOLINT
  ierr = VecSetFromOptions(snes_x);
  CHKERRQ(ierr); // NOLINT

  VecDuplicate(snes_x, &snes_f);
  VecDuplicate(snes_x, &x0);

  if (predictor) {
    // Storage for previous solution
    VecDuplicate(snes_x, &x1);
  }

  // Nonlinear solver interface (SNES)
  SNESCreate(BoutComm::get(), &snes);

  // Set the callback function
  SNESSetFunction(snes, snes_f, FormFunction, this);

  std::string snes_type = (*options)["snes_type"].withDefault("anderson");
  SNESSetType(snes, snes_type.c_str());

  // Set up the Jacobian
  MatCreateAIJ(BoutComm::get(), nlocal, nlocal,  // Local sizes
               PETSC_DETERMINE, PETSC_DETERMINE, // Global sizes
               3, // Number of nonzero entries in diagonal portion of local submatrix
               PETSC_NULL,
               0, // Number of nonzeros per row in off-diagonal portion of local submatrix
               PETSC_NULL, &Jmf);
#if PETSC_VERSION_GE(3, 4, 0)
  SNESSetJacobian(snes, Jmf, Jmf, SNESComputeJacobianDefault, this);
#else
  // Before 3.4
  SNESSetJacobian(snes, Jmf, Jmf, SNESDefaultComputeJacobian, this);
#endif
  MatSetOption(Jmf, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  // Set tolerances
  BoutReal atol =
      (*options)["atol"].doc("Absolute tolerance in SNES solve").withDefault(1e-16);
  BoutReal rtol =
      (*options)["rtol"].doc("Relative tolerance in SNES solve").withDefault(1e-10);

  int maxits = (*options)["max_nonlinear_it"]
                   .doc("Maximum number of iterations per SNES solve")
                   .withDefault(50);

  upper_its = (*options)["upper_its"]
                  .doc("Iterations above which the next timestep is reduced")
                  .withDefault(static_cast<int>(maxits * 0.8));

  lower_its = (*options)["lower_its"]
                  .doc("Iterations below which the next timestep is increased")
                  .withDefault(static_cast<int>(maxits * 0.5));

  SNESSetTolerances(snes, atol, rtol, PETSC_DEFAULT, maxits, PETSC_DEFAULT);

  // Get runtime options
  SNESSetFromOptions(snes);

  return 0;
}

int SNESSolver::run() {
  TRACE("SNESSolver::run()");
  // Set initial guess at the solution from variables
  {
    BoutReal* xdata = nullptr;
    int ierr = VecGetArray(snes_x, &xdata);
    CHKERRQ(ierr); // NOLINT
    save_vars(xdata);
    ierr = VecRestoreArray(snes_x, &xdata);
    CHKERRQ(ierr); // NOLINT
  }

  for (int s = 0; s < nsteps; s++) {
    BoutReal target = simtime + out_timestep;

    bool looping = true;
    do {
      // Copy the state (snes_x) into initial values (x0)
      VecCopy(snes_x, x0);

      // Set the timestep
      dt = timestep;
      looping = true;
      if (simtime + dt >= target) {
        // Ensure that the timestep goes to the next output time and then stops
        looping = false;
        dt = target - simtime;
      }

      if (predictor and (time1 > 0.0)) {
        // Use (time1, x1) and (simtime, x0) to make prediction
        // snes_x <- x0 + (dt / (simtime - time1)) * (x0 - x1)
        // snes_x <- -β * x1 + (1 + β) * snes_x
        BoutReal beta = dt / (simtime - time1);
        VecAXPBY(snes_x, -beta, (1. + beta), x1);
      }

      // Run the solver
      SNESSolve(snes, nullptr, snes_x);

      // Find out if converged
      SNESConvergedReason reason;
      SNESGetConvergedReason(snes, &reason);
      if (reason < 0) {
        // Diverged

        // Try a smaller timestep
        timestep /= 2.0;
        // Restore state
        VecCopy(x0, snes_x);

        // Check lock state
        PetscInt lock_state;
        VecLockGet(snes_x, &lock_state);
        if (lock_state > 0) {
          // Locked for read
          output.write("WARNING: snes_x locked for reading\n");
        } else if (lock_state < 0) {
          // Locked for write
          output.write("WARNING: snes_x locked for writing\n");
        }
        continue; // Try again
      }

      if (predictor) {
        // Save previous values: x1 <- x0
        VecCopy(x0, x1);
        time1 = simtime;
      }

      simtime += dt;
      int its;
      SNESGetIterationNumber(snes, &its);

      if (diagnose) {
        output.write("SNES time: %e, timestep: %e, iterations: %d\n", simtime, timestep,
                     its);
      }

      if (looping) {
        if (its <= lower_its) {
          // Increase timestep slightly
          timestep *= 1.1;
        } else if (its >= upper_its) {
          // Reduce timestep slightly
          timestep *= 0.9;
        }
      }
    } while (looping);

    // Put the result into variables
    {
      const BoutReal* xdata = nullptr;
      int ierr = VecGetArrayRead(snes_x, &xdata);
      CHKERRQ(ierr); // NOLINT
      load_vars(const_cast<BoutReal*>(xdata));
      ierr = VecRestoreArrayRead(snes_x, &xdata);
      CHKERRQ(ierr); // NOLINT
    }
    run_rhs(simtime); // Run RHS to calculate auxilliary variables

    /// Call the monitor function

    if (call_monitors(simtime, s, nsteps) != 0) {
      break; // User signalled to quit
    }
  }

  return 0;
}

// f = rhs
PetscErrorCode SNESSolver::snes_function(Vec x, Vec f) {
  // Get data from PETSc into BOUT++ fields
  const BoutReal* xdata = nullptr;
  int ierr = VecGetArrayRead(x, &xdata);
  CHKERRQ(ierr); // NOLINT
  load_vars(const_cast<BoutReal*>(xdata));
  ierr = VecRestoreArrayRead(x, &xdata);
  CHKERRQ(ierr); // NOLINT

  try {
    // Call RHS function
    run_rhs(simtime + dt);
  } catch (BoutException& e) {
    // Simulation might fail, e.g. negative densities
    // if timestep too large
    output.write("WARNING: BoutException thrown: %s\n", e.what());

    // Tell SNES that the input was out of domain
    SNESSetFunctionDomainError(snes);
    // Note: Returning non-zero error here leaves vectors in locked state
    return 0;
  }

  // Copy derivatives back
  BoutReal* fdata = nullptr;
  ierr = VecGetArray(f, &fdata);
  CHKERRQ(ierr); // NOLINT
  save_derivs(fdata);
  ierr = VecRestoreArray(f, &fdata);
  CHKERRQ(ierr); // NOLINT

  // Backward Euler
  // Set fdata = xdata - x0 - Δt*fdata
  VecAYPX(f, -dt, x);   // f <- x - Δt*f
  VecAXPY(f, -1.0, x0); // f <- f - x0

  return 0;
}

#endif // BOUT_HAS_PETSC
