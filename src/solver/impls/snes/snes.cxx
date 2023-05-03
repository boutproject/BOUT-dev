#include "bout/build_config.hxx"

#if BOUT_HAS_PETSC

#include "snes.hxx"

#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/msg_stack.hxx>
#include <bout/utils.hxx>

#include <array>
#include <cmath>
#include <vector>

#include <bout/output.hxx>

#include "petscsnes.h"
/*
 * PETSc callback function, which evaluates the nonlinear
 * function to be solved by SNES.
 *
 * This function assumes the context void pointer is a pointer
 * to an SNESSolver object.
 */
static PetscErrorCode FormFunction(SNES UNUSED(snes), Vec x, Vec f, void* ctx) {
  return static_cast<SNESSolver*>(ctx)->snes_function(x, f, false);
}

/*!
 * PETSc callback function for forming Jacobian
 *
 * This function can be a linearised form of FormFunction
 */
static PetscErrorCode FormFunctionForDifferencing(void* ctx, Vec x, Vec f) {
  return static_cast<SNESSolver*>(ctx)->snes_function(x, f, true);
}

/*!
 * SNES callback for forming Jacobian with coloring
 *
 * This can be a linearised and simplified form of FormFunction
 */
static PetscErrorCode FormFunctionForColoring(SNES UNUSED(snes), Vec x, Vec f,
                                              void* ctx) {
  return static_cast<SNESSolver*>(ctx)->snes_function(x, f, true);
}

static PetscErrorCode snesPCapply(PC pc, Vec x, Vec y) {
  int ierr;

  // Get the context
  SNESSolver* s;
  ierr = PCShellGetContext(pc, reinterpret_cast<void**>(&s));
  CHKERRQ(ierr);

  PetscFunctionReturn(s->precon(x, y));
}

SNESSolver::SNESSolver(Options* opts)
    : Solver(opts),
      timestep(
          (*options)["timestep"].doc("Initial backward Euler timestep").withDefault(1.0)),
      dt_min_reset((*options)["dt_min_reset"]
                       .doc("If dt falls below this, reset to starting dt")
                       .withDefault(1e-6)),
      max_timestep((*options)["max_timestep"].doc("Maximum timestep").withDefault(1e37)),
      snes_type((*options)["snes_type"]
                    .doc("PETSc nonlinear solver method to use")
                    .withDefault("anderson")),
      atol((*options)["atol"].doc("Absolute tolerance in SNES solve").withDefault(1e-16)),
      rtol((*options)["rtol"].doc("Relative tolerance in SNES solve").withDefault(1e-10)),
      stol((*options)["stol"]
               .doc("Convergence tolerance in terms of the norm of the change in "
                    "the solution between steps")
               .withDefault(1e-8)),
      maxits((*options)["max_nonlinear_iterations"]
                 .doc("Maximum number of nonlinear iterations per SNES solve")
                 .withDefault(50)),
      lower_its((*options)["lower_its"]
                    .doc("Iterations below which the next timestep is increased")
                    .withDefault(static_cast<int>(maxits * 0.5))),
      upper_its((*options)["upper_its"]
                    .doc("Iterations above which the next timestep is reduced")
                    .withDefault(static_cast<int>(maxits * 0.8))),
      diagnose(
          (*options)["diagnose"].doc("Print additional diagnostics").withDefault(false)),
      diagnose_failures((*options)["diagnose_failures"]
                            .doc("Print more diagnostics when SNES fails")
                            .withDefault<bool>(false)),
      equation_form(
          (*options)["equation_form"]
              .doc("Form of equation to solve: rearranged_backward_euler (default);"
                   " pseudo_transient; backward_euler; direct_newton")
              .withDefault(BoutSnesEquationForm::rearranged_backward_euler)),
      predictor((*options)["predictor"].doc("Use linear predictor?").withDefault(true)),
      use_precon((*options)["use_precon"]
                     .doc("Use user-supplied preconditioner?")
                     .withDefault<bool>(false)),
      ksp_type((*options)["ksp_type"]
                   .doc("Linear solver type. By default let PETSc decide (gmres)")
                   .withDefault("default")),
      kspsetinitialguessnonzero((*options)["kspsetinitialguessnonzero"]
                                    .doc("Set the initial guess to be non-zero")
                                    .withDefault<bool>(false)),
      maxl((*options)["maxl"].doc("Maximum number of linear iterations").withDefault(20)),
      pc_type(
          (*options)["pc_type"]
              .doc("Preconditioner type. By default lets PETSc decide (ilu or bjacobi)")
              .withDefault("default")),
      line_search_type((*options)["line_search_type"]
                           .doc("Line search type: basic, bt, l2, cp, nleqerr")
                           .withDefault("default")),
      matrix_free((*options)["matrix_free"]
                      .doc("Use matrix free Jacobian?")
                      .withDefault<bool>(false)),
      lag_jacobian((*options)["lag_jacobian"]
                       .doc("Re-use the Jacobian this number of SNES iterations")
                       .withDefault(50)),
      use_coloring((*options)["use_coloring"]
                       .doc("Use matrix coloring to calculate Jacobian?")
                       .withDefault<bool>(true)) {}

int SNESSolver::init() {

  TRACE("Initialising SNES solver");

  Solver::init();
  output << "\n\tSNES steady state solver\n";

  // Calculate number of variables
  nlocal = getLocalN();

  // Get total problem size
  int ntmp;
  if (bout::globals::mpi->MPI_Allreduce(&nlocal, &ntmp, 1, MPI_INT, MPI_SUM,
                                        BoutComm::get())) {
    throw BoutException("MPI_Allreduce failed!");
  }
  neq = ntmp;

  output_info.write("\t3d fields = {:d}, 2d fields = {:d} neq={:d}, local_N={:d}\n",
                    n3Dvars(), n2Dvars(), neq, nlocal);

  // Initialise PETSc components
  int ierr;

  // Vectors
  ierr = VecCreate(BoutComm::get(), &snes_x);
  CHKERRQ(ierr);
  ierr = VecSetSizes(snes_x, nlocal, PETSC_DECIDE);
  CHKERRQ(ierr);
  ierr = VecSetFromOptions(snes_x);
  CHKERRQ(ierr);

  VecDuplicate(snes_x, &snes_f);
  VecDuplicate(snes_x, &x0);

  if (equation_form == BoutSnesEquationForm::rearranged_backward_euler) {
    // Need an intermediate vector for rearranged Backward Euler
    VecDuplicate(snes_x, &delta_x);
  }

  if (predictor) {
    // Storage for previous solution
    VecDuplicate(snes_x, &x1);
  }

  // Nonlinear solver interface (SNES)
  SNESCreate(BoutComm::get(), &snes);

  // Set the callback function
  SNESSetFunction(snes, snes_f, FormFunction, this);

  SNESSetType(snes, snes_type.c_str());

  // Line search
  if (line_search_type != "default") {
    SNESLineSearch linesearch;
    SNESGetLineSearch(snes, &linesearch);
    SNESLineSearchSetType(linesearch, line_search_type.c_str());
  }

  // Set up the Jacobian
  if (matrix_free) {
    /*
      PETSc SNES matrix free Jacobian, using a different
      operator for differencing.

      See PETSc examples
      http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tests/ex7.c.html
      and this thread:
      http://lists.mcs.anl.gov/pipermail/petsc-users/2014-January/020075.html

     */
    MatCreateSNESMF(snes, &Jmf);

    // Set a function to be called for differencing
    // This can be a linearised form of the SNES function
    MatMFFDSetFunction(Jmf, FormFunctionForDifferencing, this);

    // Calculate Jacobian matrix free using FormFunctionForDifferencing
    SNESSetJacobian(snes, Jmf, Jmf, MatMFFDComputeJacobian, this);

  } else {
    // Calculate the Jacobian using finite differences
    if (use_coloring) {
      // Use matrix coloring
      // This greatly reduces the number of times the rhs() function needs
      // to be evaluated when calculating the Jacobian.

      // Use global mesh for now
      Mesh* mesh = bout::globals::mesh;

      //////////////////////////////////////////////////
      // Get the local indices by starting at 0
      Field3D index = globalIndex(0);

      //////////////////////////////////////////////////
      // Pre-allocate PETSc storage

      output_progress.write("Setting Jacobian matrix sizes\n");

      int localN = getLocalN(); // Number of rows on this processor
      int n2d = f2d.size();
      int n3d = f3d.size();

      // Set size of Matrix on each processor to localN x localN
      MatCreate(BoutComm::get(), &Jmf);
      MatSetSizes(Jmf, localN, localN, PETSC_DETERMINE, PETSC_DETERMINE);
      MatSetFromOptions(Jmf);

      std::vector<PetscInt> d_nnz(localN);
      std::vector<PetscInt> o_nnz(localN);

      // Set values for most points

      const auto star_pattern_2d = 5 * (n3d + n2d);
      const auto star_pattern_3d = 7 * n3d + 5 * n2d;
      const auto star_pattern = (mesh->LocalNz > 1) ? star_pattern_3d : star_pattern_2d;
      for (int i = 0; i < localN; i++) {
        // Non-zero elements on this processor
        d_nnz[i] = star_pattern;
        // Non-zero elements on neighboring processor
        o_nnz[i] = 0;
      }

      // X boundaries
      if (mesh->firstX()) {
        // Lower X boundary
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          for (int z = 0; z < mesh->LocalNz; z++) {
            int localIndex = ROUND(index(mesh->xstart, y, z));
            ASSERT2((localIndex >= 0) && (localIndex < localN));
            const int num_fields = (z == 0) ? n2d + n3d : n3d;
            for (int i = 0; i < num_fields; i++) {
              d_nnz[localIndex + i] -= (n3d + n2d);
            }
          }
        }
      } else {
        // On another processor
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          for (int z = 0; z < mesh->LocalNz; z++) {
            int localIndex = ROUND(index(mesh->xstart, y, z));
            ASSERT2((localIndex >= 0) && (localIndex < localN));
            const int num_fields = (z == 0) ? n2d + n3d : n3d;
            for (int i = 0; i < num_fields; i++) {
              d_nnz[localIndex + i] -= (n3d + n2d);
              o_nnz[localIndex + i] += (n3d + n2d);
            }
          }
        }
      }

      if (mesh->lastX()) {
        // Upper X boundary
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          for (int z = 0; z < mesh->LocalNz; z++) {
            int localIndex = ROUND(index(mesh->xend, y, z));
            ASSERT2((localIndex >= 0) && (localIndex < localN));
            const int num_fields = (z == 0) ? n2d + n3d : n3d;
            for (int i = 0; i < num_fields; i++) {
              d_nnz[localIndex + i] -= (n3d + n2d);
            }
          }
        }
      } else {
        // On another processor
        for (int y = mesh->ystart; y <= mesh->yend; y++) {
          for (int z = 0; z < mesh->LocalNz; z++) {
            int localIndex = ROUND(index(mesh->xend, y, z));
            ASSERT2((localIndex >= 0) && (localIndex < localN));
            const int num_fields = (z == 0) ? n2d + n3d : n3d;
            for (int i = 0; i < num_fields; i++) {
              d_nnz[localIndex + i] -= (n3d + n2d);
              o_nnz[localIndex + i] += (n3d + n2d);
            }
          }
        }
      }

      // Y boundaries

      for (int x = mesh->xstart; x <= mesh->xend; x++) {
        // Default to no boundary
        // NOTE: This assumes that communications in Y are to other
        //   processors. If Y is communicated with this processor (e.g. NYPE=1)
        //   then this will result in PETSc warnings about out of range allocations

        // z = 0 case
        int localIndex = ROUND(index(x, mesh->ystart, 0));
        ASSERT2(localIndex >= 0);

        // All 2D and 3D fields
        for (int i = 0; i < n2d + n3d; i++) {
          o_nnz[localIndex + i] += (n3d + n2d);
        }

        for (int z = 1; z < mesh->LocalNz; z++) {
          localIndex = ROUND(index(x, mesh->ystart, z));

          // Only 3D fields
          for (int i = 0; i < n3d; i++) {
            o_nnz[localIndex + i] += (n3d + n2d);
          }
        }

        // z = 0 case
        localIndex = ROUND(index(x, mesh->yend, 0));
        // All 2D and 3D fields
        for (int i = 0; i < n2d + n3d; i++) {
          o_nnz[localIndex + i] += (n3d + n2d);
        }

        for (int z = 1; z < mesh->LocalNz; z++) {
          localIndex = ROUND(index(x, mesh->yend, z));

          // Only 3D fields
          for (int i = 0; i < n3d; i++) {
            o_nnz[localIndex + i] += (n3d + n2d);
          }
        }
      }

      for (RangeIterator it = mesh->iterateBndryLowerY(); !it.isDone(); it++) {
        // A boundary, so no communication

        // z = 0 case
        int localIndex = ROUND(index(it.ind, mesh->ystart, 0));
        if (localIndex < 0) {
          // This can occur because it.ind includes values in x boundary e.g. x=0
          continue;
        }
        // All 2D and 3D fields
        for (int i = 0; i < n2d + n3d; i++) {
          o_nnz[localIndex + i] -= (n3d + n2d);
        }

        for (int z = 1; z < mesh->LocalNz; z++) {
          int localIndex = ROUND(index(it.ind, mesh->ystart, z));

          // Only 3D fields
          for (int i = 0; i < n3d; i++) {
            o_nnz[localIndex + i] -= (n3d + n2d);
          }
        }
      }

      for (RangeIterator it = mesh->iterateBndryUpperY(); !it.isDone(); it++) {
        // A boundary, so no communication

        // z = 0 case
        int localIndex = ROUND(index(it.ind, mesh->yend, 0));
        if (localIndex < 0) {
          continue; // Out of domain
        }

        // All 2D and 3D fields
        for (int i = 0; i < n2d + n3d; i++) {
          o_nnz[localIndex + i] -= (n3d + n2d);
        }

        for (int z = 1; z < mesh->LocalNz; z++) {
          int localIndex = ROUND(index(it.ind, mesh->yend, z));

          // Only 3D fields
          for (int i = 0; i < n3d; i++) {
            o_nnz[localIndex + i] -= (n3d + n2d);
          }
        }
      }

      output_progress.write("Pre-allocating Jacobian\n");

      // Pre-allocate
      MatMPIAIJSetPreallocation(Jmf, 0, d_nnz.data(), 0, o_nnz.data());
      MatSetUp(Jmf);
      MatSetOption(Jmf, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

      // Determine which row/columns of the matrix are locally owned
      int Istart, Iend;
      MatGetOwnershipRange(Jmf, &Istart, &Iend);

      // Convert local into global indices
      index += Istart;

      // Now communicate to fill guard cells
      mesh->communicate(index);

      //////////////////////////////////////////////////
      // Mark non-zero entries

      output_progress.write("Marking non-zero Jacobian entries\n");

      // Offsets for a 5-point pattern
      constexpr std::size_t stencil_size = 5;
      const std::array<int, stencil_size> xoffset = {0, -1, 1, 0, 0};
      const std::array<int, stencil_size> yoffset = {0, 0, 0, -1, 1};

      PetscScalar val = 1.0;

      for (int x = mesh->xstart; x <= mesh->xend; x++) {
        for (int y = mesh->ystart; y <= mesh->yend; y++) {

          int ind0 = ROUND(index(x, y, 0));

          // 2D fields
          for (int i = 0; i < n2d; i++) {
            PetscInt row = ind0 + i;

            // Loop through each point in the 5-point stencil
            for (std::size_t c = 0; c < stencil_size; c++) {
              int xi = x + xoffset[c];
              int yi = y + yoffset[c];

              if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx)
                  || (yi >= mesh->LocalNy)) {
                continue;
              }

              int ind2 = ROUND(index(xi, yi, 0));

              if (ind2 < 0) {
                continue; // A boundary point
              }

              // Depends on all variables on this cell
              for (int j = 0; j < n2d; j++) {
                PetscInt col = ind2 + j;

                MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
              }
            }
          }

          // 3D fields
          for (int z = 0; z < mesh->LocalNz; z++) {

            int ind = ROUND(index(x, y, z));

            for (int i = 0; i < n3d; i++) {
              PetscInt row = ind + i;
              if (z == 0) {
                row += n2d;
              }

              // Depends on 2D fields
              for (int j = 0; j < n2d; j++) {
                PetscInt col = ind0 + j;
                MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
              }

              // 5 point star pattern
              for (std::size_t c = 0; c < stencil_size; c++) {
                int xi = x + xoffset[c];
                int yi = y + yoffset[c];

                if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx)
                    || (yi >= mesh->LocalNy)) {
                  continue;
                }

                int ind2 = ROUND(index(xi, yi, z));
                if (ind2 < 0) {
                  continue; // Boundary point
                }

                if (z == 0) {
                  ind2 += n2d;
                }

                // 3D fields on this cell
                for (int j = 0; j < n3d; j++) {
                  PetscInt col = ind2 + j;
                  MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
                }
              }

              int nz = mesh->LocalNz;
              if (nz > 1) {
                // Multiple points in z

                int zp = (z + 1) % nz;

                int ind2 = ROUND(index(x, y, zp));
                if (zp == 0) {
                  ind2 += n2d;
                }
                for (int j = 0; j < n3d; j++) {
                  PetscInt col = ind2 + j;
                  MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
                }

                int zm = (z - 1 + nz) % nz;
                ind2 = ROUND(index(x, y, zm));
                if (zm == 0) {
                  ind2 += n2d;
                }
                for (int j = 0; j < n3d; j++) {
                  PetscInt col = ind2 + j;
                  MatSetValues(Jmf, 1, &row, 1, &col, &val, INSERT_VALUES);
                }
              }
            }
          }
        }
      }
      // Finished marking non-zero entries

      output_progress.write("Assembling Jacobian matrix\n");

      // Assemble Matrix
      MatAssemblyBegin(Jmf, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(Jmf, MAT_FINAL_ASSEMBLY);

      output_progress.write("Creating Jacobian coloring\n");

      ISColoring iscoloring;

      MatColoring coloring; // This new in PETSc 3.5
      MatColoringCreate(Jmf, &coloring);
      MatColoringSetType(coloring, MATCOLORINGSL);
      MatColoringSetFromOptions(coloring);
      // Calculate index sets
      MatColoringApply(coloring, &iscoloring);
      MatColoringDestroy(&coloring);

      // Create data structure for SNESComputeJacobianDefaultColor
      MatFDColoringCreate(Jmf, iscoloring, &fdcoloring);
      // Set the function to difference
      MatFDColoringSetFunction(
          fdcoloring, reinterpret_cast<PetscErrorCode (*)()>(FormFunctionForColoring),
          this);
      MatFDColoringSetFromOptions(fdcoloring);
      MatFDColoringSetUp(Jmf, iscoloring, fdcoloring);
      ISColoringDestroy(&iscoloring);

      SNESSetJacobian(snes, Jmf, Jmf, SNESComputeJacobianDefaultColor, fdcoloring);
    } else {
      // Brute force calculation

      MatCreateAIJ(
          BoutComm::get(), nlocal, nlocal,  // Local sizes
          PETSC_DETERMINE, PETSC_DETERMINE, // Global sizes
          3, // Number of nonzero entries in diagonal portion of local submatrix
          nullptr,
          0, // Number of nonzeros per row in off-diagonal portion of local submatrix
          nullptr, &Jmf);
#if PETSC_VERSION_GE(3, 4, 0)
      SNESSetJacobian(snes, Jmf, Jmf, SNESComputeJacobianDefault, this);
#else
      // Before 3.4
      SNESSetJacobian(snes, Jmf, Jmf, SNESDefaultComputeJacobian, this);
#endif
      MatSetOption(Jmf, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }

    // Re-use Jacobian
    SNESSetLagJacobian(snes, lag_jacobian);
    // Set Jacobian and preconditioner to persist across time steps
    SNESSetLagJacobianPersists(snes, PETSC_TRUE);
    SNESSetLagPreconditionerPersists(snes, PETSC_TRUE);
    SNESSetLagPreconditioner(snes, 1); // Rebuild when Jacobian is rebuilt
  }

  // Set tolerances
  SNESSetTolerances(snes, atol, rtol, stol, maxits, PETSC_DEFAULT);

  // Force SNES to take at least one nonlinear iteration.
  // This may prevent the solver from getting stuck in false steady state conditions
#if PETSC_VERSION_GE(3, 8, 0)
  SNESSetForceIteration(snes, PETSC_TRUE);
#endif

  // Get KSP context from SNES
  KSP ksp;
  SNESGetKSP(snes, &ksp);

  if (ksp_type != "default") {
    KSPSetType(ksp, ksp_type.c_str());
  }

  if (kspsetinitialguessnonzero) {
    // Set the initial guess to be non-zero
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  }

  KSPSetTolerances(ksp,
                   PETSC_DEFAULT, // rtol
                   PETSC_DEFAULT, // abstol
                   PETSC_DEFAULT, // dtol (divergence tolerance)
                   maxl);         // Maximum number of iterations

  // Get PC context from KSP
  PC pc;
  KSPGetPC(ksp, &pc);

  if (use_precon && hasPreconditioner()) {
    output_info.write("\tUsing user-supplied preconditioner\n");

    // Set a Shell (matrix-free) preconditioner type
    PCSetType(pc, PCSHELL);

    // Specify the preconditioner function
    PCShellSetApply(pc, snesPCapply);
    // Context used to supply object pointer
    PCShellSetContext(pc, this);
  } else if (matrix_free) {
    // Can't use preconditioner because no Jacobian matrix available
    PCSetType(pc, PCNONE);
  } else {
    // Set PC type from input
    if (pc_type != "default") {
      PCSetType(pc, pc_type.c_str());
    }
  }

  // Get runtime options
  lib.setOptionsFromInputFile(snes);

  {
    // Some reporting
    PCType pctype;
    PCGetType(pc, &pctype);
    KSPType ksptype;
    KSPGetType(ksp, &ksptype);
    SNESType snestype;
    SNESGetType(snes, &snestype);
    output_info.write("SNES Type : {}\n", snestype);
    if (ksptype) {
      output_info.write("KSP Type  : {}\n", ksptype);
    }
    if (pctype) {
      output_info.write("PC Type   : {}\n", pctype);
    }
  }

  return 0;
}

int SNESSolver::run() {
  TRACE("SNESSolver::run()");
  // Set initial guess at the solution from variables
  {
    BoutReal* xdata = nullptr;
    int ierr = VecGetArray(snes_x, &xdata);
    CHKERRQ(ierr);
    save_vars(xdata);
    ierr = VecRestoreArray(snes_x, &xdata);
    CHKERRQ(ierr);
  }

  for (int s = 0; s < getNumberOutputSteps(); s++) {
    BoutReal target = simtime + getOutputTimestep();

    bool looping = true;
    int snes_failures = 0; // Count SNES convergence failures
    int saved_jacobian_lag = 0;
    do {
      // Copy the state (snes_x) into initial values (x0)
      VecCopy(snes_x, x0);

      if (timestep < dt_min_reset) {
        // Hit the minimum timestep, probably through repeated failures

        if (saved_jacobian_lag != 0) {
          // Already tried this and it didn't work
          throw BoutException("Solver failed after many attempts");
        }

        // Try resetting the preconditioner, turn off predictor, and use a large timestep
        SNESGetLagJacobian(snes, &saved_jacobian_lag);
        SNESSetLagJacobian(snes, 1);
        timestep = getOutputTimestep();
        predictor = false; // Predictor can cause problems in near steady-state.
      }

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
      PetscErrorCode ierr = SNESSolve(snes, nullptr, snes_x);

      // Find out if converged
      SNESConvergedReason reason;
      SNESGetConvergedReason(snes, &reason);
      if ((ierr != 0) or (reason < 0)) {
        // Diverged or SNES failed

        if (diagnose_failures) {
          // Print diagnostics to help identify source of the problem

          output.write("\n======== SNES failed =========\n");
          output.write("\nReturn code: {}, reason: {}\n", ierr, reason);
          for (const auto& f : f2d) {
            output.write(
                "{} : ({} -> {}), ddt: ({} -> {})\n", f.name,
                min(*f.var, true, "RGN_NOBNDRY"), max(*f.var, true, "RGN_NOBNDRY"),
                min(*f.F_var, true, "RGN_NOBNDRY"), max(*f.F_var, true, "RGN_NOBNDRY"));
          }
          for (const auto& f : f3d) {
            output.write(
                "{} : ({} -> {}), ddt: ({} -> {})\n", f.name,
                min(*f.var, true, "RGN_NOBNDRY"), max(*f.var, true, "RGN_NOBNDRY"),
                min(*f.F_var, true, "RGN_NOBNDRY"), max(*f.F_var, true, "RGN_NOBNDRY"));
          }
        }

        ++snes_failures;

        // Try a smaller timestep
        timestep /= 2.0;
        // Restore state
        VecCopy(x0, snes_x);

        // Check lock state
        PetscInt lock_state;
        VecLockGet(snes_x, &lock_state);
        if (lock_state > 0) {
          // Locked for read
          output_warn.write("WARNING: snes_x locked for reading\n");
        } else if (lock_state < 0) {
          // Locked for write
          output_warn.write("WARNING: snes_x locked for writing\n");
        }
        looping = true;
        continue; // Try again
      }

      if (saved_jacobian_lag != 0) {
        // Following successful step, reset Jacobian lag
        // to previous value
        SNESSetLagJacobian(snes, saved_jacobian_lag);
        saved_jacobian_lag = 0;
      }

      if (predictor) {
        // Save previous values: x1 <- x0
        VecCopy(x0, x1);
        time1 = simtime;
      }

      int nl_its;
      SNESGetIterationNumber(snes, &nl_its);

      if (nl_its == 0) {
        // This can occur even with SNESSetForceIteration
        // Results in simulation state freezing and rapidly going to the end

        {
          const BoutReal* xdata = nullptr;
          int ierr = VecGetArrayRead(snes_x, &xdata);
          CHKERRQ(ierr);
          load_vars(const_cast<BoutReal*>(xdata));
          ierr = VecRestoreArrayRead(snes_x, &xdata);
          CHKERRQ(ierr);
        }
        run_rhs(simtime);

        // Copy derivatives back
        {
          BoutReal* fdata = nullptr;
          ierr = VecGetArray(snes_f, &fdata);
          CHKERRQ(ierr);
          save_derivs(fdata);
          ierr = VecRestoreArray(snes_f, &fdata);
          CHKERRQ(ierr);
        }

        // Forward Euler
        VecAXPY(snes_x, dt, snes_f);
      }

      simtime += dt;

      if (diagnose) {
        // Gather and print diagnostic information

        int lin_its;
        SNESGetLinearSolveIterations(snes, &lin_its);

        output.print("\r"); // Carriage return for printing to screen
        output.write("Time: {}, timestep: {}, nl iter: {}, lin iter: {}, reason: {}",
                     simtime, timestep, nl_its, lin_its, reason);
        if (snes_failures > 0) {
          output.write(", SNES failures: {}", snes_failures);
          snes_failures = 0;
        }
        output.write("\n");
      }

      if (looping) {
        if (nl_its <= lower_its) {
          // Increase timestep slightly
          timestep *= 1.1;

          if (timestep > max_timestep) {
            timestep = max_timestep;
          }
        } else if (nl_its >= upper_its) {
          // Reduce timestep slightly
          timestep *= 0.9;
        }
      }
    } while (looping);

    // Put the result into variables
    {
      const BoutReal* xdata = nullptr;
      int ierr = VecGetArrayRead(snes_x, &xdata);
      CHKERRQ(ierr);
      load_vars(const_cast<BoutReal*>(xdata));
      ierr = VecRestoreArrayRead(snes_x, &xdata);
      CHKERRQ(ierr);
    }
    run_rhs(simtime); // Run RHS to calculate auxilliary variables

    if (call_monitors(simtime, s, getNumberOutputSteps()) != 0) {
      break; // User signalled to quit
    }
  }

  return 0;
}

// f = rhs
PetscErrorCode SNESSolver::snes_function(Vec x, Vec f, bool linear) {
  // Get data from PETSc into BOUT++ fields
  const BoutReal* xdata = nullptr;
  int ierr = VecGetArrayRead(x, &xdata);
  CHKERRQ(ierr);
  load_vars(const_cast<BoutReal*>(xdata));
  ierr = VecRestoreArrayRead(x, &xdata);
  CHKERRQ(ierr);

  try {
    // Call RHS function
    run_rhs(simtime + dt, linear);
  } catch (BoutException& e) {
    // Simulation might fail, e.g. negative densities
    // if timestep too large
    output_warn.write("WARNING: BoutException thrown: {}\n", e.what());

    // Tell SNES that the input was out of domain
    SNESSetFunctionDomainError(snes);
    // Note: Returning non-zero error here leaves vectors in locked state
    return 0;
  }

  // Copy derivatives back
  BoutReal* fdata = nullptr;
  ierr = VecGetArray(f, &fdata);
  CHKERRQ(ierr);
  save_derivs(fdata);
  ierr = VecRestoreArray(f, &fdata);
  CHKERRQ(ierr);

  switch (equation_form) {
  case BoutSnesEquationForm::pseudo_transient: {
    // Pseudo-transient timestepping (as in UEDGE)
    // f <- f - x/Δt
    VecAXPY(f, -1. / dt, x);
    break;
  }
  case BoutSnesEquationForm::rearranged_backward_euler: {
    // Rearranged Backward Euler
    // f = (x0 - x)/Δt + f
    // First calculate x - x0 to minimise floating point issues
    VecWAXPY(delta_x, -1.0, x0, x); // delta_x = x - x0
    VecAXPY(f, -1. / dt, delta_x);  // f <- f - delta_x / dt
    break;
  }
  case BoutSnesEquationForm::backward_euler: {
    // Backward Euler
    // Set f = x - x0 - Δt*f
    VecAYPX(f, -dt, x);   // f <- x - Δt*f
    VecAXPY(f, -1.0, x0); // f <- f - x0
    break;
  }
  case BoutSnesEquationForm::direct_newton: {
    // Direct Newton solve -> don't modify f
    break;
  }
  default: {
    throw BoutException("Invalid choice of equation_form. Try 0-3");
  }
  };

  return 0;
}

/*
 * Preconditioner function
 */
PetscErrorCode SNESSolver::precon(Vec x, Vec f) {
  if (!hasPreconditioner()) {
    // No user preconditioner
    throw BoutException("No user preconditioner");
  }

  int ierr;

  // Get data from PETSc into BOUT++ fields
  Vec solution;
  SNESGetSolution(snes, &solution);
  BoutReal* soldata;
  ierr = VecGetArray(solution, &soldata);
  CHKERRQ(ierr);
  load_vars(soldata);
  ierr = VecRestoreArray(solution, &soldata);
  CHKERRQ(ierr);

  // Load vector to be inverted into ddt() variables
  const BoutReal* xdata;
  ierr = VecGetArrayRead(x, &xdata);
  CHKERRQ(ierr);
  load_derivs(const_cast<BoutReal*>(xdata)); // Note: load_derivs does not modify data
  ierr = VecRestoreArrayRead(x, &xdata);
  CHKERRQ(ierr);

  // Run the preconditioner
  runPreconditioner(simtime + dt, dt, 0.0);

  // Save the solution from F_vars
  BoutReal* fdata;
  ierr = VecGetArray(f, &fdata);
  CHKERRQ(ierr);
  save_derivs(fdata);
  ierr = VecRestoreArray(f, &fdata);
  CHKERRQ(ierr);

  return 0;
}

#endif // BOUT_HAS_PETSC
