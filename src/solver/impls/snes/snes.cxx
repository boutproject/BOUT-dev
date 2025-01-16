#include "bout/build_defines.hxx"

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

#include "petscmat.h"
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
      pc_hypre_type((*options)["pc_hypre_type"]
                        .doc("hypre preconditioner type: euclid, pilut, parasails, "
                             "boomeramg, ams, ads")
                        .withDefault("pilut")),
      line_search_type((*options)["line_search_type"]
                           .doc("Line search type: basic, bt, l2, cp, nleqerr")
                           .withDefault("default")),
      matrix_free((*options)["matrix_free"]
                      .doc("Use matrix free Jacobian?")
                      .withDefault<bool>(false)),
      matrix_free_operator((*options)["matrix_free_operator"]
                               .doc("Use matrix free Jacobian-vector operator?")
                               .withDefault<bool>(true)),
      lag_jacobian((*options)["lag_jacobian"]
                       .doc("Re-use the Jacobian this number of SNES iterations")
                       .withDefault(50)),
      use_coloring((*options)["use_coloring"]
                       .doc("Use matrix coloring to calculate Jacobian?")
                       .withDefault<bool>(true)),
      jacobian_recalculated(false),
      prune_jacobian((*options)["prune_jacobian"]
                         .doc("Remove small elements in the Jacobian?")
                         .withDefault<bool>(false)),
      prune_abstol((*options)["prune_abstol"]
                       .doc("Prune values with absolute values smaller than this")
                       .withDefault<BoutReal>(1e-16)),
      prune_fraction((*options)["prune_fraction"]
                         .doc("Prune if fraction of small elements is larger than this")
                         .withDefault<BoutReal>(0.2)),
      scale_rhs((*options)["scale_rhs"]
                    .doc("Scale time derivatives (Jacobian row scaling)?")
                    .withDefault<bool>(false)),
      scale_vars((*options)["scale_vars"]
                     .doc("Scale variables (Jacobian column scaling)?")
                     .withDefault<bool>(false)) {}

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
  output_info.write("Creating vector\n");
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
    ierr = VecDuplicate(snes_x, &delta_x);
    CHKERRQ(ierr);
  }

  if (predictor) {
    // Storage for previous solution
    ierr = VecDuplicate(snes_x, &x1);
    CHKERRQ(ierr);
  }

  if (scale_rhs) {
    // Storage for rhs factors, one per evolving variable
    ierr = VecDuplicate(snes_x, &rhs_scaling_factors);
    CHKERRQ(ierr);
    // Set all factors to 1 to start with
    ierr = VecSet(rhs_scaling_factors, 1.0);
    CHKERRQ(ierr);
    // Array to store inverse Jacobian row norms
    ierr = VecDuplicate(snes_x, &jac_row_inv_norms);
    CHKERRQ(ierr);
  }

  if (scale_vars) {
    // Storage for var factors, one per evolving variable
    ierr = VecDuplicate(snes_x, &var_scaling_factors);
    CHKERRQ(ierr);
    // Set all factors to 1 to start with
    ierr = VecSet(var_scaling_factors, 1.0);
    CHKERRQ(ierr);
    // Storage for scaled 'x' state vectors
    ierr = VecDuplicate(snes_x, &scaled_x);
    CHKERRQ(ierr);
  }

  // Nonlinear solver interface (SNES)
  output_info.write("Create SNES\n");
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
  if (matrix_free or matrix_free_operator) {
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
  }

  if (matrix_free) {
    // Use matrix free for both operator and preconditioner
    // Calculate Jacobian matrix free using FormFunctionForDifferencing
    SNESSetJacobian(snes, Jmf, Jmf, MatMFFDComputeJacobian, this);

  } else {
    // Calculate the Jacobian using finite differences.
    // The finite difference Jacobian (Jfd) may be used for both operator
    // and preconditioner or, if matrix_free_operator, in only the preconditioner.
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

      int n2d = f2d.size();
      int n3d = f3d.size();

      // Set size of Matrix on each processor to nlocal x nlocal
      MatCreate(BoutComm::get(), &Jfd);
      MatSetSizes(Jfd, nlocal, nlocal, PETSC_DETERMINE, PETSC_DETERMINE);
      MatSetFromOptions(Jfd);

      std::vector<PetscInt> d_nnz(nlocal);
      std::vector<PetscInt> o_nnz(nlocal);

      // Set values for most points
      const int ncells_x = (mesh->LocalNx > 1) ? 2 : 0;
      const int ncells_y = (mesh->LocalNy > 1) ? 2 : 0;
      const int ncells_z = (mesh->LocalNz > 1) ? 2 : 0;

      const auto star_pattern = (1 + ncells_x + ncells_y) * (n3d + n2d) + ncells_z * n3d;

      // Offsets. Start with the central cell
      std::vector<std::pair<int, int>> xyoffsets{{0, 0}};
      if (ncells_x != 0) {
        // Stencil includes points in X
        xyoffsets.push_back({-1, 0});
        xyoffsets.push_back({1, 0});
      }
      if (ncells_y != 0) {
        // Stencil includes points in Y
        xyoffsets.push_back({0, -1});
        xyoffsets.push_back({0, 1});
      }

      output_info.write("Star pattern: {} non-zero entries\n", star_pattern);
      for (int i = 0; i < nlocal; i++) {
        // Non-zero elements on this processor
        d_nnz[i] = star_pattern;
        // Non-zero elements on neighboring processor
        o_nnz[i] = 0;
      }

      // X boundaries
      if (ncells_x != 0) {
        if (mesh->firstX()) {
          // Lower X boundary
          for (int y = mesh->ystart; y <= mesh->yend; y++) {
            for (int z = 0; z < mesh->LocalNz; z++) {
              const int localIndex = ROUND(index(mesh->xstart, y, z));
              ASSERT2((localIndex >= 0) && (localIndex < nlocal));
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
              const int localIndex = ROUND(index(mesh->xstart, y, z));
              ASSERT2((localIndex >= 0) && (localIndex < nlocal));
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
              const int localIndex = ROUND(index(mesh->xend, y, z));
              ASSERT2((localIndex >= 0) && (localIndex < nlocal));
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
              const int localIndex = ROUND(index(mesh->xend, y, z));
              ASSERT2((localIndex >= 0) && (localIndex < nlocal));
              const int num_fields = (z == 0) ? n2d + n3d : n3d;
              for (int i = 0; i < num_fields; i++) {
                d_nnz[localIndex + i] -= (n3d + n2d);
                o_nnz[localIndex + i] += (n3d + n2d);
              }
            }
          }
        }
      }

      // Y boundaries
      if (ncells_y != 0) {
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
            d_nnz[localIndex + i] -= (n3d + n2d);
          }

          for (int z = 1; z < mesh->LocalNz; z++) {
            localIndex = ROUND(index(x, mesh->ystart, z));

            // Only 3D fields
            for (int i = 0; i < n3d; i++) {
              o_nnz[localIndex + i] += (n3d + n2d);
              d_nnz[localIndex + i] -= (n3d + n2d);
            }
          }

          // z = 0 case
          localIndex = ROUND(index(x, mesh->yend, 0));
          // All 2D and 3D fields
          for (int i = 0; i < n2d + n3d; i++) {
            o_nnz[localIndex + i] += (n3d + n2d);
            d_nnz[localIndex + i] -= (n3d + n2d);
          }

          for (int z = 1; z < mesh->LocalNz; z++) {
            localIndex = ROUND(index(x, mesh->yend, z));

            // Only 3D fields
            for (int i = 0; i < n3d; i++) {
              o_nnz[localIndex + i] += (n3d + n2d);
              d_nnz[localIndex + i] -= (n3d + n2d);
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
          const int localIndex = ROUND(index(it.ind, mesh->yend, 0));
          if (localIndex < 0) {
            continue; // Out of domain
          }

          // All 2D and 3D fields
          for (int i = 0; i < n2d + n3d; i++) {
            o_nnz[localIndex + i] -= (n3d + n2d);
          }

          for (int z = 1; z < mesh->LocalNz; z++) {
            const int localIndex = ROUND(index(it.ind, mesh->yend, z));

            // Only 3D fields
            for (int i = 0; i < n3d; i++) {
              o_nnz[localIndex + i] -= (n3d + n2d);
            }
          }
        }
      }

      output_progress.write("Pre-allocating Jacobian\n");

      // Pre-allocate
      MatMPIAIJSetPreallocation(Jfd, 0, d_nnz.data(), 0, o_nnz.data());
      MatSeqAIJSetPreallocation(Jfd, 0, d_nnz.data());
      MatSetUp(Jfd);
      MatSetOption(Jfd, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

      // Determine which row/columns of the matrix are locally owned
      int Istart, Iend;
      MatGetOwnershipRange(Jfd, &Istart, &Iend);

      // Convert local into global indices
      // Note: Not in the boundary cells, to keep -1 values
      for (const auto& i : mesh->getRegion3D("RGN_NOBNDRY")) {
        index[i] += Istart;
      }

      // Now communicate to fill guard cells
      mesh->communicate(index);

      //////////////////////////////////////////////////
      // Mark non-zero entries

      output_progress.write("Marking non-zero Jacobian entries\n");

      const PetscScalar val = 1.0;

      for (int x = mesh->xstart; x <= mesh->xend; x++) {
        for (int y = mesh->ystart; y <= mesh->yend; y++) {

          const int ind0 = ROUND(index(x, y, 0));

          // 2D fields
          for (int i = 0; i < n2d; i++) {
            const PetscInt row = ind0 + i;

            // Loop through each point in the 5-point stencil
            for (const auto& xyoffset : xyoffsets) {
              const int xi = x + xyoffset.first;
              const int yi = y + xyoffset.second;

              if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx)
                  || (yi >= mesh->LocalNy)) {
                continue;
              }

              const int ind2 = ROUND(index(xi, yi, 0));

              if (ind2 < 0) {
                continue; // A boundary point
              }

              // Depends on all variables on this cell
              for (int j = 0; j < n2d; j++) {
                const PetscInt col = ind2 + j;
                ierr = MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES);
                CHKERRQ(ierr);
              }
            }
          }

          // 3D fields
          for (int z = 0; z < mesh->LocalNz; z++) {

            const int ind = ROUND(index(x, y, z));

            for (int i = 0; i < n3d; i++) {
              PetscInt row = ind + i;
              if (z == 0) {
                row += n2d;
              }

              // Depends on 2D fields
              for (int j = 0; j < n2d; j++) {
                const PetscInt col = ind0 + j;
                ierr = MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES);
                CHKERRQ(ierr);
              }

              // Star pattern
              for (const auto& xyoffset : xyoffsets) {
                const int xi = x + xyoffset.first;
                const int yi = y + xyoffset.second;

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
                  const PetscInt col = ind2 + j;
                  ierr = MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES);
                  if (ierr != 0) {
                    output.write("ERROR: {} : ({}, {}) -> ({}, {}) : {} -> {}\n", row, x,
                                 y, xi, yi, ind2, ind2 + n3d - 1);
                  }
                  CHKERRQ(ierr);
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
                  const PetscInt col = ind2 + j;
                  ierr = MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES);
                  CHKERRQ(ierr);
                }

                int zm = (z - 1 + nz) % nz;
                ind2 = ROUND(index(x, y, zm));
                if (zm == 0) {
                  ind2 += n2d;
                }
                for (int j = 0; j < n3d; j++) {
                  const PetscInt col = ind2 + j;
                  ierr = MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES);
                  CHKERRQ(ierr);
                }
              }
            }
          }
        }
      }
      // Finished marking non-zero entries

      output_progress.write("Assembling Jacobian matrix\n");

      // Assemble Matrix
      MatAssemblyBegin(Jfd, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(Jfd, MAT_FINAL_ASSEMBLY);

      output_progress.write("Creating Jacobian coloring\n");
      updateColoring();

      if (prune_jacobian) {
        // Will remove small elements from the Jacobian.
        // Save a copy to recover from over-pruning
        ierr = MatDuplicate(Jfd, MAT_SHARE_NONZERO_PATTERN, &Jfd_original);
        CHKERRQ(ierr);
      }
    } else {
      // Brute force calculation
      // There is usually no reason to use this, except as a check of
      // the coloring calculation.

      MatCreateAIJ(
          BoutComm::get(), nlocal, nlocal,  // Local sizes
          PETSC_DETERMINE, PETSC_DETERMINE, // Global sizes
          3, // Number of nonzero entries in diagonal portion of local submatrix
          nullptr,
          0, // Number of nonzeros per row in off-diagonal portion of local submatrix
          nullptr, &Jfd);

      if (matrix_free_operator) {
        SNESSetJacobian(snes, Jmf, Jfd, SNESComputeJacobianDefault, this);
      } else {
        SNESSetJacobian(snes, Jfd, Jfd, SNESComputeJacobianDefault, this);
      }

      MatSetOption(Jfd, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }

    // Re-use Jacobian
    // Note: If the 'Amat' Jacobian is matrix free, SNESComputeJacobian
    //       always updates its reference 'u' vector every nonlinear iteration
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

      if (pc_type == "hypre") {
#if PETSC_HAVE_HYPRE
        // Set the type of hypre preconditioner
        PCHYPRESetType(pc, pc_hypre_type.c_str());
#else
        throw BoutException("PETSc was not configured with Hypre.");
#endif
      }
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
  int ierr;
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
    int loop_count = 0;
    do {
      if (scale_vars) {
        // Individual variable scaling
        // Note: If variables are rescaled then the Jacobian columns
        //       need to be scaled or recalculated

        if (loop_count % 100 == 0) {
          // Rescale state (snes_x) so that all quantities are around 1
          // If quantities are near zero then RTOL is used
          int istart, iend;
          VecGetOwnershipRange(snes_x, &istart, &iend);

          // Take ownership of snes_x and var_scaling_factors data
          PetscScalar* snes_x_data = nullptr;
          ierr = VecGetArray(snes_x, &snes_x_data);
          CHKERRQ(ierr);
          PetscScalar* x1_data;
          ierr = VecGetArray(x1, &x1_data);
          CHKERRQ(ierr);
          PetscScalar* var_scaling_factors_data;
          ierr = VecGetArray(var_scaling_factors, &var_scaling_factors_data);
          CHKERRQ(ierr);

          // Normalise each value in the state
          // Limit normalisation so scaling factor is never smaller than rtol
          for (int i = 0; i < iend - istart; ++i) {
            const PetscScalar norm =
                BOUTMAX(std::abs(snes_x_data[i]), rtol / var_scaling_factors_data[i]);
            snes_x_data[i] /= norm;
            x1_data[i] /= norm; // Update history for predictor
            var_scaling_factors_data[i] *= norm;
          }

          // Restore vector underlying data
          ierr = VecRestoreArray(var_scaling_factors, &var_scaling_factors_data);
          CHKERRQ(ierr);
          ierr = VecRestoreArray(x1, &x1_data);
          CHKERRQ(ierr);
          ierr = VecRestoreArray(snes_x, &snes_x_data);
          CHKERRQ(ierr);

          if (diagnose) {
            // Print maximum and minimum scaling factors
            PetscReal max_scale, min_scale;
            VecMax(var_scaling_factors, nullptr, &max_scale);
            VecMin(var_scaling_factors, nullptr, &min_scale);
            output.write("Var scaling: {} -> {}\n", min_scale, max_scale);
          }

          // Force recalculation of the Jacobian
          SNESGetLagJacobian(snes, &saved_jacobian_lag);
          SNESSetLagJacobian(snes, 1);
        }
      }
      ++loop_count;

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

      // Get number of iterations
      int nl_its;
      SNESGetIterationNumber(snes, &nl_its);
      int lin_its;
      SNESGetLinearSolveIterations(snes, &lin_its);

      if ((ierr != 0) or (reason < 0)) {
        // Diverged or SNES failed

        if (diagnose_failures) {
          // Print diagnostics to help identify source of the problem

          output.write("\n======== SNES failed =========\n");
          output.write("\nReturn code: {}, reason: {}\n", ierr, static_cast<int>(reason));
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

        // Recalculate the Jacobian
        if (jacobian_pruned and (snes_failures > 2) and (4 * lin_its > 3 * maxl)) {
          // Taking 3/4 of maximum linear iterations on average per linear step
          // May indicate a preconditioner problem.
          // Restore pruned non-zero elements
          if (diagnose) {
            output.write("\nRestoring Jacobian\n");
          }
          ierr = MatCopy(Jfd_original, Jfd, DIFFERENT_NONZERO_PATTERN);
          CHKERRQ(ierr);
          // The non-zero pattern has changed, so update coloring
          updateColoring();
          jacobian_pruned = false; // Reset flag. Will be set after pruning.
        }
        if (saved_jacobian_lag == 0) {
          SNESGetLagJacobian(snes, &saved_jacobian_lag);
          SNESSetLagJacobian(snes, 1);
        }

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

      if (nl_its == 0) {
        // This can occur even with SNESSetForceIteration
        // Results in simulation state freezing and rapidly going to the end

        if (scale_vars) {
          // scaled_x <- snes_x * var_scaling_factors
          ierr = VecPointwiseMult(scaled_x, snes_x, var_scaling_factors);
          CHKERRQ(ierr);

          const BoutReal* xdata = nullptr;
          ierr = VecGetArrayRead(scaled_x, &xdata);
          CHKERRQ(ierr);
          load_vars(const_cast<BoutReal*>(xdata));
          ierr = VecRestoreArrayRead(scaled_x, &xdata);
          CHKERRQ(ierr);
        } else {
          const BoutReal* xdata = nullptr;
          ierr = VecGetArrayRead(snes_x, &xdata);
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

        output.print("\r"); // Carriage return for printing to screen
        output.write("Time: {}, timestep: {}, nl iter: {}, lin iter: {}, reason: {}",
                     simtime, timestep, nl_its, lin_its, static_cast<int>(reason));
        if (snes_failures > 0) {
          output.write(", SNES failures: {}", snes_failures);
          snes_failures = 0;
        }
        output.write("\n");
      }

#if PETSC_VERSION_GE(3, 20, 0)
      // MatFilter and MatEliminateZeros(Mat, bool) require PETSc >= 3.20
      if (jacobian_recalculated and prune_jacobian) {
        jacobian_recalculated = false; // Reset flag

        // Remove small elements from the Jacobian and recompute the coloring
        // Only do this if there are a significant number of small elements.
        int small_elements = 0;
        int total_elements = 0;

        // Get index of rows owned by this processor
        int rstart, rend;
        MatGetOwnershipRange(Jfd, &rstart, &rend);

        PetscInt ncols;
        const PetscScalar* vals;
        for (int row = rstart; row < rend; row++) {
          MatGetRow(Jfd, row, &ncols, nullptr, &vals);
          for (int col = 0; col < ncols; col++) {
            if (std::abs(vals[col]) < prune_abstol) {
              ++small_elements;
            }
            ++total_elements;
          }
          MatRestoreRow(Jfd, row, &ncols, nullptr, &vals);
        }

        if (small_elements > prune_fraction * total_elements) {
          if (diagnose) {
            output.write("\nPruning Jacobian elements: {} / {}\n", small_elements,
                         total_elements);
          }

          // Prune Jacobian, keeping diagonal elements
          ierr = MatFilter(Jfd, prune_abstol, PETSC_TRUE, PETSC_TRUE);

          // Update the coloring from Jfd matrix
          updateColoring();

          // Mark the Jacobian as pruned. This is so that it is only restored if pruned.
          jacobian_pruned = true;
        }
      }
#endif // PETSC_VERSION_GE(3,20,0)

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
    if (scale_vars) {
      // scaled_x <- snes_x * var_scaling_factors
      int ierr = VecPointwiseMult(scaled_x, snes_x, var_scaling_factors);
      CHKERRQ(ierr);

      const BoutReal* xdata = nullptr;
      ierr = VecGetArrayRead(scaled_x, &xdata);
      CHKERRQ(ierr);
      load_vars(const_cast<BoutReal*>(xdata));
      ierr = VecRestoreArrayRead(scaled_x, &xdata);
      CHKERRQ(ierr);
    } else {
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
  if (scale_vars) {
    // scaled_x <- x * var_scaling_factors
    int ierr = VecPointwiseMult(scaled_x, x, var_scaling_factors);
    CHKERRQ(ierr);

    const BoutReal* xdata = nullptr;
    ierr = VecGetArrayRead(scaled_x, &xdata);
    CHKERRQ(ierr);
    load_vars(const_cast<BoutReal*>(
        xdata)); // const_cast needed due to load_vars API. Not writing to xdata.
    ierr = VecRestoreArrayRead(scaled_x, &xdata);
    CHKERRQ(ierr);
  } else {
    const BoutReal* xdata = nullptr;
    int ierr = VecGetArrayRead(x, &xdata);
    CHKERRQ(ierr);
    load_vars(const_cast<BoutReal*>(
        xdata)); // const_cast needed due to load_vars API. Not writing to xdata.
    ierr = VecRestoreArrayRead(x, &xdata);
    CHKERRQ(ierr);
  }

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
  int ierr = VecGetArray(f, &fdata);
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

  if (scale_rhs) {
    // f <- f * rhs_scaling_factors
    ierr = VecPointwiseMult(f, f, rhs_scaling_factors);
    CHKERRQ(ierr);
  }

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

PetscErrorCode SNESSolver::scaleJacobian(Mat Jac_new) {
  jacobian_recalculated = true;

  if (!scale_rhs) {
    return 0; // Not scaling the RHS values
  }

  int ierr;

  // Get index of rows owned by this processor
  int rstart, rend;
  MatGetOwnershipRange(Jac_new, &rstart, &rend);

  // Check that the vector has the same ownership range
  int istart, iend;
  VecGetOwnershipRange(jac_row_inv_norms, &istart, &iend);
  if ((rstart != istart) or (rend != iend)) {
    throw BoutException("Ownership ranges different: [{}, {}) and [{}, {})\n", rstart,
                        rend, istart, iend);
  }

  // Calculate the norm of each row of the Jacobian
  PetscScalar* row_inv_norm_data;
  ierr = VecGetArray(jac_row_inv_norms, &row_inv_norm_data);
  CHKERRQ(ierr);

  PetscInt ncols;
  const PetscScalar* vals;
  for (int row = rstart; row < rend; ++row) {
    MatGetRow(Jac_new, row, &ncols, nullptr, &vals);

    // Calculate a norm of this row of the Jacobian
    PetscScalar norm = 0.0;
    for (int col = 0; col < ncols; col++) {
      PetscScalar absval = std::abs(vals[col]);
      if (absval > norm) {
        norm = absval;
      }
      // Can we identify small elements and remove them?
      // so we don't need to calculate them next time
    }

    // Store in the vector as 1 / norm
    row_inv_norm_data[row - rstart] = 1. / norm;

    MatRestoreRow(Jac_new, row, &ncols, nullptr, &vals);
  }

  ierr = VecRestoreArray(jac_row_inv_norms, &row_inv_norm_data);
  CHKERRQ(ierr);

  // Modify the RHS scaling: factor = factor / norm
  ierr = VecPointwiseMult(rhs_scaling_factors, rhs_scaling_factors, jac_row_inv_norms);
  CHKERRQ(ierr);

  if (diagnose) {
    // Print maximum and minimum scaling factors
    PetscReal max_scale, min_scale;
    VecMax(rhs_scaling_factors, nullptr, &max_scale);
    VecMin(rhs_scaling_factors, nullptr, &min_scale);
    output.write("RHS scaling: {} -> {}\n", min_scale, max_scale);
  }

  // Scale the Jacobian rows by multiplying on the left by 1/norm
  ierr = MatDiagonalScale(Jac_new, jac_row_inv_norms, nullptr);
  CHKERRQ(ierr);

  return 0;
}

///
/// Input Parameters:
///   snes - nonlinear solver object
///   x1 - location at which to evaluate Jacobian
///   ctx - MatFDColoring context or NULL
///
/// Output Parameters:
///   Jac - Jacobian matrix (not altered in this routine)
///   Jac_new - newly computed Jacobian matrix to use with preconditioner (generally the same as
///   Jac)
static PetscErrorCode ComputeJacobianScaledColor(SNES snes, Vec x1, Mat Jac, Mat Jac_new,
                                                 void* ctx) {
  PetscErrorCode err = SNESComputeJacobianDefaultColor(snes, x1, Jac, Jac_new, ctx);
  CHKERRQ(err);

  if ((err != 0) or (ctx == nullptr)) {
    return err;
  }

  // Get the the SNESSolver pointer from the function call context
  SNESSolver* fctx = nullptr;
  err = MatFDColoringGetFunction(static_cast<MatFDColoring>(ctx), nullptr,
                                 reinterpret_cast<void**>(&fctx));
  CHKERRQ(err);

  // Call the SNESSolver function
  return fctx->scaleJacobian(Jac_new);
}

void SNESSolver::updateColoring() {
  // Re-calculate the coloring
  MatColoring coloring = NULL;
  MatColoringCreate(Jfd, &coloring);
  MatColoringSetType(coloring, MATCOLORINGSL);
  MatColoringSetFromOptions(coloring);

  // Calculate new index sets
  ISColoring iscoloring = NULL;
  MatColoringApply(coloring, &iscoloring);
  MatColoringDestroy(&coloring);

  // Replace the old coloring with the new one
  MatFDColoringDestroy(&fdcoloring);
  MatFDColoringCreate(Jfd, iscoloring, &fdcoloring);
  MatFDColoringSetFunction(
      fdcoloring, reinterpret_cast<PetscErrorCode (*)()>(FormFunctionForColoring), this);
  MatFDColoringSetFromOptions(fdcoloring);
  MatFDColoringSetUp(Jfd, iscoloring, fdcoloring);
  ISColoringDestroy(&iscoloring);

  // Replace the CTX pointer in SNES Jacobian
  if (matrix_free_operator) {
    // Use matrix-free calculation for operator, finite difference for preconditioner
    SNESSetJacobian(snes, Jmf, Jfd, ComputeJacobianScaledColor, fdcoloring);
  } else {
    SNESSetJacobian(snes, Jfd, Jfd, ComputeJacobianScaledColor, fdcoloring);
  }
}

#endif // BOUT_HAS_PETSC
