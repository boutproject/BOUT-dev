#include "bout/petsc_coloring.hxx"
#include "bout/globals.hxx"

/*!
 * SNES callback for forming Jacobian with coloring
 *
 * This can be a linearised and simplified form of FormFunction
 */
static PetscErrorCode FormFunctionForColoring(SNES UNUSED(snes), Vec x, Vec f,
                                              void* ctx) {
  return static_cast<PetscPreconditioner*>(ctx)->form_function(x, f);
}

PetscPreconditioner::PetscPreconditioner(Options& options, int nlocal,
                                         Field3D local_index, FormFunction form_function)
    : mesh(localIndex.getMesh()), form_function(std::move(form_function)),
      use_coloring((*options)["use_coloring"]
                       .doc("Use matrix coloring to calculate Jacobian?")
                       .withDefault<bool>(true)),
{

  if (use_coloring) {
    // Use matrix coloring
    // This greatly reduces the number of times the rhs() function needs
    // to be evaluated when calculating the Jacobian.

    //////////////////////////////////////////////////
    // Get the local indices by starting at 0
    Field3D index = localIndex;

    //////////////////////////////////////////////////
    // Pre-allocate PETSc storage

    output_progress.write("Setting Jacobian matrix sizes\n");

    const int n2d = f2d.size();
    const int n3d = f3d.size();

    // Set size of Matrix on each processor to nlocal x nlocal
    MatCreate(BoutComm::get(), &Jfd);
    MatSetSizes(Jfd, nlocal, nlocal, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetFromOptions(Jfd);
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

    // Non-zero elements on this processor
    std::vector<PetscInt> d_nnz;
    std::vector<PetscInt> o_nnz;
    auto n_square = (*options)["stencil:square"]
                        .doc("Extent of stencil (square)")
                        .withDefault<int>(0);
    auto n_taxi = (*options)["stencil:taxi"]
                      .doc("Extent of stencil (taxi-cab norm)")
                      .withDefault<int>(0);
    auto n_cross =
        (*options)["stencil:cross"].doc("Extent of stencil (cross)").withDefault<int>(0);
    // Set n_taxi 2 if nothing else is set
    // Probably a better way to do this
    if (n_square == 0 && n_taxi == 0 && n_cross == 0) {
      output_info.write("Setting solver:stencil:taxi = 2\n");
      n_taxi = 2;
    }

    auto const xy_offsets = ColoringStencil::getOffsets(n_square, n_taxi, n_cross);
    {
      // This is ugly but can't think of a better and robust way to
      // count the non-zeros for some arbitrary stencil
      // effectively the same loop as the one that sets the non-zeros below
      std::vector<std::set<int>> d_nnz_map2d(nlocal);
      std::vector<std::set<int>> o_nnz_map2d(nlocal);
      std::vector<std::set<int>> d_nnz_map3d(nlocal);
      std::vector<std::set<int>> o_nnz_map3d(nlocal);
      // Loop over every element in 2D to count the *unique* non-zeros
      for (int x = mesh->xstart; x <= mesh->xend; x++) {
        for (int y = mesh->ystart; y <= mesh->yend; y++) {

          const int ind0 = ROUND(index(x, y, 0)) - Istart;

          // 2D fields
          for (int i = 0; i < n2d; i++) {
            const PetscInt row = ind0 + i;
            // Loop through each point in the stencil
            for (const auto& [x_off, y_off] : xy_offsets) {
              const int xi = x + x_off;
              const int yi = y + y_off;
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
                if (col >= Istart && col < Iend) {
                  d_nnz_map2d[row].insert(col);
                } else {
                  o_nnz_map2d[row].insert(col);
                }
              }
            }
          }
          // 3D fields
          for (int z = 0; z < mesh->LocalNz; z++) {
            const int ind = ROUND(index(x, y, z)) - Istart;

            for (int i = 0; i < n3d; i++) {
              PetscInt row = ind + i;
              if (z == 0) {
                row += n2d;
              }

              // Depends on 2D fields
              for (int j = 0; j < n2d; j++) {
                const PetscInt col = ind0 + j;
                if (col >= Istart && col < Iend) {
                  d_nnz_map2d[row].insert(col);
                } else {
                  o_nnz_map2d[row].insert(col);
                }
              }

              // Star pattern
              for (const auto& [x_off, y_off] : xy_offsets) {
                const int xi = x + x_off;
                const int yi = y + y_off;

                if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx)
                    || (yi >= mesh->LocalNy)) {
                  continue;
                }

                int ind2 = ROUND(index(xi, yi, 0));
                if (ind2 < 0) {
                  continue; // Boundary point
                }

                if (z == 0) {
                  ind2 += n2d;
                }

                // 3D fields on this cell
                for (int j = 0; j < n3d; j++) {
                  const PetscInt col = ind2 + j;
                  if (col >= Istart && col < Iend) {
                    d_nnz_map3d[row].insert(col);
                  } else {
                    o_nnz_map3d[row].insert(col);
                  }
                }
              }
            }
          }
        }
      }

      d_nnz.reserve(nlocal);
      d_nnz.reserve(nlocal);

      for (int i = 0; i < nlocal; ++i) {
        // Assume all elements in the z direction are potentially coupled
        d_nnz.emplace_back(d_nnz_map3d[i].size() * mesh->LocalNz + d_nnz_map2d[i].size());
        o_nnz.emplace_back(o_nnz_map3d[i].size() * mesh->LocalNz + o_nnz_map2d[i].size());
      }
    }

    output_progress.write("Pre-allocating Jacobian\n");
    // Pre-allocate
    MatMPIAIJSetPreallocation(Jfd, 0, d_nnz.data(), 0, o_nnz.data());
    MatSeqAIJSetPreallocation(Jfd, 0, d_nnz.data());
    MatSetUp(Jfd);
    MatSetOption(Jfd, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

    //////////////////////////////////////////////////
    // Mark non-zero entries

    output_progress.write("Marking non-zero Jacobian entries\n");
    PetscScalar val = 1.0;
    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      for (int y = mesh->ystart; y <= mesh->yend; y++) {

        const int ind0 = ROUND(index(x, y, 0));

        // 2D fields
        for (int i = 0; i < n2d; i++) {
          const PetscInt row = ind0 + i;

          // Loop through each point in the stencil
          for (const auto& [x_off, y_off] : xy_offsets) {
            const int xi = x + x_off;
            const int yi = y + y_off;
            if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx) || (yi >= mesh->LocalNy)) {
              continue;
            }

            int ind2 = ROUND(index(xi, yi, 0));
            if (ind2 < 0) {
              continue; // A boundary point
            }

            // Depends on all variables on this cell
            for (int j = 0; j < n2d; j++) {
              PetscInt col = ind2 + j;
              ierr = MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES);
              CHKERRQ(ierr);
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
              ierr = MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES);
              CHKERRQ(ierr);
            }

            // Star pattern
            for (const auto& [x_off, y_off] : xy_offsets) {
              int xi = x + x_off;
              int yi = y + y_off;

              if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx)
                  || (yi >= mesh->LocalNy)) {
                continue;
              }
              for (int zi = 0; zi < mesh->LocalNz; ++zi) {
                int ind2 = ROUND(index(xi, yi, zi));
                if (ind2 < 0) {
                  continue; // Boundary point
                }

                if (z == 0) {
                  ind2 += n2d;
                }

                // 3D fields on this cell
                for (int j = 0; j < n3d; j++) {
                  PetscInt col = ind2 + j;
                  ierr = MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES);

                  if (ierr != 0) {
                    output.write("ERROR: {} {} : ({}, {}) -> ({}, {}) : {} -> {}\n", row,
                                 x, y, xi, yi, ind2, ind2 + n3d - 1);
                  }
                  CHKERRQ(ierr);
                }
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

    // The above will probably miss some non-zero entries at process boundaries
    // Making sure the colouring matrix is symmetric will in some/all(?)
    // of the missing non-zeros
    if (options["force_symmetric_coloring"]
            .doc("Modifies coloring matrix to force it to be symmetric")
            .withDefault<bool>(false)) {
      Mat Jfd_T;
      MatCreateTranspose(Jfd, &Jfd_T);
      MatAXPY(Jfd, 1, Jfd_T, DIFFERENT_NONZERO_PATTERN);
    }

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

    MatSetOption(Jfd, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  }
}

void PetscPreconditioner::updateColoring() {
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
