#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "bout/petsc_preconditioner.hxx"

#include "bout/assert.hxx"
#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"
#include "bout/utils.hxx"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {
class ColoringStencil {
private:
  static bool isInSquare(int i, int j, int n_square) {
    return std::abs(i) <= n_square && std::abs(j) <= n_square;
  }
  static bool isInCross(int i, int j, int n_cross) {
    if (i == 0) {
      return std::abs(j) <= n_cross;
    }
    if (j == 0) {
      return std::abs(i) <= n_cross;
    }
    return false;
  }
  static bool isInTaxi(int i, int j, int n_taxi) {
    return std::abs(i) + std::abs(j) <= n_taxi;
  }

public:
  static auto getOffsets(int n_square, int n_taxi, int n_cross) {
    ASSERT2(n_square >= 0 && n_cross >= 0 && n_taxi >= 0
            && n_square + n_cross + n_taxi > 0);
    auto inside = [&](int i, int j) {
      return isInSquare(i, j, n_square) || isInTaxi(i, j, n_taxi)
             || isInCross(i, j, n_cross);
    };
    std::vector<std::pair<int, int>> xy_offsets;
    auto loop_bound = std::max({n_square, n_taxi, n_cross});
    for (int i = -loop_bound; i <= loop_bound; ++i) {
      for (int j = -loop_bound; j <= loop_bound; ++j) {
        if (inside(i, j)) {
          xy_offsets.emplace_back(i, j);
        }
      }
    }
    return xy_offsets;
  }
};
} // namespace

void PetscPreconditioner::reset() {
  if (fdcoloring != nullptr) {
    MatFDColoringDestroy(&fdcoloring);
    fdcoloring = nullptr;
  }
  if (Jfd != nullptr) {
    MatDestroy(&Jfd);
    Jfd = nullptr;
  }
}

PetscErrorCode PetscPreconditioner::saveMatrix(Mat matrix, const std::string& filename,
                                               bout::PetscMatrixExportFormat format) {
  if (matrix == nullptr) {
    throw BoutException("Cannot save Jacobian matrix: matrix has not been created yet");
  }

  PetscViewer viewer{nullptr};
  if (format == bout::PetscMatrixExportFormat::binary) {
    PetscCall(PetscViewerBinaryOpen(BoutComm::get(), filename.c_str(), FILE_MODE_WRITE,
                                    &viewer));
  } else {
    PetscCall(PetscViewerASCIIOpen(BoutComm::get(), filename.c_str(), &viewer));
  }

  PetscCall(MatView(matrix, viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode
PetscPreconditioner::saveMatrix(const std::string& filename,
                                bout::PetscMatrixExportFormat format) const {
  return saveMatrix(Jfd, filename, format);
}

PetscErrorCode PetscPreconditioner::createJacobianPattern(Field3D& index,
                                                          Options& options,
                                                          PetscInt nlocal, int n2d,
                                                          int n3d, MPI_Comm comm) {
  reset();

  // Use global mesh for now
  Mesh* mesh = bout::globals::mesh;

  //////////////////////////////////////////////////
  // Pre-allocate PETSc storage

  output_progress.write("Setting Jacobian matrix sizes\n");

  // Set size of Matrix on each processor to nlocal x nlocal
  PetscCall(MatCreate(comm, &Jfd));
  PetscCall(MatSetOption(Jfd, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
  PetscCall(MatSetSizes(Jfd, nlocal, nlocal, PETSC_DETERMINE, PETSC_DETERMINE));
  PetscCall(MatSetFromOptions(Jfd));

  // Determine which row/columns of the matrix are locally owned
  int Istart = 0;
  int Iend = 0;
  PetscCall(MatGetOwnershipRange(Jfd, &Istart, &Iend));

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

  auto n_square =
      options["stencil:square"].doc("Extent of stencil (square)").withDefault<int>(0);
  auto n_cross =
      options["stencil:cross"].doc("Extent of stencil (cross)").withDefault<int>(0);
  // Default n_taxi=2 if nothing else is set
  auto n_taxi = options["stencil:taxi"]
                    .doc("Extent of stencil (taxi-cab norm)")
                    .withDefault<int>((n_square == 0 && n_cross == 0) ? 2 : 0);

  const auto xy_offsets = ColoringStencil::getOffsets(n_square, n_taxi, n_cross);
  {
    // Count the *unique* non-zeros for preallocation
    std::vector<std::set<int>> d_nnz_map2d(nlocal);
    std::vector<std::set<int>> o_nnz_map2d(nlocal);
    std::vector<std::set<int>> d_nnz_map3d(nlocal);
    std::vector<std::set<int>> o_nnz_map3d(nlocal);

    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      for (int y = mesh->ystart; y <= mesh->yend; y++) {
        const int ind0 = ROUND(index(x, y, 0)) - Istart;

        // 2D fields
        for (int i = 0; i < n2d; i++) {
          const PetscInt row = ind0 + i;
          for (const auto& [x_off, y_off] : xy_offsets) {
            const int xi = x + x_off;
            const int yi = y + y_off;
            if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx) || (yi >= mesh->LocalNy)) {
              continue;
            }

            const int ind2 = ROUND(index(xi, yi, 0));
            if (ind2 < 0) {
              continue; // A boundary point
            }

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
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
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
    o_nnz.reserve(nlocal);
    for (int i = 0; i < nlocal; ++i) {
      // Assume all elements in the z direction are potentially coupled
      d_nnz.emplace_back((d_nnz_map3d[i].size() * mesh->LocalNz) + d_nnz_map2d[i].size());
      o_nnz.emplace_back((o_nnz_map3d[i].size() * mesh->LocalNz) + o_nnz_map2d[i].size());
    }
  }

  output_progress.write("Pre-allocating Jacobian\n");
  PetscCall(MatMPIAIJSetPreallocation(Jfd, 0, d_nnz.data(), 0, o_nnz.data()));
  PetscCall(MatSeqAIJSetPreallocation(Jfd, 0, d_nnz.data()));
  PetscCall(MatSetUp(Jfd));
  PetscCall(MatSetOption(Jfd, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));

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

        for (const auto& [x_off, y_off] : xy_offsets) {
          const int xi = x + x_off;
          const int yi = y + y_off;
          if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx) || (yi >= mesh->LocalNy)) {
            continue;
          }

          const int ind2 = ROUND(index(xi, yi, 0));
          if (ind2 < 0) {
            continue; // A boundary point
          }

          for (int j = 0; j < n2d; j++) {
            const PetscInt col = ind2 + j;
            PetscCall(MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES));
          }
        }
      }

      // 3D fields
      for (int z = mesh->zstart; z <= mesh->zend; z++) {
        const int ind = ROUND(index(x, y, z));

        for (int i = 0; i < n3d; i++) {
          PetscInt row = ind + i;
          if (z == 0) {
            row += n2d;
          }

          // Depends on 2D fields
          for (int j = 0; j < n2d; j++) {
            const PetscInt col = ind0 + j;
            PetscCall(MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES));
          }

          for (const auto& [x_off, y_off] : xy_offsets) {
            const int xi = x + x_off;
            const int yi = y + y_off;

            if ((xi < 0) || (yi < 0) || (xi >= mesh->LocalNx) || (yi >= mesh->LocalNy)) {
              continue;
            }

            for (int zi = mesh->zstart; zi <= mesh->zend; ++zi) {
              int ind2 = ROUND(index(xi, yi, zi));
              if (ind2 < 0) {
                continue; // Boundary point
              }

              if (z == 0) {
                ind2 += n2d;
              }

              for (int j = 0; j < n3d; j++) {
                const PetscInt col = ind2 + j;
                PetscCall(MatSetValues(Jfd, 1, &row, 1, &col, &val, INSERT_VALUES));
              }
            }
          }
        }
      }
    }
  }

  output_progress.write("Assembling Jacobian matrix\n");
  PetscCall(MatAssemblyBegin(Jfd, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Jfd, MAT_FINAL_ASSEMBLY));

  {
    PetscBool symmetric = PETSC_FALSE;
    PetscCall(MatIsSymmetric(Jfd, 1e-5, &symmetric));
    if (!static_cast<bool>(symmetric)) {
      output_warn.write("Jacobian pattern is not symmetric\n");
    }
  }

  if (options["force_symmetric_coloring"]
          .doc("Modifies coloring matrix to force it to be symmetric")
          .withDefault<bool>(false)) {
    Mat Jfd_T{nullptr};
    PetscCall(MatCreateTranspose(Jfd, &Jfd_T));
    PetscCall(MatAXPY(Jfd, 1, Jfd_T, DIFFERENT_NONZERO_PATTERN));
    PetscCall(MatDestroy(&Jfd_T));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif // BOUT_HAS_PETSC
