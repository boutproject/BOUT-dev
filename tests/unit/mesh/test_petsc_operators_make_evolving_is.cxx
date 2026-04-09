#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "gtest/gtest.h"

#include "bout/array.hxx"
#include "bout/bout_types.hxx"
#include "bout/output.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/petsc_operators.hxx"
#include "bout/region.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include <algorithm>
#include <memory>
#include <vector>

#include <petscis.h>

using bout::globals::mesh;

// ============================================================================
// Helper functions shared across unit tests
// ============================================================================

namespace {

/// Build a PetscCellMapping whose evolving region has exactly the cells listed
/// in cell_number with non-negative values, plus whatever boundary cells are in
/// forward_cell_number / backward_cell_number.  total_cells must equal the
/// number of distinct non-negative entries across all three fields.
std::shared_ptr<const PetscCellMapping>
makeMappingFromFields(const Field3D& cell_number, const Field3D& forward_cell_number,
                      const Field3D& backward_cell_number, int total_cells) {
  return std::make_shared<const PetscCellMapping>(cell_number, forward_cell_number,
                                                  backward_cell_number, total_cells);
}

/// Collect the global PETSc indices that mapOwnedInteriorCells visits, in order.
std::vector<PetscInt> collectOwnedInteriorIndices(const PetscCellMapping& mapping) {
  std::vector<PetscInt> indices;
  mapping.mapOwnedInteriorCells(
      [&](PetscInt row, const Ind3D& /*i*/, int /*stored*/) { indices.push_back(row); });
  return indices;
}

/// Unwrap an IS into a sorted std::vector<PetscInt>.  The IS is not destroyed.
std::vector<PetscInt> isToSortedVector(IS is) {
  PetscInt n;
  ISGetLocalSize(is, &n);
  const PetscInt* idx;
  ISGetIndices(is, &idx);
  std::vector<PetscInt> v(idx, idx + n);
  ISRestoreIndices(is, &idx);
  std::sort(v.begin(), v.end());
  return v;
}
} // Namespace

// ============================================================================
// Fixture
// ============================================================================

// FakeMeshFixture provides a mesh with:
//   LocalNx = 5  (xstart=1, xend=3, so mxg=1, 3 interior x-points)
//   LocalNy = 4  (ystart=1, yend=2, so mgy=1, 2 interior y-points)
//   LocalNz = 2
// Interior cells: x in [1,3], y in [1,2], z in [0,1] => 3*2*2 = 12 evolving cells.
using PetscEvolvingISTest = FakeMeshFixture;

// ============================================================================
// Helper: build the standard mapping used by most tests in this suite.
//
// Stored cell numbers are assigned to interior cells only (indices 0..11).
// No forward or backward boundary cells.
// ============================================================================
namespace {
std::shared_ptr<const PetscCellMapping> buildStandardMapping(Mesh* mesh) {
  Field3D cell_number{-1.0, mesh};

  int next = 0;
  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        cell_number(x, y, z) = next++;
      }
    }
  }

  const Field3D forward_cell_number{-1.0, mesh};
  const Field3D backward_cell_number{-1.0, mesh};

  return makeMappingFromFields(cell_number, forward_cell_number, backward_cell_number,
                               next);
}
} // Namespace

// ============================================================================
// IS_size_equals_n_evolving
//
// The IS produced by makeEvolvingIS() must contain exactly one index per
// evolving interior cell.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_size_equals_n_evolving) {
  auto mapping = buildStandardMapping(mesh);

  // FakeMesh: 3 x-interior * 2 y-interior * 2 z = 12 evolving cells
  const PetscInt expected_n_evolving =
      (mesh->xend - mesh->xstart + 1) * (mesh->yend - mesh->ystart + 1) * mesh->LocalNz;

  IS is = mapping->makeEvolvingIS();

  PetscInt local_size;
  ISGetLocalSize(is, &local_size);

  EXPECT_EQ(expected_n_evolving, local_size);

  ISDestroy(&is);
}

// ============================================================================
// IS_indices_match_mapOwnedInteriorCells_order
//
// The IS indices must appear in exactly the same order as the global PETSc
// row indices produced by mapOwnedInteriorCells.  Both traversals must visit
// the same cells in the same sequence.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_indices_match_mapOwnedInteriorCells_order) {
  auto mapping = buildStandardMapping(mesh);

  const std::vector<PetscInt> from_map = collectOwnedInteriorIndices(*mapping);

  IS is = mapping->makeEvolvingIS();
  PetscInt n;
  ISGetLocalSize(is, &n);
  const PetscInt* idx;
  ISGetIndices(is, &idx);
  const std::vector<PetscInt> from_is(idx, idx + n);
  ISRestoreIndices(is, &idx);
  ISDestroy(&is);

  ASSERT_EQ(from_map.size(), from_is.size());
  for (std::size_t i = 0; i < from_map.size(); ++i) {
    EXPECT_EQ(from_map[i], from_is[i]) << "Index mismatch at position " << i;
  }
}

// ============================================================================
// IS_indices_are_in_global_petsc_range
//
// On a single MPI rank the global PETSc row range is [rowStart(), rowEnd()).
// Every index in the IS must lie in that range.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_indices_are_in_global_petsc_range) {
  auto mapping = buildStandardMapping(mesh);

  IS is = mapping->makeEvolvingIS();
  const auto v = isToSortedVector(is);
  ISDestroy(&is);

  const PetscInt row_start = mapping->rowStart();
  const PetscInt row_end = mapping->rowEnd();

  for (PetscInt idx : v) {
    EXPECT_GE(idx, row_start) << "Index " << idx << " is below rowStart()";
    EXPECT_LT(idx, row_end) << "Index " << idx << " is at or above rowEnd()";
  }
}

// ============================================================================
// IS_excludes_xin_boundary_cells
//
// Cells in the inner-X boundary region have stored indices assigned by the
// mapping but must NOT appear in the evolving IS.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_excludes_xin_boundary_cells) {
  Field3D cell_number{-1.0, mesh};
  Field3D forward_cell_number{-1.0, mesh};
  Field3D backward_cell_number{-1.0, mesh};

  // Assign evolving cells
  int next = 0;
  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        cell_number(x, y, z) = next++;
      }
    }
  }
  // Assign one inner-X boundary cell a stored index
  const int xin_stored = next++;
  cell_number(0, mesh->ystart, 0) = xin_stored;

  auto mapping =
      makeMappingFromFields(cell_number, forward_cell_number, backward_cell_number, next);

  const PetscInt xin_petsc = mapping->storedToPetsc(xin_stored);

  IS is = mapping->makeEvolvingIS();
  const auto v = isToSortedVector(is);
  ISDestroy(&is);

  EXPECT_EQ(v.end(), std::find(v.begin(), v.end(), xin_petsc))
      << "Inner-X boundary cell (PETSc index " << xin_petsc
      << ") should not appear in the evolving IS";
}

// ============================================================================
// IS_excludes_xout_boundary_cells
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_excludes_xout_boundary_cells) {
  Field3D cell_number{-1.0, mesh};
  const Field3D forward_cell_number{-1.0, mesh};
  const Field3D backward_cell_number{-1.0, mesh};

  int next = 0;
  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        cell_number(x, y, z) = next++;
      }
    }
  }
  // Outer-X boundary: x == LocalNx - 1
  const int xout_stored = next++;
  cell_number(mesh->LocalNx - 1, mesh->ystart, 0) = xout_stored;

  auto mapping =
      makeMappingFromFields(cell_number, forward_cell_number, backward_cell_number, next);

  const PetscInt xout_petsc = mapping->storedToPetsc(xout_stored);

  IS is = mapping->makeEvolvingIS();
  const auto v = isToSortedVector(is);
  ISDestroy(&is);

  EXPECT_EQ(v.end(), std::find(v.begin(), v.end(), xout_petsc))
      << "Outer-X boundary cell (PETSc index " << xout_petsc
      << ") should not appear in the evolving IS";
}

// ============================================================================
// IS_excludes_forward_boundary_cells
//
// Virtual forward-parallel boundary cells (yup) live in forward_cell_number,
// not cell_number.  They must not appear in the evolving IS.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_excludes_forward_boundary_cells) {
  Field3D cell_number{-1.0, mesh};
  Field3D forward_cell_number{-1.0, mesh};
  const Field3D backward_cell_number{-1.0, mesh};

  int next = 0;
  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        cell_number(x, y, z) = next++;
      }
    }
  }
  // One virtual forward boundary cell
  const int fwd_stored = next++;
  forward_cell_number(mesh->xstart, mesh->yend, 0) = fwd_stored;

  auto mapping =
      makeMappingFromFields(cell_number, forward_cell_number, backward_cell_number, next);

  const PetscInt fwd_petsc = mapping->storedToPetsc(fwd_stored);

  IS is = mapping->makeEvolvingIS();
  const auto v = isToSortedVector(is);
  ISDestroy(&is);

  EXPECT_EQ(v.end(), std::find(v.begin(), v.end(), fwd_petsc))
      << "Forward boundary cell (PETSc index " << fwd_petsc
      << ") should not appear in the evolving IS";
}

// ============================================================================
// IS_excludes_backward_boundary_cells
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_excludes_backward_boundary_cells) {
  Field3D cell_number{-1.0, mesh};
  const Field3D forward_cell_number{-1.0, mesh};
  Field3D backward_cell_number{-1.0, mesh};

  int next = 0;
  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        cell_number(x, y, z) = next++;
      }
    }
  }
  const int bwd_stored = next++;
  backward_cell_number(mesh->xstart, mesh->ystart, 0) = bwd_stored;

  auto mapping =
      makeMappingFromFields(cell_number, forward_cell_number, backward_cell_number, next);

  const PetscInt bwd_petsc = mapping->storedToPetsc(bwd_stored);

  IS is = mapping->makeEvolvingIS();
  const auto v = isToSortedVector(is);
  ISDestroy(&is);

  EXPECT_EQ(v.end(), std::find(v.begin(), v.end(), bwd_petsc))
      << "Backward boundary cell (PETSc index " << bwd_petsc
      << ") should not appear in the evolving IS";
}

// ============================================================================
// IS_is_empty_when_no_evolving_cells
//
// If cell_number is entirely -1, there are no evolving cells and the IS must
// have size zero.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_is_empty_when_no_evolving_cells) {
  const Field3D cell_number{-1.0, mesh};
  const Field3D forward_cell_number{-1.0, mesh};
  const Field3D backward_cell_number{-1.0, mesh};

  // total_cells = 0: no cells at all
  auto mapping =
      makeMappingFromFields(cell_number, forward_cell_number, backward_cell_number, 0);

  IS is = mapping->makeEvolvingIS();
  PetscInt local_size;
  ISGetLocalSize(is, &local_size);
  ISDestroy(&is);

  EXPECT_EQ(0, local_size);
}

// ============================================================================
// IS_covers_all_evolving_stored_indices
//
// The set of stored indices reachable from the IS (via storedToPetsc inverse)
// must exactly equal the set of stored indices of interior cells as visited by
// mapOwnedInteriorCells.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_covers_all_evolving_stored_indices) {
  auto mapping = buildStandardMapping(mesh);

  // Collect stored indices from mapOwnedInteriorCells
  std::vector<int> stored_from_map;
  mapping->mapOwnedInteriorCells([&](PetscInt /*row*/, const Ind3D& /*i*/, int stored) {
    stored_from_map.push_back(stored);
  });
  std::sort(stored_from_map.begin(), stored_from_map.end());

  // Collect PETSc indices from IS, then map back to stored via localStoredIndices
  IS is = mapping->makeEvolvingIS();
  PetscInt n;
  ISGetLocalSize(is, &n);
  const PetscInt* idx;
  ISGetIndices(is, &idx);

  // localStoredIndices[i] is the stored index for PETSc row rowStart()+i
  const auto& local_stored = mapping->localStoredIndices();
  const PetscInt rstart = mapping->rowStart();
  std::vector<int> stored_from_is;
  for (PetscInt k = 0; k < n; ++k) {
    const PetscInt local_row = idx[k] - rstart;
    ASSERT_GE(local_row, 0);
    ASSERT_LT(local_row, static_cast<PetscInt>(local_stored.size()));
    stored_from_is.push_back(local_stored[local_row]);
  }
  ISRestoreIndices(is, &idx);
  ISDestroy(&is);

  std::sort(stored_from_is.begin(), stored_from_is.end());

  EXPECT_EQ(stored_from_map, stored_from_is);
}

// ============================================================================
// IS_indices_are_unique
//
// No PETSc index should appear more than once in the IS.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_indices_are_unique) {
  auto mapping = buildStandardMapping(mesh);

  IS is = mapping->makeEvolvingIS();
  const auto sorted = isToSortedVector(is);
  ISDestroy(&is);

  for (std::size_t i = 1; i < sorted.size(); ++i) {
    EXPECT_LT(sorted[i - 1], sorted[i])
        << "Duplicate index " << sorted[i] << " found at positions " << i - 1 << " and "
        << i;
  }
}

// ============================================================================
// IS_with_mixed_boundary_cells_has_correct_size
//
// When forward, backward, and X-boundary cells are all present alongside
// evolving cells, the IS size equals only the evolving cell count.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_with_mixed_boundary_cells_has_correct_size) {
  Field3D cell_number{-1.0, mesh};
  Field3D forward_cell_number{-1.0, mesh};
  Field3D backward_cell_number{-1.0, mesh};

  int next = 0;
  // Evolving interior cells
  const int n_evolving =
      (mesh->xend - mesh->xstart + 1) * (mesh->yend - mesh->ystart + 1) * mesh->LocalNz;
  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        cell_number(x, y, z) = next++;
      }
    }
  }
  // One of each boundary type
  cell_number(0, mesh->ystart, 0) = next++;                     // xin
  cell_number(mesh->LocalNx - 1, mesh->ystart, 0) = next++;     // xout
  forward_cell_number(mesh->xstart, mesh->yend, 0) = next++;    // yup
  backward_cell_number(mesh->xstart, mesh->ystart, 0) = next++; // ydown

  auto mapping =
      makeMappingFromFields(cell_number, forward_cell_number, backward_cell_number, next);

  IS is = mapping->makeEvolvingIS();
  PetscInt local_size;
  ISGetLocalSize(is, &local_size);
  ISDestroy(&is);

  EXPECT_EQ(n_evolving, local_size);
}

// ============================================================================
// IS_can_be_used_to_create_submatrix
//
// Regression guard: the IS must be acceptable to MatCreateSubMatrix without
// error.  Uses a simple identity-pattern matrix over the full cell space.
// ============================================================================
TEST_F(PetscEvolvingISTest, IS_can_be_used_to_create_submatrix) {
  auto mapping = buildStandardMapping(mesh);

  // Build a diagonal matrix over the full cell space
  const PetscInt n = mapping->globalSize();
  const PetscInt nlocal = mapping->localSize();
  Mat A;
  MatCreate(BoutComm::get(), &A);
  MatSetSizes(A, nlocal, nlocal, n, n);
  MatSetType(A, MATMPIAIJ);
  MatMPIAIJSetPreallocation(A, 1, nullptr, 0, nullptr);
  const PetscInt rstart = mapping->rowStart();
  for (PetscInt i = 0; i < nlocal; ++i) {
    const PetscInt global = rstart + i;
    const PetscScalar val = 1.0;
    MatSetValues(A, 1, &global, 1, &global, &val, INSERT_VALUES);
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  IS is = mapping->makeEvolvingIS();

  Mat sub = nullptr;
  // This must not throw or return a PETSc error
  const PetscErrorCode ierr = MatCreateSubMatrix(A, is, is, MAT_INITIAL_MATRIX, &sub);

  EXPECT_EQ(PETSC_SUCCESS, ierr);

  // The submatrix size should equal n_evolving × n_evolving
  if (sub != nullptr) {
    PetscInt sub_rows, sub_cols;
    MatGetSize(sub, &sub_rows, &sub_cols);
    const PetscInt n_evolving =
        (mesh->xend - mesh->xstart + 1) * (mesh->yend - mesh->ystart + 1) * mesh->LocalNz;
    EXPECT_EQ(n_evolving, sub_rows);
    EXPECT_EQ(n_evolving, sub_cols);
    MatDestroy(&sub);
  }

  ISDestroy(&is);
  MatDestroy(&A);
}

#endif // BOUT_HAS_PETSC
