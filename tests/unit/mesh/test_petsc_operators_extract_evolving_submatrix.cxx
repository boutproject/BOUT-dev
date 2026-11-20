#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "gtest/gtest.h"

#include "bout/array.hxx"
#include "bout/bout_types.hxx"
#include "bout/globals.hxx"
#include "bout/petsc_operators.hxx"
#include "bout/region.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <petscis.h>

using bout::globals::mesh;

// ============================================================================
// PetscCellMapping::extractEvolvingSubmatrix()
//
// Concept: given a PetscCellOperator with a known nonzero pattern, the
// submatrix returned by extractEvolvingSubmatrix() contains exactly the rows
// and columns that correspond to evolving interior cells, with values
// preserved, and no entries from boundary cells.
//
// The fixture is FakeMeshFixture (LocalNx=5, LocalNy=4, LocalNz=2, mxg=1).
// Interior cells: x in [1,3], y in [1,2], z in [0,1] => 12 evolving cells.
// ============================================================================

// ── Helpers ─────────────────────────────────────────────────────────────────

namespace {

// Build a mapping where evolving cells have stored indices [0, n_evolving),
// with optional extra boundary cells appended after.
struct MappingWithCells {
  std::shared_ptr<const PetscCellMapping> mapping;
  int n_evolving;
  // Global PETSc index of each evolving cell, in mapOwnedInteriorCells order
  std::vector<PetscInt> evolving_petsc;
  // Stored indices of boundary cells (xin, xout, yup, ydown)
  std::vector<int> boundary_stored;
};

MappingWithCells buildMappingWithBoundaries(Mesh* mesh, bool add_xin = false,
                                            bool add_xout = false, bool add_yup = false,
                                            bool add_ydown = false) {
  Field3D cell_number{-1.0, mesh};
  Field3D forward_cell_number{-1.0, mesh};
  Field3D backward_cell_number{-1.0, mesh};

  int next = 0;
  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        cell_number(x, y, z) = next++;
      }
    }
  }
  const int n_evolving = next;

  std::vector<int> boundary_stored;
  if (add_xin) {
    cell_number(0, mesh->ystart, 0) = next;
    boundary_stored.push_back(next++);
  }
  if (add_xout) {
    cell_number(mesh->LocalNx - 1, mesh->ystart, 0) = next;
    boundary_stored.push_back(next++);
  }
  if (add_yup) {
    forward_cell_number(mesh->xstart, mesh->yend, 0) = next;
    boundary_stored.push_back(next++);
  }
  if (add_ydown) {
    backward_cell_number(mesh->xstart, mesh->ystart, 0) = next;
    boundary_stored.push_back(next++);
  }

  auto mapping = std::make_shared<const PetscCellMapping>(
      cell_number, forward_cell_number, backward_cell_number, next);

  std::vector<PetscInt> evolving_petsc;
  mapping->mapOwnedInteriorCells(
      [&](PetscInt row, const Ind3D&, int) { evolving_petsc.push_back(row); });

  return {std::move(mapping), n_evolving, std::move(evolving_petsc),
          std::move(boundary_stored)};
}

// Build a PetscCellOperator from explicit (row, col, value) triples over a
// given mapping.  Rows and cols are stored indices (not PETSc global).
// The operator is assembled from a hand-built CSR representation.
PetscCellOperator
buildOperatorFromTriples(std::shared_ptr<const PetscCellMapping> mapping,
                         const std::vector<std::tuple<int, int, BoutReal>>& triples) {

  const int n = mapping->globalSize();
  // Build dense CSR: rows array size n+1
  std::vector<int> row_counts(n, 0);
  for (const auto& [r, c, v] : triples) {
    row_counts[r]++;
  }
  // rows (CSR row pointer)
  Array<int> rows(n + 1);
  rows[0] = 0;
  for (int i = 0; i < n; ++i) {
    rows[i + 1] = rows[i] + row_counts[i];
  }
  // cols and weights
  Array<int> cols(static_cast<int>(triples.size()));
  Array<BoutReal> weights(static_cast<int>(triples.size()));
  std::vector<int> fill(n, 0);
  for (const auto& [r, c, v] : triples) {
    const int pos = rows[r] + fill[r]++;
    cols[pos] = c;
    weights[pos] = v;
  }

  return PetscCellOperator(mapping, mapping, rows, cols, weights);
}

// Return a flat vector of (global_row, global_col, value) for every nonzero
// in a matrix, using MatGetRow.  Rows are iterated over the local ownership
// range only.
struct Nonzero {
  PetscInt row, col;
  PetscScalar val;
};

std::vector<Nonzero> matNonzeros(Mat A) {
  PetscInt rstart, rend;
  MatGetOwnershipRange(A, &rstart, &rend);
  std::vector<Nonzero> nzs;
  for (PetscInt row = rstart; row < rend; ++row) {
    PetscInt ncols;
    const PetscInt* col_idx;
    const PetscScalar* vals;
    MatGetRow(A, row, &ncols, &col_idx, &vals);
    for (PetscInt k = 0; k < ncols; ++k) {
      nzs.push_back({row, col_idx[k], vals[k]});
    }
    MatRestoreRow(A, row, &ncols, &col_idx, &vals);
  }
  return nzs;
}

} // namespace

// ── Fixture ─────────────────────────────────────────────────────────────────

using PetscExtractEvolvingTest = FakeMeshFixture;

// ============================================================================
// submatrix_global_size_equals_n_evolving_squared
//
// MatGetSize on the result must report n_evolving × n_evolving.
// ============================================================================
TEST_F(PetscExtractEvolvingTest, submatrix_global_size_equals_n_evolving_squared) {
  auto [mapping, n_evolving, evolving_petsc, boundary_stored] =
      buildMappingWithBoundaries(mesh, /*xin=*/true, /*xout=*/true);

  // Diagonal operator over the full cell space (safe: every cell maps to itself)
  const auto& local = mapping->localStoredIndices();
  std::vector<std::tuple<int, int, BoutReal>> triples;
  for (int s : local) {
    if (s >= 0) {
      triples.emplace_back(s, s, 1.0);
    }
  }
  const auto op = buildOperatorFromTriples(mapping, triples);

  Mat sub = mapping->extractEvolvingSubmatrix(op);

  PetscInt sub_rows, sub_cols;
  MatGetSize(sub, &sub_rows, &sub_cols);

  EXPECT_EQ(n_evolving, sub_rows);
  EXPECT_EQ(n_evolving, sub_cols);

  MatDestroy(&sub);
}

// ============================================================================
// submatrix_local_size_equals_n_evolving_locally
//
// MatGetLocalSize must report the local evolving count on each rank.
// ============================================================================
TEST_F(PetscExtractEvolvingTest, submatrix_local_size_equals_n_evolving_locally) {
  auto [mapping, n_evolving, evolving_petsc, boundary_stored] =
      buildMappingWithBoundaries(mesh, true, true, true, true);

  std::vector<std::tuple<int, int, BoutReal>> triples;
  for (int i = 0; i < n_evolving; ++i) {
    triples.emplace_back(i, i, 1.0);
  }
  const auto op = buildOperatorFromTriples(mapping, triples);

  Mat sub = mapping->extractEvolvingSubmatrix(op);

  PetscInt local_rows, local_cols;
  MatGetLocalSize(sub, &local_rows, &local_cols);

  // On a single rank local == global
  EXPECT_EQ(n_evolving, local_rows);
  EXPECT_EQ(n_evolving, local_cols);

  MatDestroy(&sub);
}

// ============================================================================
// submatrix_values_preserved_for_evolving_block
//
// For entries (r, c) where both r and c are evolving cells, the value in the
// submatrix must equal the value in the original operator.
//
// Strategy: build an operator with known values on evolving-to-evolving
// entries, extract the submatrix, and check each value via MatGetValues.
// ============================================================================
TEST_F(PetscExtractEvolvingTest, submatrix_values_preserved_for_evolving_block) {
  auto [mapping, n_evolving, evolving_petsc, boundary_stored] =
      buildMappingWithBoundaries(mesh, /*xin=*/true, /*xout=*/true,
                                 /*yup=*/true, /*ydown=*/true);

  // Insert entries with distinct values between pairs of evolving cells.
  // Use stored indices directly (0..n_evolving-1).
  std::vector<std::tuple<int, int, BoutReal>> triples;
  for (int i = 0; i < n_evolving; ++i) {
    // Diagonal entry
    triples.emplace_back(i, i, static_cast<BoutReal>(10 + i));
    // One off-diagonal entry per row (wrap around)
    const int j = (i + 1) % n_evolving;
    triples.emplace_back(i, j, static_cast<BoutReal>(100 + i));
  }
  const auto op = buildOperatorFromTriples(mapping, triples);

  Mat sub = mapping->extractEvolvingSubmatrix(op);
  const auto nzs = matNonzeros(sub);

  // Build expected map: submatrix uses 0-based indices within the evolving
  // block.  The IS preserves mapOwnedInteriorCells order, which matches
  // stored index order (0..n_evolving-1) on a single rank.
  // Sub row/col i corresponds to evolving stored index i.
  PetscInt sub_rstart, sub_rend;
  MatGetOwnershipRange(sub, &sub_rstart, &sub_rend);

  for (const auto& nz : nzs) {
    // nz.row and nz.col are global indices in the submatrix space.
    // On a single rank these equal the stored indices of the evolving cells.
    const int stored_row = static_cast<int>(nz.row - sub_rstart);
    const int stored_col = static_cast<int>(nz.col);

    // Find the expected value from our triples
    BoutReal expected = 0.0;
    for (const auto& [r, c, v] : triples) {
      if (r == stored_row && c == stored_col) {
        expected = v;
        break;
      }
    }
    EXPECT_DOUBLE_EQ(expected, nz.val)
        << "Wrong value at submatrix (" << nz.row << ", " << nz.col << ")";
  }

  MatDestroy(&sub);
}

// ============================================================================
// submatrix_nonzero_count_matches_evolving_entries
//
// The total nonzero count in the submatrix must equal the number of entries
// in the original operator whose row AND column are both evolving cells.
// Entries coupling to boundary cells must be absent.
// ============================================================================
TEST_F(PetscExtractEvolvingTest, submatrix_nonzero_count_matches_evolving_entries) {
  auto [mapping, n_evolving, evolving_petsc, boundary_stored] =
      buildMappingWithBoundaries(mesh, /*xin=*/true, /*xout=*/false,
                                 /*yup=*/true, /*ydown=*/false);

  std::vector<std::tuple<int, int, BoutReal>> triples;
  // Entries entirely within the evolving block
  int expected_in_sub = 0;
  for (int i = 0; i < n_evolving; ++i) {
    triples.emplace_back(i, i, 1.0);
    ++expected_in_sub;
  }
  // Entries from an evolving row into a boundary column — must be excluded
  for (int bnd : boundary_stored) {
    triples.emplace_back(0, bnd, 99.0); // row 0 (evolving) -> boundary col
  }
  // Entries from a boundary row into an evolving column — must be excluded
  for (int bnd : boundary_stored) {
    triples.emplace_back(bnd, 0, 88.0); // boundary row -> col 0 (evolving)
  }

  const auto op = buildOperatorFromTriples(mapping, triples);
  Mat sub = mapping->extractEvolvingSubmatrix(op);

  PetscInt nnz;
  MatInfo info;
  MatGetInfo(sub, MAT_GLOBAL_SUM, &info);
  nnz = static_cast<PetscInt>(info.nz_used);

  EXPECT_EQ(expected_in_sub, nnz);

  MatDestroy(&sub);
}

// ============================================================================
// submatrix_excludes_boundary_columns
//
// No column index in the submatrix should correspond to a boundary cell's
// PETSc index in the original cell space.
// ============================================================================
TEST_F(PetscExtractEvolvingTest, submatrix_excludes_boundary_columns) {
  auto [mapping, n_evolving, evolving_petsc, boundary_stored] =
      buildMappingWithBoundaries(mesh, /*xin=*/true, /*xout=*/true,
                                 /*yup=*/true, /*ydown=*/true);

  std::vector<std::tuple<int, int, BoutReal>> triples;
  for (int i = 0; i < n_evolving; ++i) {
    triples.emplace_back(i, i, 1.0);
  }
  // Add evolving-row → boundary-column coupling
  for (int bnd : boundary_stored) {
    triples.emplace_back(0, bnd, 7.0);
  }

  const auto op = buildOperatorFromTriples(mapping, triples);

  // Collect the PETSc indices of the boundary cells so we can check against them.
  // The submatrix uses a re-indexed column space [0, n_evolving), so boundary
  // columns from the full space simply will not appear.
  Mat sub = mapping->extractEvolvingSubmatrix(op);

  PetscInt sub_rows, sub_cols;
  MatGetSize(sub, &sub_rows, &sub_cols);

  // The boundary-coupling entries were given the sentinel value 7.0.
  // If any survive into the submatrix, the extraction has failed to exclude them.
  const auto nzs = matNonzeros(sub);
  for (const auto& nz : nzs) {
    EXPECT_NE(7.0, nz.val)
        << "Boundary-column entry (sentinel value 7.0) survived into the "
           "submatrix at ("
        << nz.row << ", " << nz.col << ")";
  }

  MatDestroy(&sub);
}

// ============================================================================
// submatrix_excludes_boundary_rows
//
// The submatrix must have no rows for boundary cells.  Verified by checking
// that MatGetLocalSize equals n_evolving (not n_evolving + n_boundary).
// Also checked by confirming no row index outside [sub_rstart, sub_rend)
// appears in the nonzeros.
// ============================================================================
TEST_F(PetscExtractEvolvingTest, submatrix_excludes_boundary_rows) {
  auto [mapping, n_evolving, evolving_petsc, boundary_stored] =
      buildMappingWithBoundaries(mesh, /*xin=*/true, /*xout=*/false,
                                 /*yup=*/false, /*ydown=*/true);

  std::vector<std::tuple<int, int, BoutReal>> triples;
  for (int i = 0; i < n_evolving; ++i) {
    triples.emplace_back(i, i, 1.0);
  }
  // boundary-row → evolving-column entries
  for (int bnd : boundary_stored) {
    triples.emplace_back(bnd, 0, 5.0);
  }

  const auto op = buildOperatorFromTriples(mapping, triples);
  Mat sub = mapping->extractEvolvingSubmatrix(op);

  PetscInt local_rows, local_cols;
  MatGetLocalSize(sub, &local_rows, &local_cols);
  EXPECT_EQ(n_evolving, local_rows);

  PetscInt sub_rstart, sub_rend;
  MatGetOwnershipRange(sub, &sub_rstart, &sub_rend);
  const auto nzs = matNonzeros(sub);
  for (const auto& nz : nzs) {
    EXPECT_GE(nz.row, sub_rstart);
    EXPECT_LT(nz.row, sub_rend);
  }

  MatDestroy(&sub);
}

// ============================================================================
// zero_operator_gives_zero_submatrix
//
// A zero operator (no nonzeros) must produce an assembled zero submatrix of
// the correct size.  No crash, no uninitialized memory.
// ============================================================================
TEST_F(PetscExtractEvolvingTest, zero_operator_gives_zero_submatrix) {
  auto [mapping, n_evolving, evolving_petsc, boundary_stored] =
      buildMappingWithBoundaries(mesh);

  // Empty triples: zero operator
  const auto op = buildOperatorFromTriples(mapping, {});

  Mat sub = mapping->extractEvolvingSubmatrix(op);

  PetscInt sub_rows, sub_cols;
  MatGetSize(sub, &sub_rows, &sub_cols);
  EXPECT_EQ(n_evolving, sub_rows);
  EXPECT_EQ(n_evolving, sub_cols);

  MatInfo info;
  MatGetInfo(sub, MAT_GLOBAL_SUM, &info);
  EXPECT_EQ(0, static_cast<PetscInt>(info.nz_used));

  MatDestroy(&sub);
}

// ============================================================================
// submatrix_of_identity_is_identity
//
// An identity operator over the evolving cells only (no boundary entries)
// should produce an identity submatrix: diagonal entries of 1.0, nothing else.
// ============================================================================
TEST_F(PetscExtractEvolvingTest, submatrix_of_identity_is_identity) {
  auto [mapping, n_evolving, evolving_petsc, boundary_stored] =
      buildMappingWithBoundaries(mesh);

  std::vector<std::tuple<int, int, BoutReal>> triples;
  for (int i = 0; i < n_evolving; ++i) {
    triples.emplace_back(i, i, 1.0);
  }
  const auto op = buildOperatorFromTriples(mapping, triples);

  Mat sub = mapping->extractEvolvingSubmatrix(op);
  const auto nzs = matNonzeros(sub);

  PetscInt sub_rstart, sub_rend;
  MatGetOwnershipRange(sub, &sub_rstart, &sub_rend);

  ASSERT_EQ(n_evolving, static_cast<int>(nzs.size()));
  for (const auto& nz : nzs) {
    EXPECT_EQ(nz.row, nz.col) << "Off-diagonal entry in identity submatrix";
    EXPECT_DOUBLE_EQ(1.0, nz.val);
  }

  MatDestroy(&sub);
}

// ============================================================================
// submatrix_is_independent_of_original_operator
//
// Modifying the original operator's underlying matrix after extraction must
// not affect the already-returned submatrix.  MatCreateSubMatrix with
// MAT_INITIAL_MATRIX produces an independent copy.
// ============================================================================
TEST_F(PetscExtractEvolvingTest, submatrix_is_independent_of_original_operator) {
  auto [mapping, n_evolving, evolving_petsc, boundary_stored] =
      buildMappingWithBoundaries(mesh);

  std::vector<std::tuple<int, int, BoutReal>> triples;
  for (int i = 0; i < n_evolving; ++i) {
    triples.emplace_back(i, i, 2.0);
  }
  auto op = buildOperatorFromTriples(mapping, triples);

  Mat sub = mapping->extractEvolvingSubmatrix(op);

  // Scale the original operator's matrix in place
  MatScale(op.raw(), 0.0);

  // The submatrix must still hold the original values
  const auto nzs = matNonzeros(sub);
  for (const auto& nz : nzs) {
    EXPECT_DOUBLE_EQ(2.0, nz.val)
        << "Submatrix was affected by modification of the source operator";
  }

  MatDestroy(&sub);
}

#endif // BOUT_HAS_PETSC
