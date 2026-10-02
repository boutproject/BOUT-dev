// ============================================================================
// addOperatorSparsity(Mat Jfd, Mat sub, int out_var, int in_var)
//
// Concept: the nonzero pattern of `sub` is scattered into the variable block
// (out_var, in_var) of the pre-allocated Jacobian `Jfd`, with the correct
// index stride.
//
// `sub` is a square matrix of size n_evolving x n_evolving.
// `Jfd` is a square matrix of size (n_evolving * nvars) x (n_evolving * nvars).
// nvars is inferred internally as Jfd_global / sub_global.
//
// A nonzero at (r, c) in sub produces an entry at
//   (r * nvars + out_var,  c * nvars + in_var)
// in Jfd.
//
// ============================================================================

#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include "bout/petsc_jacobian.hxx"

#include <petscmat.h>

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

// ── Helpers ──────────────────────────────────────────────────────────────────

namespace {

// Build and assemble a square MPIAIJ matrix of global size n x n,
// with the given (global_row, global_col) nonzero positions (value = 1.0).
// nlocal is the number of locally owned rows on this rank.
// rstart is the first globally owned row.
Mat buildPatternMat(PetscInt n, PetscInt nlocal, PetscInt rstart,
                    const std::vector<std::pair<PetscInt, PetscInt>>& entries) {
  Mat A;
  MatCreate(BoutComm::get(), &A);
  MatSetSizes(A, nlocal, nlocal, n, n);
  MatSetType(A, MATMPIAIJ);
  // Generous preallocation: allow up to n nonzeros per row
  MatMPIAIJSetPreallocation(A, n, nullptr, n, nullptr);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  const PetscScalar one = 1.0;
  for (const auto& [r, c] : entries) {
    MatSetValues(A, 1, &r, 1, &c, &one, INSERT_VALUES);
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  return A;
}

// Collect all (global_row, global_col) nonzero positions in a matrix,
// iterating only over locally owned rows.
std::vector<std::pair<PetscInt, PetscInt>> matNonzeroPositions(Mat A) {
  PetscInt rstart, rend;
  MatGetOwnershipRange(A, &rstart, &rend);
  std::vector<std::pair<PetscInt, PetscInt>> result;
  for (PetscInt row = rstart; row < rend; ++row) {
    PetscInt ncols;
    const PetscInt* cols;
    const PetscScalar* vals;
    MatGetRow(A, row, &ncols, &cols, &vals);
    for (PetscInt k = 0; k < ncols; ++k) {
      result.emplace_back(row, cols[k]);
    }
    MatRestoreRow(A, row, &ncols, &cols, &vals);
  }
  return result;
}

// Build a Jfd of size (n_sub * nvars) x (n_sub * nvars), with nlocal_sub
// rows per variable per rank, and rstart_sub the sub-matrix row offset.
// Returns the matrix plus the Jfd local size and row start.
struct JfdInfo {
  Mat jfd;
  PetscInt nlocal;
  PetscInt rstart;
};

JfdInfo buildEmptyJfd(PetscInt n_sub, PetscInt nlocal_sub, PetscInt rstart_sub,
                      int nvars) {
  const PetscInt n = n_sub * nvars;
  const PetscInt nlocal = nlocal_sub * nvars;

  Mat jfd;
  MatCreate(BoutComm::get(), &jfd);
  MatSetSizes(jfd, nlocal, nlocal, n, n);
  MatSetType(jfd, MATMPIAIJ);
  MatMPIAIJSetPreallocation(jfd, n, nullptr, n, nullptr);
  MatSetOption(jfd, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  PetscInt jfd_rstart, jfd_rend;
  MatAssemblyBegin(jfd, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(jfd, MAT_FINAL_ASSEMBLY);
  MatGetOwnershipRange(jfd, &jfd_rstart, &jfd_rend);
  return {jfd, nlocal, jfd_rstart};
}

// Call addOperatorSparsity and assemble Jfd.
void applyAndAssemble(Mat jfd, Mat sub, int out_var, int in_var) {
  addOperatorSparsity(jfd, sub, out_var, in_var);
  MatAssemblyBegin(jfd, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(jfd, MAT_FINAL_ASSEMBLY);
}

} // namespace

// ── Fixture ──────────────────────────────────────────────────────────────────

using PetscAddSparsityTest = FakeMeshFixture;

// ============================================================================
// single_variable_diagonal_sub_produces_diagonal_jfd
//
// With nvars=1 a diagonal sub must produce a diagonal Jfd.
// The simplest case: no variable interleaving.
// ============================================================================
TEST_F(PetscAddSparsityTest, single_variable_diagonal_sub_produces_diagonal_jfd) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub; // single rank
  const PetscInt rstart_sub = 0;
  const int nvars = 1;

  // Diagonal sub: (0,0), (1,1), (2,2)
  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {{0, 0}, {1, 1}, {2, 2}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, rstart_sub, sub_entries);

  auto [jfd, jfd_nlocal, jfd_rstart] =
      buildEmptyJfd(n_sub, nlocal_sub, rstart_sub, nvars);

  applyAndAssemble(jfd, sub, /*out_var=*/0, /*in_var=*/0);

  const auto nzs = matNonzeroPositions(jfd);
  ASSERT_EQ(3U, nzs.size());
  for (const auto& [r, c] : nzs) {
    EXPECT_EQ(r, c) << "Expected diagonal entry, got (" << r << "," << c << ")";
  }

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// correct_block_is_filled_diagonal_coupling
//
// With nvars=3, out_var=in_var=1 (self-coupling of variable 1),
// a diagonal sub must populate only the (1,1) variable block of Jfd.
// ============================================================================
TEST_F(PetscAddSparsityTest, correct_block_is_filled_diagonal_coupling) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 3;
  const int out_var = 1;
  const int in_var = 1;

  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {{0, 0}, {1, 1}, {2, 2}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  applyAndAssemble(jfd, sub, out_var, in_var);

  const auto nzs = matNonzeroPositions(jfd);
  ASSERT_EQ(3U, nzs.size());

  for (const auto& [r, c] : nzs) {
    // Row must be in the out_var=1 slot: r % nvars == 1
    EXPECT_EQ(out_var, r % nvars)
        << "Entry (" << r << "," << c << ") is in wrong row variable";
    // Col must be in the in_var=1 slot: c % nvars == 1
    EXPECT_EQ(in_var, c % nvars)
        << "Entry (" << r << "," << c << ") is in wrong col variable";
    // Cell index must match: r / nvars == c / nvars for diagonal sub
    EXPECT_EQ(r / nvars, c / nvars)
        << "Entry (" << r << "," << c << ") cell indices don't match";
  }

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// correct_block_is_filled_off_diagonal_coupling
//
// With nvars=3, out_var=0, in_var=2 (cross-variable coupling),
// all Jfd entries must be in the (0,2) variable block.
// ============================================================================
TEST_F(PetscAddSparsityTest, correct_block_is_filled_off_diagonal_coupling) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 3;
  const int out_var = 0;
  const int in_var = 2;

  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {{0, 0}, {1, 1}, {2, 2}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  applyAndAssemble(jfd, sub, out_var, in_var);

  const auto nzs = matNonzeroPositions(jfd);
  ASSERT_EQ(3U, nzs.size());

  for (const auto& [r, c] : nzs) {
    EXPECT_EQ(out_var, r % nvars)
        << "Entry (" << r << "," << c << ") is in wrong row variable";
    EXPECT_EQ(in_var, c % nvars)
        << "Entry (" << r << "," << c << ") is in wrong col variable";
  }

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// stride_is_correct_for_nvars_2
//
// With nvars=2 and a sub containing entries (0,0) and (0,1), the Jfd
// entry positions must follow the stride formula:
//   Jfd row = sub_row * nvars + out_var
//   Jfd col = sub_col * nvars + in_var
// ============================================================================
TEST_F(PetscAddSparsityTest, stride_is_correct_for_nvars_2) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 2;
  const int out_var = 0;
  const int in_var = 1;

  // Sub entries: (0,0), (1,2), (2,1)  — arbitrary non-diagonal pattern
  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {{0, 0}, {1, 2}, {2, 1}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  applyAndAssemble(jfd, sub, out_var, in_var);

  // Build expected Jfd positions using the stride formula
  std::vector<std::pair<PetscInt, PetscInt>> expected;
  for (const auto& [sr, sc] : sub_entries) {
    expected.emplace_back((sr * nvars) + out_var, (sc * nvars) + in_var);
  }
  std::sort(expected.begin(), expected.end());

  auto nzs = matNonzeroPositions(jfd);
  std::sort(nzs.begin(), nzs.end());

  ASSERT_EQ(expected.size(), nzs.size());
  for (std::size_t k = 0; k < expected.size(); ++k) {
    EXPECT_EQ(expected[k].first, nzs[k].first) << "Row mismatch at entry " << k;
    EXPECT_EQ(expected[k].second, nzs[k].second) << "Col mismatch at entry " << k;
  }

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// no_other_blocks_written
//
// With nvars=3 and only (out_var=1, in_var=1) populated, the other eight
// variable blocks must remain empty.
// ============================================================================
TEST_F(PetscAddSparsityTest, no_other_blocks_written) {
  const PetscInt n_sub = 4;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 3;
  const int out_var = 1;
  const int in_var = 1;

  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {
      {0, 0}, {1, 1}, {2, 2}, {3, 3}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  applyAndAssemble(jfd, sub, out_var, in_var);

  const auto nzs = matNonzeroPositions(jfd);
  ASSERT_EQ(sub_entries.size(), nzs.size()) << "Wrong number of entries in Jfd";
  for (const auto& [r, c] : nzs) {
    EXPECT_EQ(out_var, r % nvars) << "Unexpected entry in row-variable " << r % nvars
                                  << " at (" << r << "," << c << ")";
    EXPECT_EQ(in_var, c % nvars) << "Unexpected entry in col-variable " << c % nvars
                                 << " at (" << r << "," << c << ")";
  }

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// multiple_calls_union_patterns
//
// Two calls with different (out_var, in_var) pairs must union their patterns:
// the total nonzero count equals the sum of each operator's count, and the
// entries from each call are in the correct block.
// ============================================================================
TEST_F(PetscAddSparsityTest, multiple_calls_union_patterns) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 2;

  // Sub A: diagonal
  const std::vector<std::pair<PetscInt, PetscInt>> entries_a = {{0, 0}, {1, 1}, {2, 2}};
  Mat sub_a = buildPatternMat(n_sub, nlocal_sub, 0, entries_a);

  // Sub B: one off-diagonal entry per row
  const std::vector<std::pair<PetscInt, PetscInt>> entries_b = {{0, 1}, {1, 2}, {2, 0}};
  Mat sub_b = buildPatternMat(n_sub, nlocal_sub, 0, entries_b);

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  // Insert both without assembling between calls
  addOperatorSparsity(jfd, sub_a, /*out_var=*/0, /*in_var=*/0);
  addOperatorSparsity(jfd, sub_b, /*out_var=*/1, /*in_var=*/1);
  MatAssemblyBegin(jfd, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(jfd, MAT_FINAL_ASSEMBLY);

  const auto nzs = matNonzeroPositions(jfd);

  // Total count: 3 from sub_a + 3 from sub_b = 6
  EXPECT_EQ(6U, nzs.size());

  // Check each entry is in the correct variable block
  for (const auto& [r, c] : nzs) {
    const int rv = static_cast<int>(r % nvars);
    const int cv = static_cast<int>(c % nvars);
    const bool in_a_block = (rv == 0 && cv == 0);
    const bool in_b_block = (rv == 1 && cv == 1);
    EXPECT_TRUE(in_a_block || in_b_block)
        << "Entry (" << r << "," << c << ") is in unexpected variable block (" << rv
        << "," << cv << ")";
  }

  MatDestroy(&sub_a);
  MatDestroy(&sub_b);
  MatDestroy(&jfd);
}

// ============================================================================
// empty_sub_adds_no_entries
//
// A zero-pattern sub must not modify Jfd.
// ============================================================================
TEST_F(PetscAddSparsityTest, empty_sub_adds_no_entries) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 2;

  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, {});

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  applyAndAssemble(jfd, sub, 0, 0);

  const auto nzs = matNonzeroPositions(jfd);
  EXPECT_EQ(0U, nzs.size());

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// nonzero_count_equals_sub_nonzero_count
//
// For a single call the total Jfd nonzero count must equal the sub nonzero
// count exactly — no extras, no missing.
// ============================================================================
TEST_F(PetscAddSparsityTest, nonzero_count_equals_sub_nonzero_count) {
  const PetscInt n_sub = 4;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 3;

  // Arbitrary non-trivial pattern
  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {
      {0, 0}, {0, 1}, {1, 1}, {1, 3}, {2, 2}, {3, 0}, {3, 3}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  applyAndAssemble(jfd, sub, 0, 0);

  MatInfo info;
  MatGetInfo(jfd, MAT_GLOBAL_SUM, &info);
  EXPECT_EQ(static_cast<PetscInt>(sub_entries.size()),
            static_cast<PetscInt>(info.nz_used));

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// first_and_last_variable_blocks_correct
//
// Explicitly check out_var=0,in_var=0 and out_var=nvars-1,in_var=nvars-1
// to confirm the stride formula is correct at the boundaries.
// ============================================================================
TEST_F(PetscAddSparsityTest, first_and_last_variable_blocks_correct) {
  const PetscInt n_sub = 2;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 4;

  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {
      {0, 0}, {0, 1}, {1, 0}, {1, 1}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  // Test first block (0,0)
  {
    auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);
    applyAndAssemble(jfd, sub, 0, 0);

    const auto nzs = matNonzeroPositions(jfd);
    ASSERT_EQ(sub_entries.size(), nzs.size()) << "First block: wrong number of entries";
    for (const auto& [r, c] : nzs) {
      EXPECT_EQ(0, r % nvars) << "First block row variable wrong";
      EXPECT_EQ(0, c % nvars) << "First block col variable wrong";
    }
    MatDestroy(&jfd);
  }

  // Test last block (nvars-1, nvars-1)
  {
    auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);
    applyAndAssemble(jfd, sub, nvars - 1, nvars - 1);

    const auto nzs = matNonzeroPositions(jfd);
    ASSERT_EQ(sub_entries.size(), nzs.size()) << "Last block: wrong number of entries";
    for (const auto& [r, c] : nzs) {
      EXPECT_EQ(nvars - 1, r % nvars) << "Last block row variable wrong";
      EXPECT_EQ(nvars - 1, c % nvars) << "Last block col variable wrong";
    }
    MatDestroy(&jfd);
  }

  MatDestroy(&sub);
}

// ============================================================================
// all_pairs_single_variable
//
// With nvars=1 the all-pairs overload is identical to addOperatorSparsity
// with out_var=0, in_var=0.
// ============================================================================
TEST_F(PetscAddSparsityTest, all_pairs_single_variable) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 1;

  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {{0, 0}, {1, 1}, {2, 2}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  addOperatorSparsity(jfd, sub);
  MatAssemblyBegin(jfd, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(jfd, MAT_FINAL_ASSEMBLY);

  const auto nzs = matNonzeroPositions(jfd);
  ASSERT_EQ(sub_entries.size(), nzs.size());

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// all_pairs_fills_every_variable_block
//
// With nvars=3 every one of the 9 variable blocks must contain exactly the
// same pattern as sub.  The total nonzero count must equal
// nvars * nvars * nnz(sub).
// ============================================================================
TEST_F(PetscAddSparsityTest, all_pairs_fills_every_variable_block) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 3;

  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {{0, 0}, {1, 1}, {2, 2}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  addOperatorSparsity(jfd, sub);
  MatAssemblyBegin(jfd, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(jfd, MAT_FINAL_ASSEMBLY);

  const auto nzs = matNonzeroPositions(jfd);
  const std::size_t expected_total =
      static_cast<std::size_t>(nvars * nvars) * sub_entries.size();
  ASSERT_EQ(expected_total, nzs.size());

  // Every (out_var, in_var) pair must appear exactly once per sub entry.
  // Count how many Jfd entries land in each variable block.
  std::vector<std::vector<int>> block_counts(nvars, std::vector<int>(nvars, 0));
  for (const auto& [r, c] : nzs) {
    block_counts[r % nvars][c % nvars]++;
  }
  for (int ov = 0; ov < nvars; ++ov) {
    for (int iv = 0; iv < nvars; ++iv) {
      EXPECT_EQ(static_cast<int>(sub_entries.size()), block_counts[ov][iv])
          << "Block (" << ov << "," << iv << ") has wrong entry count";
    }
  }

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// all_pairs_entries_match_single_pair_calls
//
// The all-pairs overload must produce exactly the same Jfd nonzero positions
// as calling the single-pair overload for every (out_var, in_var).
// ============================================================================
TEST_F(PetscAddSparsityTest, all_pairs_entries_match_single_pair_calls) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 2;

  // Non-trivial pattern so the stride mapping is exercised
  const std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {
      {0, 0}, {0, 2}, {1, 1}, {2, 0}, {2, 2}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  // Build expected positions by calling the single-pair version for each block
  std::vector<std::pair<PetscInt, PetscInt>> expected;
  {
    auto [jfd_ref, nlocal_ref, rstart_ref] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);
    for (int ov = 0; ov < nvars; ++ov) {
      for (int iv = 0; iv < nvars; ++iv) {
        addOperatorSparsity(jfd_ref, sub, ov, iv);
      }
    }
    MatAssemblyBegin(jfd_ref, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jfd_ref, MAT_FINAL_ASSEMBLY);
    expected = matNonzeroPositions(jfd_ref);
    MatDestroy(&jfd_ref);
  }

  // Build actual positions using the all-pairs overload
  std::vector<std::pair<PetscInt, PetscInt>> actual;
  {
    auto [jfd_act, nlocal_act, rstart_act] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);
    addOperatorSparsity(jfd_act, sub);
    MatAssemblyBegin(jfd_act, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jfd_act, MAT_FINAL_ASSEMBLY);
    actual = matNonzeroPositions(jfd_act);
    MatDestroy(&jfd_act);
  }

  std::sort(expected.begin(), expected.end());
  std::sort(actual.begin(), actual.end());

  ASSERT_EQ(expected.size(), actual.size());
  for (std::size_t k = 0; k < expected.size(); ++k) {
    EXPECT_EQ(expected[k], actual[k]) << "Mismatch at position " << k;
  }

  MatDestroy(&sub);
}

// ============================================================================
// all_pairs_empty_sub_fills_nothing
//
// An empty sub must leave Jfd empty regardless of nvars.
// ============================================================================
TEST_F(PetscAddSparsityTest, all_pairs_empty_sub_fills_nothing) {
  const PetscInt n_sub = 4;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 3;

  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, {});

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  addOperatorSparsity(jfd, sub);
  MatAssemblyBegin(jfd, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(jfd, MAT_FINAL_ASSEMBLY);

  const auto nzs = matNonzeroPositions(jfd);
  EXPECT_EQ(0U, nzs.size());

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

// ============================================================================
// all_pairs_no_entries_outside_expected_positions
//
// Every entry in Jfd must lie at a position consistent with the stride formula
// applied to some (out_var, in_var) pair.  Specifically, for each entry
// (r, c) in Jfd there must exist a (sr, sc) in sub such that
//   r == sr * nvars + (r % nvars)  and  c == sc * nvars + (c % nvars).
// ============================================================================
TEST_F(PetscAddSparsityTest, all_pairs_no_entries_outside_expected_positions) {
  const PetscInt n_sub = 3;
  const PetscInt nlocal_sub = n_sub;
  const int nvars = 2;

  std::vector<std::pair<PetscInt, PetscInt>> sub_entries = {
      {0, 1}, {1, 0}, {1, 2}, {2, 2}};
  Mat sub = buildPatternMat(n_sub, nlocal_sub, 0, sub_entries);

  // Pre-compute sub cell pairs for fast lookup
  const std::set<std::pair<PetscInt, PetscInt>> sub_set(sub_entries.begin(),
                                                        sub_entries.end());

  auto [jfd, jfd_nlocal, jfd_rstart] = buildEmptyJfd(n_sub, nlocal_sub, 0, nvars);

  addOperatorSparsity(jfd, sub);
  MatAssemblyBegin(jfd, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(jfd, MAT_FINAL_ASSEMBLY);

  const auto nzs = matNonzeroPositions(jfd);
  for (const auto& [r, c] : nzs) {
    const PetscInt sr = r / nvars;
    const PetscInt sc = c / nvars;
    EXPECT_TRUE(sub_set.count({sr, sc}) > 0)
        << "Entry (" << r << "," << c << ") maps to sub cell (" << sr << "," << sc
        << ") which is not in sub";
  }

  MatDestroy(&sub);
  MatDestroy(&jfd);
}

#endif // BOUT_HAS_PETSC
