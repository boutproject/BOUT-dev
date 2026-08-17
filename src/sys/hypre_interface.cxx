
#include "bout/build_defines.hxx"

#if BOUT_HAS_HYPRE

#include "bout/hypre_interface.hxx"

#include <numeric>
#include <unordered_map>

namespace bout {

BoundaryElimination::BoundaryElimination(HYPRE_Int nrows, HYPRE_Int* ncols,
                                         HYPRE_BigInt* rows, HYPRE_Int** row_indexes_ptr,
                                         HYPRE_BigInt* cols, HYPRE_Complex* values,
                                         HYPRE_Int nb, HYPRE_Int* bi_array)
    : nb(nb) {
  HYPRE_Int* row_indexes;
  const auto find_local_row = [nrows, rows](HYPRE_BigInt row) -> HYPRE_Int {
    auto row_position = std::lower_bound(rows, rows + nrows, row);
    if ((row_position == rows + nrows) || (*row_position != row)) {
      throw BoutException("Could not find local row {} while constructing boundary "
                          "elimination data",
                          row);
    }
    return static_cast<HYPRE_Int>(std::distance(rows, row_position));
  };

  // Create the row_indexes array
  HypreMalloc(row_indexes, sizeof(HYPRE_Int) * nrows);
  row_indexes[0] = 0;
  for (HYPRE_Int i = 1; i < nrows; i++) {
    row_indexes[i] = row_indexes[i - 1] + ncols[i - 1];
  }

  // Allocate arrays
  HypreMalloc(binum_array, sizeof(HYPRE_Int) * nb);
  HypreMalloc(bjnum_array, sizeof(HYPRE_Int) * nb);
  HypreMalloc(bdep_array, sizeof(HYPRE_Int) * nb);
  HypreMalloc(bii_array, sizeof(HYPRE_Complex) * nb);
  HypreMalloc(bij_array, sizeof(HYPRE_Complex) * nb);

  std::vector<HYPRE_Int> boundary_at_row(nrows, -1);
  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    boundary_at_row[find_local_row(bi_array[bnum])] = bnum;
  }

  struct BoundaryOccurrence {
    HYPRE_Int rownum;
    HYPRE_Int row_offset;
  };
  std::unordered_map<HYPRE_BigInt, HYPRE_Int> boundary_number_for_row;
  boundary_number_for_row.reserve(nb);
  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    boundary_number_for_row.emplace(bi_array[bnum], bnum);
  }

  std::vector<std::vector<BoundaryOccurrence>> occurrences_by_boundary(nb);
  for (HYPRE_Int rownum = 0; rownum < nrows; rownum++) {
    const HYPRE_Int row_start = row_indexes[rownum];
    for (HYPRE_Int m = 0; m < ncols[rownum]; m++) {
      const auto boundary_position = boundary_number_for_row.find(cols[row_start + m]);
      if (boundary_position != boundary_number_for_row.end()) {
        occurrences_by_boundary[boundary_position->second].push_back({rownum, m});
      }
    }
  }

  std::vector<HYPRE_BigInt> bjrow_array(nb, -1);
  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    // Get boundary equation information and adjust boundary equations
    HYPRE_Int i = bi_array[bnum];
    const HYPRE_Int binum = find_local_row(i);
    HYPRE_Int bcoeffnum = row_indexes[binum];
    HYPRE_Complex bii{0.0}, bij{0.0};
    HYPRE_BigInt j = -1;

    for (HYPRE_Int m = 0; m < ncols[binum]; m++) {
      if (cols[bcoeffnum + m] == i) {
        bii = values[bcoeffnum + m];
        values[bcoeffnum + m] = -1.0; // Identity equation (negative definite matrix)
      } else {
        j = cols[bcoeffnum + m];
        bij = values[bcoeffnum + m];
        values[bcoeffnum + m] = 0.0; // Identity equation
      }
    }
    if (j < 0) {
      throw BoutException("Boundary row {} does not contain a retained neighbour", i);
    }
    if (bii == 0.0) {
      throw BoutException("Boundary row {} has zero diagonal coefficient", i);
    }
    ncols[binum] = 1; // Identity equation

    // Update arrays
    binum_array[bnum] = binum;
    bjnum_array[bnum] = find_local_row(j);
    bdep_array[bnum] = boundary_at_row[bjnum_array[bnum]];
    bii_array[bnum] = bii;
    bij_array[bnum] = bij;
    bjrow_array[bnum] = j;
  }

  std::vector<HYPRE_Int> aoffsets(nb + 1, 0);
  std::vector<HYPRE_Int> aknums;
  std::vector<HYPRE_Complex> akis;
  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    const HYPRE_Int binum = binum_array[bnum];
    const HYPRE_Int bjnum = bjnum_array[bnum];
    const HYPRE_BigInt j = bjrow_array[bnum];

    for (const auto& occurrence : occurrences_by_boundary[bnum]) {
      const HYPRE_Int aknum = occurrence.rownum;
      if (aknum == binum) {
        continue;
      }
      if ((boundary_at_row[aknum] >= 0) && (aknum != bjnum)) {
        continue;
      }

      const HYPRE_Int acoeffnum = row_indexes[aknum];
      const HYPRE_Int aki_position = acoeffnum + occurrence.row_offset;
      const HYPRE_Complex aki = values[aki_position];
      if (aki == 0.0) {
        continue;
      }

      HYPRE_Int j_position = -1;
      for (HYPRE_Int m = 0; m < ncols[aknum]; m++) {
        if (cols[acoeffnum + m] == j) {
          j_position = acoeffnum + m;
          break;
        }
      }
      if (j_position < 0) {
        continue;
      }

      values[aki_position] = 0.0; // Eliminate coupling to boundary equation
      values[j_position] -= aki * bij_array[bnum] / bii_array[bnum];

      aknums.push_back(aknum);
      akis.push_back(aki);
    }
    aoffsets[bnum + 1] = static_cast<HYPRE_Int>(aknums.size());
  }

  na = static_cast<HYPRE_Int>(aknums.size());
  HypreMalloc(aoffset_array, sizeof(HYPRE_Int) * (nb + 1));
  std::copy(aoffsets.begin(), aoffsets.end(), aoffset_array);
  HypreMalloc(aknum_array, sizeof(HYPRE_Int) * na);
  std::copy(aknums.begin(), aknums.end(), aknum_array);
  HypreMalloc(aki_array, sizeof(HYPRE_Complex) * na);
  std::copy(akis.begin(), akis.end(), aki_array);

  std::vector<HYPRE_Int> dependency_depth(nb, -1);
  const auto get_depth = [&dependency_depth, this](const auto& self,
                                                   HYPRE_Int bnum) -> HYPRE_Int {
    if (dependency_depth[bnum] >= 0) {
      return dependency_depth[bnum];
    }
    const HYPRE_Int dep = bdep_array[bnum];
    if (dep < 0) {
      dependency_depth[bnum] = 0;
    } else {
      dependency_depth[bnum] = 1 + self(self, dep);
    }
    return dependency_depth[bnum];
  };
  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    get_depth(get_depth, bnum);
  }

  expansion_order.resize(nb);
  std::iota(expansion_order.begin(), expansion_order.end(), 0);
  std::sort(expansion_order.begin(), expansion_order.end(),
            [&dependency_depth](HYPRE_Int lhs, HYPRE_Int rhs) {
              if (dependency_depth[lhs] != dependency_depth[rhs]) {
                return dependency_depth[lhs] < dependency_depth[rhs];
              }
              return lhs < rhs;
            });
  reduction_order = expansion_order;
  std::reverse(reduction_order.begin(), reduction_order.end());

  // Set return arguments
  *row_indexes_ptr = row_indexes;
}

BCValuesPtr
BoundaryElimination::copyBoundaryRowValues(const HYPRE_Complex* values) const {
  BCValuesPtr boundary_values = std::make_shared<HypreComplexArray>(nb);

  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    boundary_values->data[bnum] = values[binum_array[bnum]];
  }

  return boundary_values;
}

BCValuesPtr
BoundaryElimination::evaluateBoundaryEquations(const HYPRE_Complex* values) const {
  BCValuesPtr boundary_values = std::make_shared<HypreComplexArray>(nb);

  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    boundary_values->data[bnum] = (bii_array[bnum] * values[binum_array[bnum]])
                                  + (bij_array[bnum] * values[bjnum_array[bnum]]);
  }

  return boundary_values;
}

BCValuesPtr BoundaryElimination::reduceRightHandSideInPlace(HYPRE_Complex* rhs) const {

  // Allocate array to store boundary row values
  BCValuesPtr brhs = copyBoundaryRowValues(rhs);

  for (HYPRE_Int bnum : reduction_order) {
    for (HYPRE_Int anum = aoffset_array[bnum]; anum < aoffset_array[bnum + 1]; anum++) {
      HYPRE_Int aknum = aknum_array[anum];
      rhs[aknum] -= aki_array[anum] * brhs->data[bnum] / bii_array[bnum];
    }
    if (bdep_array[bnum] >= 0) {
      brhs->data[bdep_array[bnum]] = rhs[bjnum_array[bnum]];
    }
  }

  return brhs;
}

void BoundaryElimination::expandSolutionInPlace(BCValuesPtr brhs,
                                                HYPRE_Complex* solution) const {

  for (HYPRE_Int bnum : expansion_order) {
    HYPRE_Int binum = binum_array[bnum];
    HYPRE_Int bjnum = bjnum_array[bnum];
    solution[binum] =
        (brhs->data[bnum] - bij_array[bnum] * solution[bjnum]) / bii_array[bnum];
  }
}

void BoundaryElimination::expandMatvecResultInPlace(BCValuesPtr boundary_operator_values,
                                                    BCValuesPtr full_boundary_values,
                                                    HYPRE_Complex* result) const {
  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    for (HYPRE_Int anum = aoffset_array[bnum]; anum < aoffset_array[bnum + 1]; anum++) {
      HYPRE_Int aknum = aknum_array[anum];
      result[aknum] +=
          aki_array[anum] * boundary_operator_values->data[bnum] / bii_array[bnum];
    }
  }

  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    result[binum_array[bnum]] = full_boundary_values->data[bnum];
  }
}

} // namespace bout

#endif // BOUT_HAS_HYPRE
