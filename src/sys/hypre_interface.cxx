
#include "bout/build_defines.hxx"

#if BOUT_HAS_HYPRE

#include "bout/hypre_interface.hxx"

#include <numeric>

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

  // Assume just one interior equation coupled to each boundary equation
  na = nb;

  // Allocate arrays
  HypreMalloc(binum_array, sizeof(HYPRE_Int) * nb);
  HypreMalloc(bjnum_array, sizeof(HYPRE_Int) * nb);
  HypreMalloc(bdep_array, sizeof(HYPRE_Int) * nb);
  HypreMalloc(bii_array, sizeof(HYPRE_Complex) * nb);
  HypreMalloc(bij_array, sizeof(HYPRE_Complex) * nb);
  HypreMalloc(aknum_array, sizeof(HYPRE_Int) * na);
  HypreMalloc(aki_array, sizeof(HYPRE_Complex) * na);

  std::vector<HYPRE_Int> boundary_at_row(nrows, -1);
  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    boundary_at_row[find_local_row(bi_array[bnum])] = bnum;
  }

  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    // Get boundary equation information and adjust boundary equations
    HYPRE_Int i = bi_array[bnum];
    const HYPRE_Int binum = find_local_row(i);
    HYPRE_Int bcoeffnum = row_indexes[binum];
    HYPRE_Complex bii{0.0}, bij{0.0};
    HYPRE_Int j = 0;

    for (HYPRE_Int m = 0; m < 2; m++) { // Assume only two boundary equation coefficients
      if (cols[bcoeffnum + m] == i) {
        bii = values[bcoeffnum + m];
        values[bcoeffnum + m] = -1.0; // Identity equation (negative definite matrix)
      } else {
        j = cols[bcoeffnum + m];
        bij = values[bcoeffnum + m];
        values[bcoeffnum + m] = 0.0; // Identity equation
      }
    }
    ncols[binum] = 1; // Identity equation

    /* Get interior equation information and adjust interior equations */
    HYPRE_Int k = j; // Assume equation k = j
    const HYPRE_Int aknum = find_local_row(k);
    HYPRE_Int acoeffnum = row_indexes[aknum];

    HYPRE_Int mkj = 0;
    HYPRE_Complex aki{0.0};
    for (HYPRE_Int m = 0; m < ncols[aknum]; m++) {
      if (cols[acoeffnum + m] == j) {
        mkj = m; // Save for update of akj value below
      }
      if (cols[acoeffnum + m] == i) {
        aki = values[acoeffnum + m];
        values[acoeffnum + m] = 0.0; // Eliminate coupling to boundary equation
      }
    }
    values[acoeffnum + mkj] -= aki * bij / bii; // Update akj value

    // Update arrays
    HYPRE_Int anum = bnum; // Assume only one interior equation k
    binum_array[bnum] = binum;
    bjnum_array[bnum] = aknum; // Assume only one interior equation k
    bdep_array[bnum] = boundary_at_row[aknum];
    bii_array[bnum] = bii;
    bij_array[bnum] = bij;
    aknum_array[anum] = aknum;
    aki_array[anum] = aki;
  }

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
    HYPRE_Int anum = bnum; // Assume only one interior equation per boundary equation
    HYPRE_Int aknum = aknum_array[anum];
    rhs[aknum] -= aki_array[anum] * brhs->data[bnum] / bii_array[bnum];
    if (bdep_array[bnum] >= 0) {
      brhs->data[bdep_array[bnum]] = rhs[aknum];
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
  for (HYPRE_Int anum = 0; anum < na; anum++) {
    HYPRE_Int bnum = anum; // Assume only one interior equation per boundary equation
    HYPRE_Int aknum = aknum_array[anum];
    result[aknum] +=
        aki_array[anum] * boundary_operator_values->data[bnum] / bii_array[bnum];
  }

  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    result[binum_array[bnum]] = full_boundary_values->data[bnum];
  }
}

} // namespace bout

#endif // BOUT_HAS_HYPRE
