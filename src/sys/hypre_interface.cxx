
#include "bout/build_defines.hxx"

#if BOUT_HAS_HYPRE

#include "bout/hypre_interface.hxx"

namespace bout {

BCMatrixEquations::BCMatrixEquations(HYPRE_Int nrows, HYPRE_Int* ncols,
                                     HYPRE_BigInt* rows, HYPRE_Int** row_indexes_ptr,
                                     HYPRE_BigInt* cols, HYPRE_Complex* values,
                                     HYPRE_Int nb, HYPRE_Int* bi_array)
    : nb(nb) {
  HYPRE_Int* row_indexes;

  // Create the row_indexes array
  row_indexes = (HYPRE_Int*)malloc(sizeof(HYPRE_Int) * nrows);
  row_indexes[0] = 0;
  for (HYPRE_Int i = 1; i < nrows; i++) {
    row_indexes[i] = row_indexes[i - 1] + ncols[i - 1];
  }

  // Assume just one interior equation coupled to each boundary equation
  na = nb;

  // Allocate arrays
  HypreMalloc(binum_array, sizeof(HYPRE_Int) * nb);
  HypreMalloc(bjnum_array, sizeof(HYPRE_Int) * nb);
  HypreMalloc(bii_array, sizeof(HYPRE_Complex) * nb);
  HypreMalloc(bij_array, sizeof(HYPRE_Complex) * nb);
  HypreMalloc(aknum_array, sizeof(HYPRE_Int) * na);
  HypreMalloc(aki_array, sizeof(HYPRE_Complex) * na);

  HYPRE_Int binum = 0;
  HYPRE_Int aknum = 0;
  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    // Get boundary equation information and adjust boundary equations
    // Find row i in rows array (assume i increases and rows is sorted)
    HYPRE_Int i = bi_array[bnum];
    for (; binum < nrows; binum++) {
      if (i == rows[binum]) {
        break; // Found row i in rows array
      }
    }
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
    /* Find row k in rows array (assume k increases and rows is sorted) */
    HYPRE_Int k = j; // Assume equation k = j
    for (; aknum < nrows; aknum++) {
      if (k == rows[aknum]) {
        break; // Found row k in rows array
      }
    }
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
    bii_array[bnum] = bii;
    bij_array[bnum] = bij;
    aknum_array[anum] = aknum;
    aki_array[anum] = aki;
  }

  // Set return arguments
  *row_indexes_ptr = row_indexes;
}

BCValuesPtr BCMatrixEquations::adjustBCRightHandSideEquations(HYPRE_Complex* rhs) {

  // Allocate array to store boundary row values
  BCValuesPtr brhs = std::make_shared<HypreComplexArray>(nb);

  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    HYPRE_Int binum = binum_array[bnum];
    brhs->data[bnum] = rhs[binum];
  }

  for (HYPRE_Int anum = 0; anum < na; anum++) {
    HYPRE_Int bnum = anum; // Assume only one interior equation per boundary equation
    HYPRE_Int aknum = aknum_array[anum];
    rhs[aknum] -= aki_array[anum] * brhs->data[bnum] / bii_array[bnum];
  }

  return brhs;
}

void BCMatrixEquations::adjustBCSolutionEquations(BCValuesPtr brhs,
                                                  HYPRE_Complex* solution) {

  for (HYPRE_Int bnum = 0; bnum < nb; bnum++) {
    HYPRE_Int binum = binum_array[bnum];
    HYPRE_Int bjnum = bjnum_array[bnum];
    solution[binum] =
        (brhs->data[bnum] - bij_array[bnum] * solution[bjnum]) / bii_array[bnum];
  }
}

} // namespace bout

#endif // BOUT_HAS_HYPRE
