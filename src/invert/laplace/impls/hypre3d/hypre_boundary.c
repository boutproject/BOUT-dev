
/*
 * This function modifies the input for the HYPRE_IJMatrixSetValues() routine to
 * eliminate the boundary condition equations (see below for details on how the
 * equations are adjusted).  It modifies the arrays ncols, rows, cols, and
 * values.  It also returns a row_indexes array.  This can then be passed to the
 * HYPRE_IJMatrixSetValues2() routine to set up the matrix in hypre.
 *
 * The arguments nb and bi_array indicate the boundary equations.  The routine
 * returns info needed to adjust the right-hand-side and solution vector through
 * the functions AdjustRightHandSideEquations and AdjustSolutionEquations.
 *
 * NOTE: It may make sense from an organizational standpoint to collect many of
 * these arguments in a structure of some sort.
 *
 * Notation, assumptions, and other details:
 *
 * - Boundary equation i is assumed to have two coefficients
 *
 *      b_ii * u_i + b_ij * u_j = rhs_i
 *
 * - We also assume that each boundary equation has only one interior equation k
 *   coupled to it (such that k = j) with coupling coefficient a_ki
 *
 *      a_ki * u_i + a_kj * u_j + ... = rhs_k
 *
 * - Each equation k is adjusted as follows:
 *
 *      a_kj = a_kj - a_ki * b_ij / b_ii
 *      a_ki = 0
 *
 * - Boundary equations are adjusted to be identity equations in the matrix, but
 *   the boundary coefficients (b_ii, b_ij) are returned for use later
 *
 * - Right-hand-side equations are adjusted in AdjustRightHandSideEquations() as
 *   follows: rhs_k = rhs_k - a_ki * rhs_i / b_ii
 *
 * - Solution unknowns are adjusted at boundaries in AdjustSolutionEquations as
 *   follows: u_i = (rhs_i - b_ij * u_j) / b_ii
 *
 * - Naming conventions: Arrays starting with 'b' are boundary equation arrays
 *   indexed by 'bnum', and arrays starting with 'a' are non-boundary arrays
 *   (interior matrix equations) indexed by 'anum'.  When 'num' is prefixed with
 *   a row or column number 'i', 'j', or 'k', the array holds the corresponding
 *   local data index for that row or column (e.g., an index into the local
 *   solution vector).  Matrix coefficients are named as above, e.g., 'bij' is
 *   the coefficient for b_ij.
 */

void
AdjustBCMatrixEquations(
   HYPRE_Int       nrows,
   HYPRE_Int      *ncols,
   HYPRE_BigInt   *rows,
   HYPRE_Int     **row_indexes_ptr,
   HYPRE_BigInt   *cols,
   HYPRE_Complex  *values,
   HYPRE_Int       nb,              // number of boundary equations
   HYPRE_Int      *bi_array,        // row i for each boundary equation
   HYPRE_Int     **binum_array_ptr, // data index for row i (for each boundary equation)
   HYPRE_Int     **bjnum_array_ptr, // data index for col j (for each boundary equation)
   HYPRE_Complex **bii_array_ptr,   // coefficient b_ii (for each boundary equation)
   HYPRE_Complex **bij_array_ptr,   // coefficient b_ij (for each boundary equation)
   HYPRE_Int      *na_ptr,          // number of interior equations to adjust
   HYPRE_Int     **aknum_array_ptr, // data index for row k (for each interior equation)
   HYPRE_Complex **aki_array_ptr)   // coefficient a_ki (for each interior equation)
{
   HYPRE_Int     *row_indexes;
   HYPRE_Int      na, *binum_array, *bjnum_array, *aknum_array;
   HYPRE_Complex *bii_array, *bij_array, *aki_array;
   HYPRE_Int      i, j, k, m, mkj, anum, bnum, acoeffnum, bcoeffnum;
   HYPRE_Int      binum, aknum;
   HYPRE_Complex  bii, bij, aki;

   /* Create the row_indexes array */
   row_indexes = (HYPRE_Int *)malloc(sizeof(HYPRE_Int) * nrows);
   row_indexes[0] = 0;
   for (i = 1; i < nrows; i++)
   {
      row_indexes[i] = row_indexes[i-1] + ncols[i-1];
   }

   /* Assume just one interior equation coupled to each boundary equation */
   na = nb;

   /* Allocate return arrays */
   HypreMalloc(binum_array, sizeof(HYPRE_Int) * nb);
   HypreMalloc(bjnum_array, sizeof(HYPRE_Int) * nb);
   HypreMalloc(bii_array,   sizeof(HYPRE_Complex) * nb);
   HypreMalloc(bij_array,   sizeof(HYPRE_Complex) * nb);
   HypreMalloc(aknum_array, sizeof(HYPRE_Int) * na);
   HypreMalloc(aki_array,   sizeof(HYPRE_Complex) * na);

   binum = 0;
   aknum = 0;
   for (bnum = 0; bnum < nb; bnum++)
   {
      /* Get boundary equation information and adjust boundary equations */
      /* Find row i in rows array (assume i increases and rows is sorted) */
      i = bi_array[bnum];
      for (; binum < nrows; binum++)
      {
         if (i == rows[binum])
         {
            break;   // Found row i in rows array
         }
      }
      bcoeffnum = row_indexes[binum];
      for (m = 0; m < 2; m++)               // Assume only two boundary equation coefficients
      {
         if (cols[bcoeffnum + m] == i)
         {
            bii = values[bcoeffnum + m];
            values[bcoeffnum + m] = -1.0;   // Identity equation (negative definite matrix)
         }
         else
         {
            j = cols[bcoeffnum + m];
            bij = values[bcoeffnum + m];
            values[bcoeffnum + m] = 0.0;   // Identity equation
         }
      }
      ncols[binum] = 1;                  // Identity equation

      /* Get interior equation information and adjust interior equations */
      /* Find row k in rows array (assume k increases and rows is sorted) */
      k = j;         // Assume equation k = j
      for (; aknum < nrows; aknum++)
      {
         if (k == rows[aknum])
         {
            break;   // Found row k in rows array
         }
      }
      acoeffnum = row_indexes[aknum];
      for (m = 0; m < ncols[aknum]; m++)
      {
         if (cols[acoeffnum + m] == j)
         {
            mkj = m;                       // Save for update of akj value below
         }
         if (cols[acoeffnum + m] == i)
         {
            aki = values[acoeffnum + m];
            values[acoeffnum + m] = 0.0;   // Eliminate coupling to boundary equation
         }
      }
      values[acoeffnum + mkj] -= aki * bij / bii;   // Update akj value

      /* Update return arrays */
      anum = bnum;                // Assume only one interior equation k
      binum_array[bnum] = binum;
      bjnum_array[bnum] = aknum;  // Assume only one interior equation k
      bii_array[bnum]   = bii;
      bij_array[bnum]   = bij;
      aknum_array[anum] = aknum;
      aki_array[anum]   = aki;
   }

   /* Set return arguments */
   *row_indexes_ptr = row_indexes;
   *binum_array_ptr = binum_array;
   *bjnum_array_ptr = bjnum_array;
   *bii_array_ptr   = bii_array;
   *bij_array_ptr   = bij_array;
   *na_ptr          = na;
   *aknum_array_ptr = aknum_array;
   *aki_array_ptr   = aki_array;
}

void
AdjustBCRightHandSideEquations(
   HYPRE_Complex  *rhs,
   HYPRE_Int       nb,
   HYPRE_Int      *binum_array,
   HYPRE_Complex  *bii_array,
   HYPRE_Complex  *bij_array,
   HYPRE_Complex **brhs_array_ptr,
   HYPRE_Int       na,
   HYPRE_Int      *aknum_array,
   HYPRE_Complex  *aki_array)
{
   HYPRE_Complex *brhs_array;
   HYPRE_Int      anum, bnum, binum, aknum;

   HypreMalloc(brhs_array, sizeof(HYPRE_Complex) * nb);

   for (bnum = 0; bnum < nb; bnum++)
   {
      binum = binum_array[bnum];
      brhs_array[bnum] = rhs[binum];
   }

   for (anum = 0; anum < na; anum++)
   {
      bnum  = anum;   // Assume only one interior equation per boundary equation
      aknum = aknum_array[anum];
      rhs[aknum] -= aki_array[anum] * bij_array[bnum] / bii_array[bnum];
   }

   *brhs_array_ptr = brhs_array;
}

void
AdjustBCSolutionEquations(
   HYPRE_Complex  *solution,
   HYPRE_Int       nb,
   HYPRE_Int      *binum_array,
   HYPRE_Int      *bjnum_array,
   HYPRE_Complex  *bii_array,
   HYPRE_Complex  *bij_array,
   HYPRE_Complex  *brhs_array)
{
   HYPRE_Int  bnum, binum, bjnum;

   for (bnum = 0; bnum < nb; bnum++)
   {
      binum = binum_array[bnum];
      bjnum = bjnum_array[bnum];
      solution[binum] = (brhs_array[bnum] - bij_array[bnum] * solution[bjnum]) / bii_array[bnum];
   }
}

