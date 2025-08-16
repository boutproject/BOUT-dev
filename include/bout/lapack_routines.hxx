/// \file
///
/// Serial code to invert a complex tridiagonal system
///
/// Solves a banded matrix given the matrix in compact form:
///
///     a[0...(n-1)][0...(m1+m2)]
///
/// and the rhs vector:
///
///     b[0...(n-1)]
///
/// ``a`` is overwritten, and ``b`` is replaced by the solution

/**************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef BOUT_LAPACK_ROUTINES_H
#define BOUT_LAPACK_ROUTINES_H

#include "bout/dcomplex.hxx"
#include "bout/utils.hxx"

/// Tridiagonal inversion
///
/// \param a Left of diagonal (so a[0] not used)
/// \param b diagonal terms
/// \param c Right of diagonal (so c[n-1] not used)
/// \param r RHS vector
/// \param u Result vector
/// \param n Length of vectors
///
///     b[0]   c[0]
///     a[1]   b[1]   c[1]
///             .       .    .
///
/// Uses LAPACK routine ZGTSV
int tridag(const dcomplex* a, const dcomplex* b, const dcomplex* c, const dcomplex* r,
           dcomplex* u, int n);
bool tridag(const BoutReal* a, const BoutReal* b, const BoutReal* c, const BoutReal* r,
            BoutReal* x, int n);

/// Cyclic tridiagonal
///
/// Uses Sherman-Morrison formula
void cyclic_tridag(BoutReal* a, BoutReal* b, BoutReal* c, BoutReal* r, BoutReal* x,
                   int n);
void cyclic_tridag(dcomplex* a, dcomplex* b, dcomplex* c, dcomplex* r, dcomplex* x,
                   int n);

/// Complex band matrix solver using ZGBSV
void cband_solve(Matrix<dcomplex>& a, int n, int m1, int m2, Array<dcomplex>& b);

#endif // BOUT_LAPACK_ROUTINES_H
