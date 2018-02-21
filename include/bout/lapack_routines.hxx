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

#ifndef __LAPACK_ROUTINES_H__
#define __LAPACK_ROUTINES_H__

/* Tridiagonal inversion
 *
 * a = Left of diagonal (so a[0] not used)
 * b = diagonal terms
 * c = Right of diagonal (so c[n-1] not used)
 * r = RHS vector
 * u = Result vector
 * n = Length of vectors
 *
 * 
 * b[0]   c[0]
 * a[1]   b[1]   c[1]
 *         .       .    . 
 * 
 */

// Tri-diagonal solvers
int tridag(const dcomplex *a, const dcomplex *b, const dcomplex *c, const dcomplex *r, dcomplex *u, int n);
bool tridag(const BoutReal *a, const BoutReal *b, const BoutReal *c, const BoutReal *r, BoutReal *x, int n);

// Cyclic tridiagonal
void cyclic_tridag(BoutReal *a, BoutReal *b, BoutReal *c, BoutReal *r, BoutReal *x, int n);
void cyclic_tridag(dcomplex *a, dcomplex *b, dcomplex *c, dcomplex *r, dcomplex *x, int n);

/// Complex band matrix solver
void cband_solve(dcomplex **a, int n, int m1, int m2, dcomplex *b);

#endif // __LAPACK_ROUTINES_H__

