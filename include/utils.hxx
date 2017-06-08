/*!*************************************************************************
 * \file utils.hxx
 *
 * A mix of short utilities for memory management, strings, and some
 * simple but common calculations
 *
 **************************************************************************
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

#ifndef __UTILS_H__
#define __UTILS_H__

#include "bout_types.hxx"
#include "dcomplex.hxx"
#include "boutexception.hxx"

#include "bout/deprecated.hxx"

#include <cmath>

using std::abs;
using std::swap;

/*!
 * Allocates an array of \p size BoutReals
 */
DEPRECATED(BoutReal *rvector(int size));

/*!
 * Resizes an array of BoutReals to \p newsize
 */ 
DEPRECATED(BoutReal *rvresize(BoutReal *v, int newsize));

/*!
 * Frees an array of BoutReals
 */
DEPRECATED(void rvfree(BoutReal *r));

/*!
 * Allocates an array of \p size ints
 */
DEPRECATED(int *ivector(int size));

/*!
 * Resizes an array of ints to \p newsize
 */
DEPRECATED(int *ivresize(int *v, int newsize));

/*!
 * Frees an array of ints
 */
DEPRECATED(void ivfree(int *v));

/*!
 * Allocate a 2D array of \p xsize by \p ysize BoutReals
 */
DEPRECATED(BoutReal **rmatrix(int xsize, int ysize));

/*!
 * Allocate a 2D array of \p xsize by \p ysize ints
 */
DEPRECATED(int **imatrix(int xsize, int ysize));

/*!
 * Create a 2D array of \p xsize by \p ysize 
 * This is allocated as two blocks of data so that
 * the values are in a contiguous array.
 * 
 * Note: This returns C-style pointers, and makes
 * no effort to manage memory. Prefer other methods
 * (like standard containers) over this if possible.
 * 
 * Example
 * -------
 * 
 * BoutReal **m = matrix<BoutReal>(nx, ny);
 */
template <class T>
T **matrix(int xsize, int ysize) {
  long i;
  T **m;

  if(xsize == 0)
     xsize = 1;
  if(ysize == 0)
     ysize = 1;

  if((m = new T*[xsize]) == NULL)
    throw BoutException("Error: could not allocate memory:%d\n", xsize);
  
  if((m[0] = new T[xsize*ysize]) == NULL)
    throw BoutException("Error: could not allocate memory\n");

  for(i=1;i<xsize;i++) {
    m[i] = m[i-1] + ysize;
  }
  return m;
}

/*!
 * Free a 2D array of BoutReals, assumed to have been allocated using rmatrix()
 */
DEPRECATED(void free_rmatrix(BoutReal **m));

/*!
 * Free a 2D array of ints, assumed to have been allocated using imatrix()
 */
DEPRECATED(void free_imatrix(int **m));

/*!
 * Free a matrix, assumed to have been allocated using matrix()
 *
 * @param[in] T  The matrix to free
 *
 * Example
 * -------
 *
 *     BoutReal **m = matrix<BoutReal>(nx, ny);
 *     ...
 *     free_matrix(m);
 */ 
template <class T>
void free_matrix(T **m) {
  delete[] m[0];
  delete[] m;
}

/*!
 * Allocate a 3D BoutReal array of size \p nrow x \p ncol \p ndep
 
 * Note: Prefer other methods like standard containers
 */ 
BoutReal ***r3tensor(int nrow, int ncol, int ndep);

/*!
 * Free a 3D BoutReal array, assumed to have been created
 * by r3tensor()
 *
 */
void free_r3tensor(BoutReal ***m);

/*!
 * Allocate a 3D int array of size \p nrow x \p ncol \p ndep
 
 * Note: Prefer other methods like standard containers
 */ 
int ***i3tensor(int nrow, int ncol, int ndep);

/*!
 * Free a 3D int array, assumed to have been created
 * by i3tensor()
 */
void free_i3tensor(int ***m);

/*!
 * Allocate a 2D array of \p nrow by \p ncol dcomplex objects
 */
DEPRECATED(dcomplex **cmatrix(int nrow, int ncol));

/*!
 * Free a 2D array, assumed to have been allocated using cmatrix
 */
DEPRECATED(void free_cmatrix(dcomplex** cm));

/*!
 * Get Random number between 0 and 1
 */
inline BoutReal randomu() {
  return ((BoutReal) rand()) / ((BoutReal) RAND_MAX);
}

/*!
 * Calculate the square of a variable \p t
 * i.e. t * t
 */
template <typename T>
T SQ(T t){
  return t*t;
}

/*!
 * Round \p x to the nearest integer
 */
inline int ROUND(BoutReal x){
  return (x > 0.0) ? (int) (x + 0.5) : (int) (x - 0.5);
}

/*!
 * Calculate the maximum of a list of values
 * using a > b operator
 */
template <typename T>
T BOUTMAX(T a){
  return a;
}
template <typename T, typename... Args>
T BOUTMAX(T a,T b,Args... args){
  T c = BOUTMAX(b,args...);
  return c > a ? c : a;
}

/*!
 * Calculate the minimum of a list of values
 * using the a < b operator
 */
template <typename T>
T BOUTMIN(T a){
  return a;
}
template <typename T, typename... Args>
T BOUTMIN(T a,T b,Args... args){
  T c = BOUTMIN(b,args...);
  return c < a ? c : a;
}

/*!
 * Check if a number is a power of 2
 */ 
inline bool is_pow2(int x) {
  return x && !((x-1) & x);
}

/*!
 * Return the sign of a number \p a
 * by testing if a > 0 
 */
template <typename T>
T SIGN(T a) { // Return +1 or -1 (0 -> +1)
  return a > 0 ? 1 : -1;
}

/*!
 * The minimum absolute value of \p a and \p b
 *
 * if \p a and \p b have opposite signs, return zero
 *
 * if |a| < |b| then return a, otherwise return b
 */
inline BoutReal MINMOD(BoutReal a, BoutReal b) {
  return 0.5*(SIGN(a) + SIGN(b)) * BOUTMIN(fabs(a), fabs(b));
}

#endif // __UTILS_H__
