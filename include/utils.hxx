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

#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/deprecated.hxx"
#include "msg_stack.hxx"
#include "unused.hxx"

#include <string>
#include <list>
#include <cmath>
#include <algorithm>

using std::abs;
using std::swap;

/// Helper class for 2D arrays
///
/// Allows bounds checking through `operator()` with CHECK > 1
template <typename T>
class Matrix {
public:
  typedef T data_type;
  Matrix() : n1(0), n2(0){};
  Matrix(unsigned int n1, unsigned int n2) : n1(n1), n2(n2) {
    data = Array<T>(n1*n2);
  }
  Matrix(const Matrix &other) : n1(other.n1), n2(other.n2), data(other.data) {
    // Prevent copy on write for Matrix
    data.ensureUnique();
  }

  Matrix& operator=(const Matrix &other) {
    n1 = other.n1;
    n2 = other.n2;
    data = other.data;
    // Prevent copy on write for Matrix
    data.ensureUnique();
    return *this;
  }
  
  T& operator()(unsigned int i1, unsigned int i2) {
    ASSERT2(0<=i1 && i1<n1);
    ASSERT2(0<=i2 && i2<n2);
    return data[i1*n2+i2];
  }
  const T& operator()(unsigned int i1, unsigned int i2) const {
    ASSERT2(0<=i1 && i1<n1);
    ASSERT2(0<=i2 && i2<n2);
    return data[i1*n2+i2];
  }

  Matrix& operator=(const T&val){
    for (auto &i: data) {
      i = val;
    };
    return *this;
  };
  
  // To provide backwards compatibility with matrix to be removed
  DEPRECATED(T* operator[](unsigned int i1)) {
    ASSERT2(0<=i1 && i1<n1);
    return &(data[i1*n2]);
  }
  // To provide backwards compatibility with matrix to be removed
  DEPRECATED(const T* operator[](unsigned int i1) const) {
    ASSERT2(0<=i1 && i1<n1);
    return &(data[i1*n2]);
  }

  T* begin() { return std::begin(data);};
  const T* begin() const { return std::begin(data);};
  T* end() { return std::end(data);};
  const T* end() const { return std::end(data);};

  std::tuple<unsigned int, unsigned int> shape() { return std::make_tuple(n1, n2);};

  bool empty(){
    return n1*n2 == 0;
  }

  /*!
   * Ensures that this Matrix does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() {
    data.ensureUnique();
  }
  
private:
  unsigned int n1, n2;
  Array<T> data;
};

// For backwards compatibility with old matrix -- to be removed
template <typename T>
DEPRECATED(void free_matrix(Matrix<T> UNUSED(m)));
template <typename T>
void free_matrix(Matrix<T> UNUSED(m)) {};

/// Helper class for 3D arrays
///
/// Allows bounds checking through `operator()` with CHECK > 1
template <typename T>
class Tensor {
public:
  typedef T data_type;
  Tensor() : n1(0), n2(0), n3(0) {};
  Tensor(unsigned int n1, unsigned int n2, unsigned int n3) : n1(n1), n2(n2), n3(n3) {
    data = Array<T>(n1*n2*n3);
  }
  Tensor(const Tensor &other) : n1(other.n1), n2(other.n2), n3(other.n3), data(other.data) {
    // Prevent copy on write for Tensor
    data.ensureUnique();
  }

  Tensor& operator=(const Tensor &other) {
    n1 = other.n1;
    n2 = other.n2;
    n3 = other.n3;
    data = other.data;
    // Prevent copy on write for Tensor
    data.ensureUnique();
    return *this;
  }

  T& operator()(unsigned int i1, unsigned int i2, unsigned int i3) {
    ASSERT2(0<=i1 && i1<n1);
    ASSERT2(0<=i2 && i2<n2);
    ASSERT2(0<=i3 && i3<n3);
    return data[(i1*n2+i2)*n3 + i3];
  }
  const T& operator()(unsigned int i1, unsigned int i2, unsigned int i3) const {
    ASSERT2(0<=i1 && i1<n1);
    ASSERT2(0<=i2 && i2<n2);
    ASSERT2(0<=i3 && i3<n3);
    return data[(i1*n2+i2)*n3 + i3];
  }

  Tensor& operator=(const T&val){
    for(auto &i: data){
      i = val;
    };
    return *this;
  };
  
  T* begin() { return std::begin(data);};
  const T* begin() const { return std::begin(data);};
  T* end() { return std::end(data);};
  const T* end() const { return std::end(data);};
  
  std::tuple<unsigned int, unsigned int, unsigned int> shape() { return std::make_tuple(n1, n2, n3);};
  
  bool empty(){
    return n1*n2*n3 == 0;
  }
  
  /*!
   * Ensures that this Tensor does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() {
    data.ensureUnique();
  }
 
private:
  unsigned int n1, n2, n3;
  Array<T> data;
};

/**************************************************************************
 * Matrix routines
 **************************************************************************/
/// Explicit inversion of a 3x3 matrix \p a
///
/// The input \p small determines how small the determinant must be for
/// us to throw due to the matrix being singular (ill conditioned);
/// If small is less than zero then instead of throwing we return 1.
/// This is ugly but can be used to support some use cases.
template <typename T> int invert3x3(Matrix<T> &a, BoutReal small = 1.0e-15) {
  TRACE("invert3x3");

  // Calculate the first co-factors
  T A = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1);
  T B = a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2);
  T C = a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0);

  // Calculate the determinant
  T det = a(0, 0) * A + a(0, 1) * B + a(0, 2) * C;

  if (abs(det) < abs(small)) {
    if (small >=0 ){
      throw BoutException("Determinant of matrix < %e --> Poorly conditioned", small);
    } else {
      return 1;
    }      
  }

  // Calculate the rest of the co-factors
  T D = a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2);
  T E = a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0);
  T F = a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1);
  T G = a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1);
  T H = a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2);
  T I = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);

  // Now construct the output, overwrites input
  T detinv = 1.0 / det;

  a(0, 0) = A * detinv;
  a(0, 1) = D * detinv;
  a(0, 2) = G * detinv;
  a(1, 0) = B * detinv;
  a(1, 1) = E * detinv;
  a(1, 2) = H * detinv;
  a(2, 0) = C * detinv;
  a(2, 1) = F * detinv;
  a(2, 2) = I * detinv;

  return 0;
};

// Give signature here as not able to mark implementation below as DEPRECATED
template <class T>
DEPRECATED(T **matrix(int xsize, int ysize));

/*!
 * Create a 2D array of \p xsize by \p ysize 
 * This is allocated as two blocks of data so that
 * the values are in a contiguous array.
 * 
 * Note: This returns C-style pointers, and makes
 * no effort to manage memory. Prefer other methods
 * (like standard containers) over this if possible.
 * 
 * \deprecated
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

  if((m = new T*[xsize]) == nullptr)
    throw BoutException("Error: could not allocate memory:%d\n", xsize);
  
  if((m[0] = new T[xsize*ysize]) == nullptr)
    throw BoutException("Error: could not allocate memory\n");

  for(i=1;i<xsize;i++) {
    m[i] = m[i-1] + ysize;
  }
  return m;
}

template <class T>
DEPRECATED(void free_matrix(T **m));
/*!
 * Free a matrix, assumed to have been allocated using matrix()
 *
 * @param[in] m  The matrix to free
 * @deprecated
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


/// Allocate a 3D BoutReal array of size \p nrow x \p ncol \p ndep
///
/// \deprecated Prefer other methods like standard containers
DEPRECATED(BoutReal ***r3tensor(int nrow, int ncol, int ndep));

/// Free a 3D BoutReal array, assumed to have been created by r3tensor
///
/// \deprecated
DEPRECATED(void free_r3tensor(BoutReal ***m));

/// Allocate a 3D int array of size \p nrow x \p ncol \p ndep
///
/// \deprecated Prefer other methods like standard containers
DEPRECATED(int ***i3tensor(int nrow, int ncol, int ndep));

/// Free a 3D int array, assumed to have been created by i3tensor()
///
/// \deprecated
DEPRECATED(void free_i3tensor(int ***m));

/*!
 * Get Random number between 0 and 1
 */
inline BoutReal randomu() {
  return static_cast<BoutReal>(rand()) / static_cast<BoutReal>(RAND_MAX);
}

/*!
 * Calculate the square of a variable \p t
 * i.e. t * t
 */
template <typename T>
T SQ(const T &t){
  return t*t;
}

/*!
 * Round \p x to the nearest integer
 */
inline int ROUND(BoutReal x){
  return (x > 0.0) ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
}

/// Calculate the maximum of a list of values
/// using a > b operator
template <typename T>
T BOUTMAX(T a) {
  return a;
}
template <typename T, typename... Args>
T BOUTMAX(T a, T b, Args... args) {
  T c = BOUTMAX(b, args...);
  return c > a ? c : a;
}

/// Calculate the minimum of a list of values
/// using the a < b operator
template <typename T>
T BOUTMIN(T a) {
  return a;
}
template <typename T, typename... Args>
T BOUTMIN(T a, T b, Args... args) {
  T c = BOUTMIN(b, args...);
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
  return a < 0 ? -1 : +1;
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

#if CHECK > 0
/// Throw an exception if \p f is not finite
inline void checkData(BoutReal f) {
  if (!finite(f)) {
    throw BoutException("BoutReal: Operation on non-finite data");
  }
}
#else
/// Ignored with disabled CHECK; Throw an exception if \p f is not finite
inline void checkData(BoutReal UNUSED(f)){};
#endif

/*!
 * Allocate memory and copy string \p s
 */ 
char* copy_string(const char* s);

/*!
 * Convert a value to a string
 * by writing to a stringstream
 */
template <class T>
const string toString(const T& val) {
  std::stringstream ss;
  ss << val;
  return ss.str();
}

/*!
 * Convert a string to lower case
 */
const string lowercase(const string &str);

/*!
 * Convert to lower case, except inside quotes (" or ')
 */
const string lowercasequote(const string &str);

/*!
 * Convert a string to a BoutReal
 * Throws BoutException if can't be done
 */ 
BoutReal stringToReal(const std::string &s);

/*!
 * Convert a string to an int
 * 
 * Throws BoutException if can't be done
 */
int stringToInt(const std::string &s);

/*!
 * Split a string on a given delimiter
 *
 * @param[in] s     The string to split (not modified by call)
 * @param[in] delim The delimiter to split on (single char)
 * @param[in, out] elems  A list to which the pieces will be appended using push_back
 */
std::list<std::string> &strsplit(const std::string &s, char delim, std::list<std::string> &elems);

/*!
 * Split a string on a given delimiter
 * 
 * @param[in] s     The string to split (not modified by call)
 * @param[in] delim The delimiter to split on (single char)
 */
std::list<std::string> strsplit(const std::string &s, char delim);

/*!
 * Strips leading and trailing spaces from a string
 * 
 * @param[in] s   The string to trim (not modified)
 * @param[in] c   Collection of characters to remove
 */
std::string trim(const std::string &s, const std::string &c=" \t\r");

/*!
 * Strips leading spaces from a string
 * 
 * @param[in] s   The string to trim (not modified)
 * @param[in] c   Collection of characters to remove
 */
std::string trimLeft(const std::string &s, const std::string &c=" \t");

/*!
 * Strips leading spaces from a string
 * 
 * @param[in] s   The string to trim (not modified)
 * @param[in] c   Collection of characters to remove
 */
std::string trimRight(const std::string &s, const std::string &c=" \t\r");

/*!
 * Strips the comments from a string
 * 
 * @param[in] s   The string to trim (not modified)
 * @param[in] c   Collection of characters to remove
 */
std::string trimComments(const std::string &s, const std::string &c="#;");

/// the bout_vsnprintf macro:
/// The first argument is an char * buffer of length len.
/// It needs to have been allocated with new[], as it may be
/// reallocated.
/// len: the length of said buffer. May be changed, mussn't be const.
/// fmt: the const char * descriping the format.
/// note that fmt should be the first argument of the function of type
/// const char * and has to be directly followed by the variable arguments.
#define bout_vsnprintf(buf,len,fmt) {                   \
    va_list va;                                         \
    va_start(va, fmt);                                  \
    int _vsnprintflen=vsnprintf(buf,len,fmt,va);        \
    va_end(va);                                         \
    if ( _vsnprintflen + 1 > int(len)) {                \
      _vsnprintflen+=1;                                 \
      delete[] buf;                                     \
      buf = new char[_vsnprintflen];                    \
      len = _vsnprintflen;                              \
      va_start(va,fmt);                                 \
      vsnprintf(buf,len,fmt,va);                        \
      va_end(va);                                       \
    }                                                   \
  }

#endif // __UTILS_H__
