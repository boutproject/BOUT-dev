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
#include "msg_stack.hxx"
#include "unused.hxx"

#include <string>
#include <list>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <memory>

namespace bout {
namespace utils {
#ifndef __cpp_lib_make_unique
// Provide our own make_unique if the stl doesn't give us one
// Implementation from https://isocpp.org/files/papers/N3656.txt
// i.e. what's already in the stl
template <class T>
struct _Unique_if {
  using _Single_object = std::unique_ptr<T>;
};

template <class T>
struct _Unique_if<T[]> {
  using _Unknown_bound = std::unique_ptr<T[]>;
};

template <class T, size_t N>
struct _Unique_if<T[N]> {
  using _Known_bound = void;
};

template <class T, class... Args>
typename _Unique_if<T>::_Single_object make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <class T>
typename _Unique_if<T>::_Unknown_bound make_unique(size_t n) {
  using U = typename std::remove_extent<T>::type;
  return std::unique_ptr<T>(new U[n]());
}

template <class T, class... Args>
typename _Unique_if<T>::_Known_bound make_unique(Args&&...) = delete;
#else
using std::make_unique;
#endif

template <typename T>
struct function_traits;

/// Traits class to get the types of function arguments for function pointers
///
/// Use like:
///
//      // A function signature we'd like to check:
///     using some_function = int(*)(int, double, std::string);
///     // Get the type of the first argument:
///     using first_argument_type =
///         bout::utils::function_traits<some_function>::arg<1>::type;
///     // The following prints "true":
///     std::cout << std::boolalpha
///         << std::is_same<double, first_argument_type>::value;
///
/// Adapted from https://stackoverflow.com/a/9065203/2043465
template <typename R, typename... Args>
struct function_traits<R (*)(Args...)> {
  /// Total number of arguments
  static constexpr size_t nargs = sizeof...(Args);

  using result_type = R;

  template <size_t i>
  struct arg {
    using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
  };
};
} // namespace utils
} // namespace bout

/// Helper class for 2D arrays
///
/// Allows bounds checking through `operator()` with CHECK > 1
///
/// If either \p n1 or \p n2 are 0, the Matrix is empty and should not
/// be indexed
template <typename T>
class Matrix {
public:
  using data_type = T;
  using size_type = int;
  
  Matrix() : n1(0), n2(0){};
  Matrix(size_type n1, size_type n2) : n1(n1), n2(n2) {
    ASSERT2(n1 >= 0);
    ASSERT2(n2 >= 0);

    data.reallocate(n1 * n2);
  }
  Matrix(const Matrix &other) : n1(other.n1), n2(other.n2), data(other.data) {
    // Prevent copy on write for Matrix
    data.ensureUnique();
  }

  /// Reallocate the Matrix to shape \p new_size_1 by \p new_size_2
  ///
  /// Note that this invalidates the existing data!
  void reallocate(size_type new_size_1, size_type new_size_2) {
    ASSERT2(new_size_1 >= 0);
    ASSERT2(new_size_2 >= 0);

    n1 = new_size_1;
    n2 = new_size_2;
    data.reallocate(new_size_1 * new_size_2);
  }

  Matrix& operator=(const Matrix &other) {
    n1 = other.n1;
    n2 = other.n2;
    data = other.data;
    // Prevent copy on write for Matrix
    data.ensureUnique();
    return *this;
  }
  
  inline T& operator()(size_type i1, size_type i2) {
    ASSERT2(0<=i1 && i1<n1);
    ASSERT2(0<=i2 && i2<n2);
    return data[i1*n2+i2];
  }
  inline const T& operator()(size_type i1, size_type i2) const {
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

  T* begin() { return std::begin(data);};
  const T* begin() const { return std::begin(data);};
  T* end() { return std::end(data);};
  const T* end() const { return std::end(data);};

  std::tuple<size_type, size_type> shape() const { return std::make_tuple(n1, n2); };

  bool empty() const { return n1 * n2 == 0; }

  /*!
   * Ensures that this Matrix does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() {
    data.ensureUnique();
  }
  
  /// Access the underlying storage
  Array<T>& getData() { return data; }
  const Array<T>& getData() const { return data; }

private:
  size_type n1, n2;
  /// Underlying 1D storage array
  Array<T> data;
};

/// Helper class for 3D arrays
///
/// Allows bounds checking through `operator()` with CHECK > 1
///
/// If any of \p n1, \p n2 or \p n3 are 0, the Tensor is empty and
/// should not be indexed
template <typename T>
class Tensor {
public:
  using data_type = T;
  using size_type = int;

  Tensor() : n1(0), n2(0), n3(0) {};
  Tensor(size_type n1, size_type n2, size_type n3) : n1(n1), n2(n2), n3(n3) {
    ASSERT2(n1 >= 0);
    ASSERT2(n2 >= 0);
    ASSERT2(n3 >= 0);
    data.reallocate(n1 * n2 * n3);
  }
  Tensor(const Tensor &other) : n1(other.n1), n2(other.n2), n3(other.n3), data(other.data) {
    // Prevent copy on write for Tensor
    data.ensureUnique();
  }

  /// Reallocate the Tensor with shape \p new_size_1 by \p new_size_2 by \p new_size_3
  ///
  /// Note that this invalidates the existing data!
  void reallocate(size_type new_size_1, size_type new_size_2, size_type new_size_3) {
    ASSERT2(new_size_1 >= 0);
    ASSERT2(new_size_2 >= 0);
    ASSERT2(new_size_3 >= 0);

    n1 = new_size_1;
    n2 = new_size_2;
    n3 = new_size_3;
    data.reallocate(new_size_1 * new_size_2 * new_size_3);
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

  T& operator()(size_type i1, size_type i2, size_type i3) {
    ASSERT2(0<=i1 && i1<n1);
    ASSERT2(0<=i2 && i2<n2);
    ASSERT2(0<=i3 && i3<n3);
    return data[(i1*n2+i2)*n3 + i3];
  }
  const T& operator()(size_type i1, size_type i2, size_type i3) const {
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

  std::tuple<size_type, size_type, size_type> shape() const {
    return std::make_tuple(n1, n2, n3);
  };

  bool empty() const { return n1 * n2 * n3 == 0; }

  /*!
   * Ensures that this Tensor does not share data with another
   * This should be called before performing any write operations
   * on the data.
   */
  void ensureUnique() {
    data.ensureUnique();
  }

  /// Access the underlying storage
  Array<T>& getData() { return data; }
  const Array<T>& getData() const { return data; }

private:
  size_type n1, n2, n3;
  /// Underlying 1D storage array
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

  if (std::abs(det) < std::abs(small)) {
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
  return 0.5*(SIGN(a) + SIGN(b)) * BOUTMIN(std::abs(a), std::abs(b));
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


/// Convert a value to a string
/// by writing to a stringstream
template <class T>
std::string toString(const T& val) {
  std::stringstream ss;
  ss << val;
  return ss.str();
}

/// Simple case where input is already a string
/// This is so that toString can be used in templates
/// where the type may be std::string.
inline std::string toString(const std::string& val) {
  return val;
}

template <>
inline std::string toString<>(const Array<BoutReal>& UNUSED(val)) {
  return "<Array>";
}

template <>
inline std::string toString<>(const Matrix<BoutReal>& UNUSED(val)) {
  return "<Matrix>";
}

template <>
inline std::string toString<>(const Tensor<BoutReal>& UNUSED(val)) {
  return "<Tensor>";
}

/// Convert a bool to "true" or "false"
inline std::string toString(const bool& val) {
  if (val) {
    return "true";
  }
  return "false";
}

/// Convert a time stamp to a string
/// This uses std::localtime and std::put_time
std::string toString(const time_t& time);

/*!
 * Convert a string to lower case
 */
const std::string lowercase(const std::string &str);

/*!
 * Convert a string to upper case
 */
const std::string uppercase(const std::string &str);

/*!
 * Convert to lower case, except inside quotes (" or ')
 */
const std::string lowercasequote(const std::string &str);

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

/// Convert pointer or reference to pointer
/// This allows consistent handling of both in macros, templates
template <typename T> T *pointer(T *val) { return val; }
template <typename T> T *pointer(T &val) { return &val; }

#endif // __UTILS_H__
