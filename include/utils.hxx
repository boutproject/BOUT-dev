/**************************************************************************
 * Memory allocation routines
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

#include <string>
#include <list>
#include <cmath>
#include <algorithm>

using std::abs;
using std::swap;

BoutReal *rvector(int size);
BoutReal *rvresize(BoutReal *v, int newsize);
void rvfree(BoutReal *r);

int *ivector(int size);
int *ivresize(int *v, int newsize);
void ivfree(int *v);

BoutReal **rmatrix(int xsize, int ysize);
int **imatrix(int xsize, int ysize);

template <class T>
T **matrix(int xsize, int ysize);

void free_rmatrix(BoutReal **m);
void free_imatrix(int **m);

template <class T>
void free_matrix(T **m);

BoutReal ***r3tensor(int nrow, int ncol, int ndep);
void free_r3tensor(BoutReal ***m);

int ***i3tensor(int nrow, int ncol, int ndep);
void free_i3tensor(int ***m);

dcomplex **cmatrix(int nrow, int ncol);
void free_cmatrix(dcomplex** cm);

BoutReal randomu(); // Get Random number between 0 and 1

template <typename T>
T SQ(T t){
  return t*t;
}

inline int ROUND(BoutReal x){
  return (x > 0.0) ? (int) (x + 0.5) : (int) (x - 0.5);
}

template <typename T>
T BOUTMAX(T a){
  return a;
}
template <typename T, typename... Args>
T BOUTMAX(T a,T b,Args... args){
  T c = BOUTMAX(b,args...);
  return c > a ? c : a;
}

template <typename T>
T BOUTMIN(T a){
  return a;
}
template <typename T, typename... Args>
T BOUTMIN(T a,T b,Args... args){
  T c = BOUTMIN(b,args...);
  return c < a ? c : a;
}

inline bool is_pow2(int x) { // Check if a number is a power of 2
  return x && !((x-1) & x);
}

template <typename T>
T SIGN(T a) { // Return +1 or -1 (0 -> +1)
  return a > 0 ? 1 : -1;
}
BoutReal MINMOD(BoutReal a, BoutReal b);

char* copy_string(const char* s);

// Convert a value to a string
template <class T>
const string toString(const T& val);

/// Convert a string to lower case
const string lowercase(const string &str);
const string lowercasequote(const string &str);

// Convert a string to a BoutReal
BoutReal stringToReal(const std::string &s);
int stringToInt(const std::string &s);

/// Split a string on a given delimiter
std::list<std::string> &strsplit(const std::string &s, char delim, std::list<std::string> &elems);
std::list<std::string> strsplit(const std::string &s, char delim);

string trim(const string &, const string &c=" \t\r");
string trimLeft(const string &, const string &c=" \t");
string trimRight(const string &, const string &c=" \t\r");
string trimComments(const string &, const string &c="#;");

/// the myvsnprintf macro:
/// The first argument is an char * buffer of length len.
/// It needs to have been allocated with new[], as it may be
/// reallocated.
/// len: the length of said buffer. May be changed, mussn't be const.
/// fmt: the const char * descriping the format.
/// va: A dummy argument. will
/// note that fmt should be the first argument of the function of type
/// const char * and has to be directly followed by the variable arguments.
#define myvsnprintf(buf,len,fmt,va) {                     \
    va_list va;                                         \
    va_start(va, fmt);                                  \
    int _vsnprintflen=vsnprintf(buf,len,fmt,va);        \
    va_end(va);                                         \
    if ( _vsnprintflen+1 > len){                        \
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
