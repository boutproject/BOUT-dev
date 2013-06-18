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

#include <utils.hxx>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>

BoutReal *rvector(int size) {
  return (BoutReal*) malloc(sizeof(BoutReal)*size);
}

BoutReal *rvresize(BoutReal *v, int newsize) {
  return (BoutReal*) realloc(v, sizeof(BoutReal)*newsize);
}

void rvfree(BoutReal *r) {
  free(r);
}

int *ivector(int size) {
  return (int*) malloc(sizeof(int)*size);
}

int *ivresize(int *v, int newsize) {
  return (int*) realloc(v, sizeof(int)*newsize);
}

void ivfree(int *v) {
  free(v);
}

BoutReal **rmatrix(int xsize, int ysize) {
  long i;
  BoutReal **m;

  if((m = (BoutReal**) malloc(xsize*sizeof(BoutReal*))) == (BoutReal**) NULL) {
    printf("Error: could not allocate memory:%d\n", xsize);
    exit(1);
  }

  if((m[0] = (BoutReal*) malloc(xsize*ysize*sizeof(BoutReal))) == (BoutReal*) NULL) {
    printf("Error: could not allocate memory\n");
    exit(1);
  }
  for(i=1;i<xsize;i++) {
    m[i] = m[i-1] + ysize;
  }

  return(m);
}

int **imatrix(int xsize, int ysize) {
  long i;
  int **m;

  if((m = (int**) malloc(xsize*sizeof(int*))) == (int**) NULL) {
    printf("Error: could not allocate memory:%d\n", xsize);
    exit(1);
  }

  if((m[0] = (int*) malloc(xsize*ysize*sizeof(int))) == (int*) NULL) {
    printf("Error: could not allocate memory\n");
    exit(1);
  }
  for(i=1;i<xsize;i++) {
    m[i] = m[i-1] + ysize;
  }

  return(m);
}

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

// Need explicit instantiation of some types
template BoutReal **matrix<BoutReal>(int,int);
template dcomplex **matrix<dcomplex>(int,int);

void free_rmatrix(BoutReal **m) {
  free(m[0]);
  free(m);
}

void free_imatrix(int **m) {
  free(m[0]);
  free(m);
}

template <class T>
void free_matrix(T **m) {
  delete[] m[0];
  delete[] m;
}

template void free_matrix<BoutReal>(BoutReal**);
template void free_matrix<dcomplex>(dcomplex**);

BoutReal ***r3tensor(int nrow, int ncol, int ndep) {
  int i,j;
  BoutReal ***t;

  /* allocate pointers to pointers to rows */
  t=(BoutReal ***) malloc((size_t)(nrow*sizeof(BoutReal**)));

  /* allocate pointers to rows and set pointers to them */
  t[0]=(BoutReal **) malloc((size_t)(nrow*ncol*sizeof(BoutReal*)));

  /* allocate rows and set pointers to them */
  t[0][0]=(BoutReal *) malloc((size_t)(nrow*ncol*ndep*sizeof(BoutReal)));

  for(j=1;j!=ncol;j++) t[0][j]=t[0][j-1]+ndep;
  for(i=1;i!=nrow;i++) {
    t[i]=t[i-1]+ncol;
    t[i][0]=t[i-1][0]+ncol*ndep;
    for(j=1;j!=ncol;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_r3tensor(BoutReal ***m) {
  free(m[0][0]);
  free(m[0]);
  free(m);
}

dcomplex **cmatrix(int nrow, int ncol) {
  dcomplex **m;
  int i;

  m = new dcomplex*[nrow];
  m[0] = new dcomplex[nrow*ncol];
  for(i=1;i<nrow;i++)
    m[i] = m[i-1] + ncol;

  return m;
}

void free_cmatrix(dcomplex** m) {
  delete[] m[0];
  delete[] m;
}

BoutReal randomu() {
  return ((BoutReal) rand()) / ((BoutReal) RAND_MAX);
}

BoutReal SQ(BoutReal x) {
  return(x*x);
}

int ROUND(BoutReal x) {
  return (x > 0.0) ? (int) (x + 0.5) : (int) (x - 0.5);
}

int BOUTMAX(int a, int b) {
  return (a > b) ? a : b;
}

BoutReal BOUTMAX(BoutReal a, BoutReal b) {
  return (a > b) ? a : b;
}

BoutReal BOUTMAX(BoutReal a, BoutReal b, BoutReal c) {
  return BOUTMAX(BOUTMAX(a,b), c);
}

BoutReal BOUTMIN(BoutReal a, BoutReal b) {
  return (a < b) ? a : b;
}
BoutReal BOUTMIN(BoutReal a, BoutReal b, BoutReal c) {
  return BOUTMIN(BOUTMIN(a,b),c);
}

bool is_pow2(int x) {
  return x && !((x-1) & x);
}

BoutReal SIGN(BoutReal a) {
  return (a >= 0) ? 1.0 : -1.0;
}

BoutReal MINMOD(BoutReal a, BoutReal b) {
  return 0.5*(SIGN(a) + SIGN(b)) * BOUTMIN(fabs(a), fabs(b));
}

/*
// integer power
BoutReal operator^(BoutReal lhs, int n)
{
  BoutReal result;
  
  if(n == 0)
    return 1.0;
  
  if(n < 0) {
    lhs = 1.0 / lhs;
    n *= -1;
  }
  
  result = 1.0;
  
  while(n > 1) {
    if( (n & 1) == 0 ) {
      lhs *= lhs;
      n /= 2;
    }else {
      result *= lhs;
      n--;
    }
  }
  result *= lhs;

  return result;
}

BoutReal operator^(BoutReal lhs, const BoutReal &rhs)
{
  return pow(lhs, rhs);
}
*/

/**************************************************************************
 * String routines
 **************************************************************************/

/// Allocate memory for a copy of given string
char* copy_string(const char* s) {
  char *s2;
  int n;

  if(s == NULL)
    return NULL;

  n = strlen(s);
  s2 = (char*) malloc(n+1);
  strcpy(s2, s);
  return s2;
}

/// Concatenate a string. This is a mildly nasty hack, and not thread-safe.
/// Simplifies some code though.
char *strconcat(const char* left, const char *right) {
  static char buffer[128];

  snprintf(buffer, 128, "%s%s", left, right);
  return buffer;
}

// Convert a value to a string
template <class T>
const string toString(const T& val) {
  std::stringstream ss;
  ss << val;
  return ss.str();
}

template const string toString(const int& val);
template const string toString(const BoutReal& val);

/// Convert a string to lower case
const string lowercase(const string &str) {
  string strlow(str);

  std::transform(strlow.begin(), strlow.end(), strlow.begin(), ::tolower);
  return strlow;
}

/// Convert to lowercase, except for inside strings
const string lowercasequote(const string &str) {
  string strlow(str);

  bool quote = false, dquote = false;
  for(int i=0;i<strlow.length(); i++) {
    if(strlow[i] == '\'') {
      quote ^= true;
    }else if(strlow[i] == '"') {
      dquote ^= true;
    }else if( (!quote) && (!dquote) ){
      strlow[i] = tolower(strlow[i]);
    }
  }
  return strlow;
}

BoutReal stringToReal(const std::string &s) {
  std::stringstream ss(s);
  BoutReal val;
  ss >> val;
  return val;
}

int stringToInt(const std::string &s) {
  std::stringstream ss(s);
  int val;
  ss >> val;
  return val;
}

std::list<std::string> &strsplit(const std::string &s, char delim, std::list<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::list<std::string> strsplit(const std::string &s, char delim) {
    std::list<std::string> elems;
    return strsplit(s, delim, elems);
}

// Strips leading and trailing spaces from a string
std::string trim(const std::string &s, const std::string &c) {
  return trimLeft(trimRight(s, c), c);
}

std::string trimRight(const std::string &s, const std::string &c) {
  std::string str(s);
  return str.erase(s.find_last_not_of(c)+1);
}

std::string trimLeft(const std::string &s, const std::string &c) {
  std::string str(s);
  return str.erase(0, s.find_first_not_of(c));
}

// Strips the comments from a string
// This is the compliment of trimLeft
std::string trimComments(const std::string &s, const std::string &c) {
  std::string str(s);
  return str.substr(0, s.find_first_of(c));
}

