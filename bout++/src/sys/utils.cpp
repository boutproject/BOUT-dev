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

#include "utils.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <algorithm>
#include <sstream>

BoutReal *rvector(int size)
{
  return (BoutReal*) malloc(sizeof(BoutReal)*size);
}

BoutReal *rvresize(BoutReal *v, int newsize)
{
  return (BoutReal*) realloc(v, sizeof(BoutReal)*newsize);
}

int *ivector(int size)
{
  return (int*) malloc(sizeof(int)*size);
}

int *ivresize(int *v, int newsize)
{
  return (int*) realloc(v, sizeof(int)*newsize);
}

BoutReal **rmatrix(int xsize, int ysize)
{
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
  for(i=1;i!=xsize;i++) {
    m[i] = m[i-1] + ysize;
  }

  return(m);
}

int **imatrix(int xsize, int ysize)
{
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
  for(i=1;i!=xsize;i++) {
    m[i] = m[i-1] + ysize;
  }

  return(m);
}

void free_rmatrix(BoutReal **m)
{
  free(m[0]);
  free(m);
}

void free_imatrix(int **m)
{
  free(m[0]);
  free(m);
}

BoutReal ***r3tensor(int nrow, int ncol, int ndep)
{
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

void free_r3tensor(BoutReal ***m)
{
  free(m[0][0]);
  free(m[0]);
  free(m);
}

dcomplex **cmatrix(int nrow, int ncol)
{
  dcomplex **m;
  int i;

  m = new dcomplex*[nrow];
  m[0] = new dcomplex[nrow*ncol];
  for(i=1;i<nrow;i++)
    m[i] = m[i-1] + ncol;

  return m;
}

void free_cmatrix(dcomplex** m)
{
  delete[] m[0];
  delete[] m;
}

BoutReal SQ(BoutReal x)
{
  return(x*x);
}

int ROUND(BoutReal x)
{
  return (x > 0.0) ? (int) (x + 0.5) : (int) (x - 0.5);
}

void SWAP(BoutReal &a, BoutReal &b)
{
  BoutReal tmp;
  
  tmp = a;
  a = b;
  b = tmp;
}

void SWAP(dcomplex &a, dcomplex &b)
{
  dcomplex tmp;
  
  tmp = a;
  a = b;
  b = tmp;
}

void SWAP(int &a, int &b)
{
  int tmp;
  tmp = a;
  a = b;
  b = tmp;
}

int MAX(int a, int b)
{
  return (a > b) ? a : b;
}

bool is_pow2(int x)
{
  return x && !((x-1) & x);
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

///.allocate memory for a copy of given string
char* copy_string(const char* s)
{
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
char *strconcat(const char* left, const char *right)
{
  static char buffer[128];
  
  snprintf(buffer, 128, "%s%s", left, right);
  return buffer;
}

/// Convert a string to lower case
const string lowercase(const string &str)
{
  string strlow(str);
  
  std::transform(strlow.begin(), strlow.end(), strlow.begin(), ::tolower);
  return strlow;
}

/// Convert to lowercase, except for inside strings
const string lowercasequote(const string &str)
{
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

std::list<std::string> &strsplit(const std::string &s, char delim, std::list<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::list<std::string> strsplit(const std::string &s, char delim)
{
    std::list<std::string> elems;
    return strsplit(s, delim, elems);
}
