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

void free_rmatrix(BoutReal **m) {
  free(m[0]);
  free(m);
}

void free_imatrix(int **m) {
  free(m[0]);
  free(m);
}

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

int ***i3tensor(int nrow, int ncol, int ndep) {
  int i,j;
  int ***t;

  /* allocate pointers to pointers to rows */
  t=(int ***) malloc((size_t)(nrow*sizeof(int**)));

  /* allocate pointers to rows and set pointers to them */
  t[0]=(int **) malloc((size_t)(nrow*ncol*sizeof(int*)));

  /* allocate rows and set pointers to them */
  t[0][0]=(int *) malloc((size_t)(nrow*ncol*ndep*sizeof(int)));

  for(j=1;j!=ncol;j++) t[0][j]=t[0][j-1]+ndep;
  for(i=1;i!=nrow;i++) {
    t[i]=t[i-1]+ncol;
    t[i][0]=t[i-1][0]+ncol*ndep;
    for(j=1;j!=ncol;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_i3tensor(int ***m) {
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
