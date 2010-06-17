
#include "utils.h"

#include <stdlib.h>
#include <stdio.h>

float *fvector(int size)
{
  return (float*) malloc(sizeof(float)*size);
}

void free_fvector(float *m)
{
  free(m);
}

float **fmatrix(int xsize, int ysize)
{
  long i;
  float **m;
  
  if((m = (float**) malloc(xsize*sizeof(float*))) == (float**) NULL) {
    printf("Error: could not allocate memory:%d\n", xsize);
    exit(1);
  }

  if((m[0] = (float*) malloc(xsize*ysize*sizeof(float))) == (float*) NULL) {
    printf("Error: could not allocate memory\n");
    exit(1);
  }
  for(i=1;i!=xsize;i++) {
    m[i] = m[i-1] + ysize;
  }

  return(m);
}

void free_fmatrix(float **m)
{
  free(m[0]);
  free(m);
}

float ***f3tensor(int nrow, int ncol, int ndep)
{
  int i,j;
  float ***t;

  /* allocate pointers to pointers to rows */
  t=(float ***) malloc((size_t)(nrow*sizeof(float**)));

  /* allocate pointers to rows and set pointers to them */
  t[0]=(float **) malloc((size_t)(nrow*ncol*sizeof(float*)));

  /* allocate rows and set pointers to them */
  t[0][0]=(float *) malloc((size_t)(nrow*ncol*ndep*sizeof(float)));

  for(j=1;j!=ncol;j++) t[0][j]=t[0][j-1]+ndep;
  for(i=1;i!=nrow;i++) {
    t[i]=t[i-1]+ncol;
    t[i][0]=t[i-1][0]+ncol*ndep;
    for(j=1;j!=ncol;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_f3tensor(float ***m)
{
  free(m[0][0]);
  free(m[0]);
  free(m);
}

float ****f4tensor(int nrow, int ncol, int ndep, int ne)
{
  int i;
  float ****t;

  t         = (float ****) malloc((size_t)(nrow*sizeof(float***)));
  t[0]      = (float ***)  malloc((size_t)(nrow*ncol*sizeof(float**)));
  t[0][0]   = (float **)   malloc((size_t)(nrow*ncol*ndep*sizeof(float*)));
  t[0][0][0]= (float *)    malloc((size_t)(nrow*ncol*ndep*ne*sizeof(float*)));
  
  for(i=1;i<nrow;i++)
    t[i] = t[i-1] + ncol;
  
  for(i=1;i<nrow*ncol;i++)
    t[0][i] = t[0][i-1] + ndep;
  
  for(i=1;i<nrow*ncol*ndep;i++)
    t[0][0][i] = t[0][0][i-1] + ne;
  
  return t;
}

void free_f4tensor(float ****m)
{
  free(m[0][0][0]);
  free(m[0][0]);
  free(m[0]);
  free(m);
}
