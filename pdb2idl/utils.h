/***************************************************************
 * Miscellaneous utilities for PDB2IDL
 ***************************************************************/

#ifndef __UTILS_H__
#define __UTILS_H__

float *fvector(int size);
void free_fvector(float *m);

float **fmatrix(int xsize, int ysize);
void free_fmatrix(float **m);

float ***f3tensor(int nrow, int ncol, int ndep);
void free_f3tensor(float ***m);

float ****f4tensor(int nrow, int ncol, int ndep, int n);
void free_f4tensor(float ****m);


#endif // __UTILS_H__
