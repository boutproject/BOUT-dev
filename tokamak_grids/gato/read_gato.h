/***********************************************************************************
 * Reads a GATO format grid file
 * Adapted from ELITE FORTRAN code by B.Dudson, Nov 2007
 ***********************************************************************************/

#ifndef __READ_GATO_H__
#define __READ_GATO_H__

typedef struct {
  int npsi, ntheta, isym;  /* Number of mesh points */
  /* If isym == 1, up-down symmetric */

  double rcnt, xma, zma, btor;
  double curtot,eaxe,dnorm;
  
  double *psiflux; /* Poloidal flux */
  double *fnorm;   /* f */
  double *ffpnorm; /* ff' */
  double *ponly;   /* mu0 p */
  double *pponly;  /* mu0 p' */
  double *qsf;     /* Safety factor */
  double *d;       /* Density set = d */
  
  double *dpdz;    /* dpsi / dz on boundary */
  double *dpdr;    /* dpsi / dr on boundary */

  double **xnorm;  /* R */
  double **znorm;  /* Z */

  int b_present;   /* 1 if Br and Bz are present */
  double **Br, **Bz;

  double **Bpol, **Btor; /* Poloidal and toroidal magnetic fields */
}gato_data;

int read_gato(char *filename, gato_data* data);
void free_gato(gato_data *data);

double* dvector(int n);
double** dmatrix(int nx, int ny);
void free_dmatrix(double **m);

#endif /* __READ_GATO_H__ */
