/***********************************************************************************
 * Reads a GATO format grid file
 * Adapted from ELITE FORTRAN code by B.Dudson, Nov 2007
 ***********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "read_gato.h"

#define BUF_LEN 256

#define READ_NEXT() \
  if(fgets(buffer, BUF_LEN-1, fp) == (char*) NULL) {                       \
    fprintf(stderr, "Error: Could not read line from '%s'\n", filename);   \
    return(2);                                                             \
  }

#define READ_ARRAY(arr, n, name) \
  for(i=0;i<(n);i++) {                         \
    if(fscanf(fp, "%le", &(arr[i])) != 1) {    \
      fprintf(stderr, "Error: Could not read %s[%d] from '%s'\n", name, i, filename); \
      return(2); \
    } \
  }

#define SQ(x) ((x)*(x))

#define PI 3.1415926
#define MU (4.0*PI*1.0e-7)

double* dvector(int n)
{
  double *a;
  a = (double*) malloc(sizeof(double)*n);
  if(a == (double*) NULL) {
    fprintf(stderr, "Error: Could not allocate memory\n");
    exit(1);
  }
  return(a);
}

double** dmatrix(int nx, int ny)
{
  double **m;
  int i;

  m = (double**) malloc(sizeof(double*)*nx);
  if(m == (double**) NULL) {
    fprintf(stderr, "Error: Could not allocate memory\n");
    exit(1);
  }
  m[0] = (double*) malloc(sizeof(double)*nx*ny);
  if(m[0] == (double*) NULL) {
    fprintf(stderr, "Error: Could not allocate memory\n");
    exit(1);
  }
  for(i=1;i<nx;i++)
    m[i] = m[i-1] + ny;

  return(m);
}

void free_dmatrix(double **m)
{
  free(m[0]);
  free(m);
}

int read_gato(char *filename, gato_data* data)
{
  FILE *fp;
  char buffer[BUF_LEN];
  int i, j;
  int tsize, t0;
  double dx, dpsi;

  if((fp = fopen(filename, "rt")) == (FILE*) NULL) {
    fprintf(stderr, "Error: Could not open file '%s'\n", filename);
    return(1);
  }

  /* First line is the date */

  if(fgets(buffer, BUF_LEN-1, fp) == (char*) NULL) {
    fprintf(stderr, "Error: Could not read date from '%s'\n", filename);
    return(2);
  }
  printf("Date: %s", buffer);
  
  /* Second line is a description */

  if(fgets(buffer, BUF_LEN-1, fp) == (char*) NULL) {
    fprintf(stderr, "Error: Could not read description from '%s'\n", filename);
    return(2);
  }
  printf("Description: %s", buffer);
  
  /* Third line contains number of mesh points */

  READ_NEXT();

  if(sscanf(buffer, "%5d %5d %5d", &(data->npsi), &(data->ntheta), &(data->isym)) != 3) {
    fprintf(stderr, "Error: Could not parse number of mesh points\n");
    return(3);
  }
  printf("Grid size: %d, %d, %d\n", data->npsi, data->ntheta, data->isym);

  /* Check values */
  if((data->isym < 0) || (data->isym > 1)) {
    fprintf(stderr, "Error: isym must be either 0 or 1\n");
    return(4);
  }
  if((data->npsi < 1) || (data->ntheta < 1)) {
    fprintf(stderr, "Error: npsi and ntheta must be > 0\n");
    return(4);
  }

  if(data->isym) {
    tsize = 2*(data->ntheta - 1) + 1;
  }else {
    tsize = data->ntheta;
  }

  /* Read normalisation factors */
  READ_NEXT();
  if(sscanf(buffer, "%le %le %le %le", 
	    &(data->rcnt), &(data->xma), &(data->zma), &(data->btor)) != 4) {
    fprintf(stderr, "Error: Could not parse normalisation factors\n");
    return(3);
  }
  printf("rcnt = %e, xma = %e, zma = %e, btor = %e\n", 
	 data->rcnt, data->xma, data->zma, data->btor);

  READ_NEXT();
  if(sscanf(buffer, "%le %le %le", &(data->curtot), &(data->eaxe), &(data->dnorm)) != 3) {
    fprintf(stderr, "Error: Could not parse normalisation factors\n");
    return(3);
  }
  printf("curtot = %e, eaxe = %e, dnorm = %e\n", data->curtot, data->eaxe, data->dnorm);

  /* Allocate memory */
  data->psiflux = dvector(data->npsi);
  data->fnorm   = dvector(data->npsi);
  data->ffpnorm = dvector(data->npsi);
  data->ponly   = dvector(data->npsi);
  data->pponly  = dvector(data->npsi);
  data->qsf     = dvector(data->npsi);
  data->d       = dvector(data->npsi);
  
  data->dpdz = dvector(tsize);
  data->dpdr = dvector(tsize);

  data->xnorm = dmatrix(tsize, data->npsi);
  data->znorm = dmatrix(tsize, data->npsi);

  /* Read 1D arrays */
  READ_ARRAY(data->psiflux, data->npsi, "psiflux");
  READ_ARRAY(data->fnorm,   data->npsi, "fnorm");
  READ_ARRAY(data->ffpnorm, data->npsi, "ffpnorm");
  READ_ARRAY(data->ponly,   data->npsi, "ponly");
  READ_ARRAY(data->pponly,  data->npsi, "pponly");
  READ_ARRAY(data->qsf,     data->npsi, "qsf");
  READ_ARRAY(data->d,       data->npsi, "dnorm");

  READ_ARRAY(data->dpdz,    data->ntheta, "dpdz");
  READ_ARRAY(data->dpdr,    data->ntheta, "dpdr");

  /* Read 2D arrays */
  for(j=0;j<data->ntheta;j++) {
    sprintf(buffer, "xnorm[%d]", j);
    READ_ARRAY(data->xnorm[j], data->npsi, buffer);
  }
  for(j=0;j<data->ntheta;j++) {
    sprintf(buffer, "znorm[%d]", j);
    READ_ARRAY(data->znorm[j], data->npsi, buffer);
  }

  /* Try to read Br and Bz (may be present) */

  data->Br = dmatrix(tsize, data->npsi);
  data->b_present = 1;
  for(i=0;i<data->ntheta;i++) {
    for(j=0;j<data->npsi;j++) {
      if(fscanf(fp, "%le", &(data->Br[i][j])) != 1) {
	printf("File does not contain Br and Bz\n");
	data->b_present = 0;
	free_dmatrix(data->Br);
	i = data->ntheta;
	j = data->npsi;
      } 
    }
  }
  
  if(data->b_present) {
    data->Bz = dmatrix(tsize, data->npsi);
    for(j=0;j<data->ntheta;j++) {
      sprintf(buffer, "Bz[%d]", j);
      READ_ARRAY(data->Bz[j], data->npsi, buffer);
    }
    printf("File contains Br and Bz\n");
  }

  /* Check psi convention */

  if(data->psiflux[0] > data->psiflux[data->npsi-1]) {
    
    if(data->btor > 0.0) {
      i = 1;
    }else if(data->btor < 0.0) {
      i = -1;
    }else {
      fprintf(stderr, "WARNING: Btor = 0 in '%s'\n", filename);
      i = 0;
    }

    if(i != 0) {
      if(i > 0) {
	printf("Reversing sign of curtot\n");
	data->curtot *= -1.0;
      }else {
	printf("Reversing sign of btor");
	data->btor *= -1.0;
      }

      for(j=0;j<data->npsi;j++) {
	data->psiflux[j] *= -1.0;
	data->ffpnorm[j] *= -1.0;
	data->pponly[j]  *= -1.0;
	if((i < 0) && (data->fnorm[0] < 0.0))
	  data->fnorm[j] *= -1.0;
      }
    }
  }
  
  /* Fill in values for up-down symmetric case */

  if(data->isym) {
    printf("Grid is up-down symmetric. Reflecting grid about midplane\n");

    for(i=data->ntheta;i<tsize;i++) {
      t0 = tsize - 1 - i;

      for(j=0;j<data->npsi;j++) {
	data->xnorm[i][j] = data->xnorm[t0][j];
	data->znorm[i][j] = 2.0*data->zma - data->znorm[t0][j]; /* Reflect about zma */
	
	if(data->b_present) {
	  data->Br[i][j] = -1.0*data->Br[t0][j]; /* Br reverses */
	  data->Bz[i][j] = data->Bz[t0][j];      /* Bz is the same */
	}
      }

      data->dpdz[i] = -1.0*data->dpdz[t0]; /* dp/dz reverses */
      data->dpdr[i] = data->dpdr[t0];      /* dp/dr is the same */
    }

    data->ntheta = tsize; /* Entire grid now stored */
  }

  /* Calculate poloidal field */

  data->Bpol = dmatrix(data->ntheta, data->npsi);
  if(data->b_present) {
    for(i=0;i<data->ntheta;i++)
      for(j=0;j<data->npsi;j++)
	data->Bpol[i][j] = sqrt(SQ(data->Br[i][j]) + SQ(data->Bz[i][j]));
  }else {
    printf("Calculating poloidal field from psi\n");
    /* Use dpsi = 2*PI*R*Bp dx (for now) */
    
    for(i=0;i<data->ntheta;i++) {
      for(j=0;j<data->npsi;j++) {
	if(j == 0) {
	  /* Inner edge */
	  dx = sqrt(SQ(data->xnorm[i][1] - data->xnorm[i][0]) + SQ(data->znorm[i][1] - data->znorm[i][0]));
	  dpsi = data->psiflux[1] - data->psiflux[0];
	}else if(j == (data->npsi - 1)) {
	  dx = sqrt(SQ(data->xnorm[i][j] - data->xnorm[i][j-1]) + SQ(data->znorm[i][j] - data->znorm[i][j-1]));
	  dpsi = data->psiflux[j] - data->psiflux[j-1];
	}else {
	  dx = sqrt(SQ(data->xnorm[i][j+1] - data->xnorm[i][j-1]) + SQ(data->znorm[i][j+1] - data->znorm[i][j-1]));
	  dpsi = data->psiflux[j+1] - data->psiflux[j-1];
	}
	
	data->Bpol[i][j] = dpsi / (dx * data->xnorm[i][j]);
      }
    }
  }

  /* Calculate toroidal magnetic field */

  data->Btor = dmatrix(data->ntheta, data->npsi);
  for(i=0;i<data->ntheta;i++) {
    for(j=0;j<data->npsi;j++) {
      data->Btor[i][j] = data->fnorm[j]/data->xnorm[i][j];
    }
  }

  fclose(fp);
  return(0);
}

/* Free gato data structure */
void free_gato(gato_data *data)
{
  free(data->psiflux);
  free(data->fnorm);
  free(data->ffpnorm);
  free(data->ponly);
  free(data->pponly);
  free(data->qsf);
  free(data->d);

  free(data->dpdz);
  free(data->dpdr);
  
  free_dmatrix(data->xnorm);
  free_dmatrix(data->znorm);

  if(data->b_present) {
    free_dmatrix(data->Br);
    free_dmatrix(data->Bz);
  }

  free_dmatrix(data->Bpol);
  free_dmatrix(data->Btor);
}

/*
int main()
{
  gato_data data;

  read_gato("krbm1.dskgato", &data);
  free_gato(&data);
  
  return(0);
}
*/
