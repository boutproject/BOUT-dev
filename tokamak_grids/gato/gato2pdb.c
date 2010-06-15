/***********************************************************************************
 * Converts a GATO file into a PDB file which is (almost) ready to be input to BOUT
 * B.Dudson, Nov 2007
 ***********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "read_gato.h"

#include "pdb.h"

char* DEFAULT_INPUT = "dskgato";

#define WRITE1D(name, var) \
  if(!PD_write_alt(fp, name, "double", var, 1, ind)) {                     \
    fprintf(stderr, "Error: Could not write to file '%s'\n", outfile);     \
    return(3);  \
  }

/* NOTE: Need to reverse the indices, and flip theta direction */
#define WRITE2D(name, var)                                           \
  for(i=0;i<data.npsi;i++)                                           \
    for(j=0;j<data.ntheta;j++)                                       \
      m[i][j] = var[(data.ntheta-1) - j][i];                         \
  if(!PD_write_alt(fp, name, "double", *m, 2, ind)) {                \
    fprintf(stderr, "Error: Could not write to file %s\n", outfile); \
  }

#define SQ(x) ((x)*(x))

int main(int argc, char **argv)
{
  char *infile;
  char *outfile;
  gato_data data;
  PDBfile *fp;
  long ind[6];
  double **m;
  int i, j, ip, im;
  double dl, dr, dz;
  
  /* Set input file name */
  infile = DEFAULT_INPUT;
  if(argc > 1) {
    /* First argument is the input file */
    infile = argv[1];
  }
  
  /* Set output file name */
  if(argc > 2) {
    outfile = argv[2];
  }else {
    /* Default output is input + ".pdb" */
    outfile = (char*) malloc(strlen(infile) + 5);
    sprintf(outfile, "%s.pdb", infile);
  }

  /* Read the GATO data file */
  if(read_gato(infile, &data)) {
    return(1);
  }

  if(!data.b_present) {
    /* Calculate Br and Bz from Bpol */
    data.Br = dmatrix(data.ntheta, data.npsi);
    data.Bz = dmatrix(data.ntheta, data.npsi);
    data.b_present = 1;

    for(i=0;i<data.ntheta;i++) {
      ip = (i == (data.ntheta-1)) ? 0 : i+1;
      im = (i == 0) ? data.ntheta-1 : i-1;

      for(j=0;j<data.npsi;j++) {
	dl = sqrt(SQ(data.xnorm[ip][j] - data.xnorm[im][j]) + SQ(data.znorm[ip][j] - data.znorm[im][j]));
	dr = data.xnorm[ip][j] - data.xnorm[im][j];
	dz = data.znorm[ip][j] - data.znorm[im][j];

	data.Br[i][j] = dr * data.Bpol[i][j] / dl;
	data.Bz[i][j] = dz * data.Bpol[i][j] / dl;
      }
    }
  }
  
  /* Open PDB file */
  
  if((fp = PD_open(outfile, "w")) == (PDBfile*) NULL) {
    fprintf(stderr, "Error: Could not open '%s' for writing\n", outfile);
    return(2);
  }
  
  /* Write to PDB file */
  
  ind[0] = 0L; ind[1] = data.npsi-1; ind[2] = 1L;
  ind[3] = 0L; ind[4] = data.ntheta-1; ind[5] = 1L;

  /* Size of the grid */
  if(!PD_write(fp, "nx", "integer", &data.npsi)) {
    fprintf(stderr, "Error: Could not write to file '%s'\n", outfile);
    return(3);
  }
  if(!PD_write(fp, "ny", "integer", &data.ntheta)) {
    fprintf(stderr, "Error: Could not write to file '%s'\n", outfile);
    return(3);
  }

  /* Write 1D flux surface functions */
  WRITE1D("psi",       data.psiflux);
  WRITE1D("f",         data.fnorm);
  WRITE1D("ffprime",   data.ffpnorm);
  WRITE1D("mu0p",      data.ponly);
  WRITE1D("mu0pprime", data.pponly);
  WRITE1D("qsafe",     data.qsf);
  WRITE1D("Ni",        data.d);
  
  /* Write 2D variables (Need to reverse indices) */
  m = dmatrix(data.npsi, data.ntheta);
  
  WRITE2D("Rxy", data.xnorm);
  WRITE2D("Zxy", data.znorm);

  WRITE2D("Brxy", data.Br);
  WRITE2D("Bzxy", data.Bz);
  WRITE2D("Bpxy", data.Bpol);
  WRITE2D("Btxy", data.Btor);

  PD_close(fp);

  free_gato(&data); /* Free data */

  return(0);
}
