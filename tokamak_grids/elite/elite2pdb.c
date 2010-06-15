/*************************************************************
 * Converts a .eqin ELITE input file into a PDB format
 * for input to the BOUT++ grid generator
 *************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "pdb.h"

#define DEFAULT_INPUT "bout.eqin"
#define MAXLINELEN 256

int read_int(FILE *f)
{
  int var;
  
  if(fscanf(f, "%d", &var) != 1) {
    printf("Error: Premature end of file\n");
    exit(1);
  }
  return(var);
}

double read_double(FILE *f)
{
  double var;
  
  if(fscanf(f, "%le", &var) != 1) {
    printf("Error: Premature end of file\n");
    exit(1);
  }
  return(var);
}

double **dmatrix(int nx, int ny)
{
  double **m;
  int i;

  m = (double**) malloc(sizeof(double*)*nx);
  m[0] = (double*) malloc(sizeof(double)*nx*ny);

  for(i=1;i<nx;i++)
    m[i] = m[i-1] + ny;

  return(m);
}

#define WRITE1D(name, var) \
  if(!PD_write_alt(fout, name, "double", var, 1, ind)) {         \
    fprintf(stderr, "Error: Could not write to file '%s' (%s)\n", \
            outfile, PD_get_error());                             \
    return(1);  \
  }

/* NOTE: Need to reverse the indices, and flip theta direction */
#define WRITE2D(name, var)                                           \
  if(!PD_write_alt(fout, name, "double", *var, 2, ind)) {            \
    fprintf(stderr, "Error: Could not write to file %s (%s)\n",      \
	    outfile, PD_get_error());                                \
    return(1);                                                       \
  }

int main(int argc, char **argv)
{
  char *infile, *outfile; /* Input and output filenames */
  FILE *fin;
  char buffer[MAXLINELEN];
  int npsi, npol; /* Size in psi and poloidal directions */
  int i, j;

  /* Values read from the file */
  double *psi     = (double*) NULL;
  double *pprime  = (double*) NULL;
  double *f       = (double*) NULL;
  double *ffprime = (double*) NULL;
  double *q       = (double*) NULL;
  double **R      = (double**) NULL;
  double **Z      = (double**) NULL;
  double **Bp     = (double**) NULL;
  double *ne      = (double*) NULL;
  double *Te      = (double*) NULL;
  double *Ti      = (double*) NULL;

  /* Derived quantities */
  double *p;   /* Pressure */
  double **Bt; /* Toroidal field */

  PDBfile *fout; /* PDB output file */
  long ind[6];   /* Indices */


  /* Default input filename */
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

  /* Open input file */
  if((fin = fopen(infile, "rt")) == (FILE*) NULL) {
    printf("Error opening input file '%s'\n", infile);
    return(1);
  }

  /* Read the EQIN data file */

  fgets(buffer, MAXLINELEN-1, fin); /* First description line */
  printf("First line: %s\n", buffer);

  npsi = read_int(fin);
  npol = read_int(fin);
  
  printf("Size of grid: %d by %d\n", npsi, npol);

  do {
    /* Get description line */
    fgets(buffer, MAXLINELEN-1, fin);
    
    if(strncasecmp(buffer, "psi", 3) == 0) {
      /* Read in PSI */
      
      if(psi != (double*) NULL) {
	printf("Error: Psi defined twice!\n");
	return(1);
      }
      psi = (double*) malloc(sizeof(double)*npsi);
      
      for(i=0;i<npsi;i++)
	psi[i] = read_double(fin);
      
    }else if(strncasecmp(buffer, "pprime", 6) == 0) {
      /* PPRIME (dPdpsi) */
      
      if(pprime != (double*) NULL) {
	printf("Error: pprime defined twice!\n");
	return(1);
      }

      pprime = (double*) malloc(sizeof(double)*npsi);
      
      for(i=0;i<npsi;i++)
	pprime[i] = read_double(fin);
      
    }else if(strncasecmp(buffer, "f(psi)", 6) == 0) {
      /* F  ( = R*Bt) */
      
      if(f != (double*) NULL) {
	printf("Error: f defined twice!\n");
	return(1);
      }

      f = (double*) malloc(sizeof(double)*npsi);
      
      for(i=0;i<npsi;i++)
	f[i] = read_double(fin);
      
    }else if(strncasecmp(buffer, "ffprime", 7) == 0) {
      /* F*Fprime */
      
      if(ffprime != (double*) NULL) {
	printf("Error: ffprime defined twice!\n");
	return(1);
      }

      ffprime = (double*) malloc(sizeof(double)*npsi);
      
      for(i=0;i<npsi;i++)
	ffprime[i] = read_double(fin);
    }else if(strncasecmp(buffer, "q", 1) == 0) {
      /* Safety factor */
      
      if(q != (double*) NULL) {
	printf("Error: q defined twice!\n");
	return(1);
      }

      q = (double*) malloc(sizeof(double)*npsi);
      
      for(i=0;i<npsi;i++)
	q[i] = read_double(fin);
    }else if(strncasecmp(buffer, "R", 1) == 0) {
      /* Major radius */
      
      if(R != (double**) NULL) {
	printf("Error: R defined twice!\n");
	return(1);
      }

      R = dmatrix(npsi, npol);
      
      for(j=0;j<npol;j++)
	for(i=0;i<npsi;i++)
	  R[i][j] = read_double(fin);

    }else if(strncasecmp(buffer, "Z", 1) == 0) {
      /* Height */
      
      if(Z != (double**) NULL) {
	printf("Error: Z defined twice!\n");
	return(1);
      }

      Z = dmatrix(npsi, npol);
      
      for(j=0;j<npol;j++)
	for(i=0;i<npsi;i++)
	  Z[i][j] = read_double(fin);
    }else if(strncasecmp(buffer, "Bp", 1) == 0) {
      /* Poloidal field */
      
      if(Bp != (double**) NULL) {
	printf("Error: Bp defined twice!\n");
	return(1);
      }

      Bp = dmatrix(npsi, npol);
      
      for(j=0;j<npol;j++)
	for(i=0;i<npsi;i++)
	  Bp[i][j] = read_double(fin);
    }else if(strncasecmp(buffer, "ne ", 3) == 0) {
      /* Density */
      
      if(ne != (double*) NULL) {
	printf("Error: ne defined twice!\n");
	return(1);
      }

      ne = (double*) malloc(sizeof(double)*npsi);
      
      for(i=0;i<npsi;i++)
	ne[i] = read_double(fin);
    }else if(strncasecmp(buffer, "Te ", 3) == 0) {
      /* Electron temperature */
      
      if(Te != (double*) NULL) {
	printf("Error: Te defined twice!\n");
	return(1);
      }

      Te = (double*) malloc(sizeof(double)*npsi);
      
      for(i=0;i<npsi;i++)
	Te[i] = read_double(fin);
    }else if(strncasecmp(buffer, "Ti ", 3) == 0) {
      /* Ion temperature */
      
      if(Ti != (double*) NULL) {
	printf("Error: Ti defined twice!\n");
	return(1);
      }

      Ti = (double*) malloc(sizeof(double)*npsi);
      
      for(i=0;i<npsi;i++)
	Ti[i] = read_double(fin);
    }else {
      /* Not matched any */
      if(!isspace(buffer[0]))
	printf("Ignoring %s\n", buffer);
    }
  }while(!feof(fin));

  /* Close the input file */
  fclose(fin);

  /* Open output PDB file */

  if((fout = PD_open(outfile, "w")) == (PDBfile*) NULL) {
    fprintf(stderr, "Error: Could not open '%s' for writing (%s)\n", outfile, PD_get_error());
    return(2);
  }

  ind[0] = 0L; ind[1] = npsi-1; ind[2] = 1L;
  ind[3] = 0L; ind[4] = npol-1; ind[5] = 1L;

  /* Size of the grid */
  if(!PD_write(fout, "nx", "integer", &npsi)) {
    fprintf(stderr, "Error: Could not write to file '%s'\n", outfile);
    return(3);
  }
  if(!PD_write(fout, "ny", "integer", &npol)) {
    fprintf(stderr, "Error: Could not write to file '%s'\n", outfile);
    return(3);
  }

  /* Check which variables are present */
  
  /* Vital variables */
  if((psi == (double*) NULL) ||
     (f == (double*) NULL) ||
     (ffprime == (double*) NULL) ||
     (pprime == (double*) NULL) ||
     (q == (double*) NULL) ||
     (R == (double**) NULL) ||
     (Z == (double**) NULL)) {
    
    printf("Missing data: Require psi, f, ffprime, pprime, q, R and Z\n");
    return(1);
  }

  

  /* Calculate toroidal field */
  Bt = dmatrix(npsi, npol);
  for(i=0;i<npsi;i++)
    for(j=0;j<npol;j++)
      Bt[i][j] = f[i] / R[i][j];

  /* Calculate pressure */
  p = (double*) malloc(sizeof(double)*npsi);
  
  for(i=0;i<npsi;i++) {
    p[i] = ne[i] * (Te[i] + Ti[i]) * 1.602e-19 * (4.0*3.14159*1e-7);
  }

  WRITE1D("psi",   psi);
  WRITE1D("f",     f);
  WRITE1D("ffprime", ffprime);
  WRITE1D("mu0p", p);
  WRITE1D("mu0pprime", pprime);
  WRITE1D("qsafe", q);
  WRITE1D("Ni", ne);
  WRITE1D("Te", Te);
  WRITE1D("Ti", Ti);

  /* Write 2D variables */
  WRITE2D("Rxy", R);
  WRITE2D("Zxy", Z);

  WRITE2D("Bpxy", Bp);
  WRITE2D("Btxy", Bt);

  /* Close PDB file */
  PD_close(fout);

  return(0);
}
