/*******************************************************
 * Data compression for BOUT dump files
 * Long-term archiving of results. Only updates when
 * variable changes by a significant amount
 *******************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "pdb.h"
#include "sdclib.h"

#include <sys/time.h>

#define FILE_MAGIC "BCD 1.0"

int get_integer(char *prompt)
{
  char buffer[256];
  int val, retry;
  retry = 1;
  do {
    printf("%s: ", prompt);
    fgets(buffer, 255, stdin);

    if(sscanf(buffer, "%d", &val) == 1)
      retry = 0;
  }while(retry);

  return(val);
}

#define DEFAULT_ABSTOL 1.0e-3
#define DEFAULT_RELTOL 1.0e-6
#define DEFAULT_ETA 1.0e-16

static double now()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

int main(int argc, char**argv)
{
  FILE *fp;
  SDCfile *out;
  PDBfile *in;
  
  char *name;
  char *infile;

  long ind4d[12];
  syment *ep;
  dimdes* dims;
  int nd, size, mint, maxt, nt, t, fnr;
  int N, order, reset, nregions;
  float *data;
  int n;
  char c;
  double run_time;
  int i;
  int tfirst;

  if(argc < 3) {
    printf("Useage: %s [options] variable file1 [file2 ... fileN]\n", argv[0]);
    return(1);
  }
  
  tfirst = 0;

  i = 1;
  while(argv[i][0] == '-') {
    // This is an option
    if(strncmp(argv[i], "-t", 2) == 0) {
      // Time is the first index
      printf("Using first index as time\n");
      tfirst = 1;
    }else {
      printf("Unrecognised option '%s'\n", argv[i]);
      exit(1);
    }
    i++;
  }

  name = argv[i]; i++;

  N = 0;

  run_time = now();

  for(fnr=i;fnr < argc;fnr++) {
    infile = argv[fnr];
    printf("Opening file %s\n", infile);

    /* Open input file */

    if((in = PD_open(infile, "r")) == NULL) {
      printf("Error: Could not open input file '%s'\n", infile);
      return(1);
    }

    /* Query size of variable */
  
    if((ep = PD_query_entry(in, name, NULL)) == (syment*) NULL) {
      printf("Error querying variable %s\n", name);
      return(1);
    }
    dims = PD_entry_dimensions(ep);
    nd = 0; /* Count number of dimensions */
    size = 1;
    while(dims != (dimdes*) NULL) {
      ind4d[3*nd] = dims->index_min;
      ind4d[3*nd+1] = dims->index_max;
      size *= dims->index_max - dims->index_min + 1;
      ind4d[3*nd+2] = 1L;
      nd++;
      dims = dims->next;
    }
    
    if(tfirst) {
      // T is the first index
      mint = ind4d[0];
      maxt = ind4d[1];
    }else {
      mint = ind4d[3*(nd-1)];
      maxt = ind4d[3*(nd-1)+1];
    }
    nt = maxt - mint + 1;

    size /= nt;

    if(N == 0) {
      /* First file */

      N = size;

      printf("Number of dimensions: %d (plus time)\n", nd-1);
      printf("Number of data points: %d\n", size);
      printf("Number of time points: %d\n", nt);
      
      order = get_integer("Maximum interpolation order");
      reset = get_integer("Reset period");
      nregions = get_integer("Number of regions");

      /* Allocate memory */
      data = (float*) malloc(sizeof(float)*size);

      /* Open output file */
      fp = fopen("output.sdc", "wb");

      /* Write header */
      n = strlen(FILE_MAGIC);
      WRITE_VAR(FILE_MAGIC, 1, n, fp, 1);
      
      c = (char) sizeof(int);
      WRITE_VAR(&c, 1, 1, fp, 1);
      c = (char) sizeof(float);
      WRITE_VAR(&c, 1, 1, fp, 1);
      c = (char) sizeof(long);
      WRITE_VAR(&c, 1, 1, fp, 1);
      
      WRITE_STRING(name, fp, 1);
      WRITE_VAR(&nd, sizeof(int), 1, fp, 1);
      if(tfirst) {
	WRITE_VAR(ind4d+3, sizeof(long), 3*(nd-1), fp, 1);
      }else
	WRITE_VAR(ind4d, sizeof(long), 3*(nd-1), fp, 1);
      WRITE_VAR(&size, sizeof(int), 1, fp, 1);

      /* Create SDC header */
      out = sdc_newfile(fp, N, order, reset, nregions);
      sdc_set_tol(out, DEFAULT_ABSTOL, DEFAULT_RELTOL, DEFAULT_ETA);
    }else {
      /* Check size matches */
      if(size != N) {
	printf("ERROR: This file has a different size (%d), expected %d\n",
	       size, N);
	return(1);
      }
    }

    for(t=mint;t<=maxt;t++) {
      if(tfirst) {
	ind4d[0] = t;
	ind4d[1] = t;
      }else {
	ind4d[3*(nd-1)] = t;
	ind4d[3*(nd-1)+1] = t;
      }
      
      printf("\rReading time-slice %4d", t);
      fflush(stdout);
      
      if (PD_read_as_alt(in, name, "float", data, ind4d) == FALSE) {
	printf("Error reading input file\n");
	return(1);
      }
      
      if(sdc_write(out, data)) {
	printf("Error writing data\n");
	return(2);
      }
    }
    printf("\n");
    
    PD_close(in);
  }

  sdc_close(out);
  fclose(fp);

  printf("Run time: %e seconds\n", now() - run_time);

  return(0);
}
