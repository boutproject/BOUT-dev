/*************************************************************************************
 * Simulation Data Compression library
 * Compresses output from simulations, particularly those
 * where the state isn't changing all that fast
 * Designed with BOUT output in mind
 *
 * MIT LICENSE:
 *
 * Copyright (c) 2007 B.Dudson, University of York
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the Software
 * is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *************************************************************************************/

#ifndef __SDCLIB_H__
#define __SDCLIB_H__

#define SDC_MAGIC "SDC 0.1"

typedef struct {
  int n;
  
  int npoints; /* Number of points which are written to the file */
  int nextra;  /* Number of points since last point */
  int *time; /* Times for each of the points */

  float **data; /* Pointers to the data block */
  float *dptr; /* Pointer to start of data block */

  long last; /* Location to write the next address to */
}SDCregion;

typedef struct {
  FILE *fp;
  long header; /* Location of header */

  int N;      /* Size of the data array */
  int order;  /* Maximum order for the predictor */
  int reset;  /* Time-points between iframes */

  int nt; /* Number of time-points */

  float abstol;
  float reltol;
  float eta;

  /* I-frame table */
  long ifpos;   /* Location of the i-frame table */
  int niframes; /* Number of i-frames */
  int *time;    /* Time for each iframe */
  long *iframe; /* Pointers to data */
  
  /* Region descriptions */
  int nregions;
  SDCregion *region;

  int writing; /* Set to 1 if writing, 0 if just reading */
}SDCfile;

#define SDC_NULL ((SDCfile*) NULL)

/* For creating an SDC file.
   Specify size of array N, interpolation order and points between iframes */
SDCfile* sdc_newfile(FILE *fp, int N, int order, int reset, int nregions);
void sdc_set_tol(SDCfile *f, float abstol, float reltol, float eta);
int sdc_write(SDCfile *f, float *data);

/* For reading an SDC file */
SDCfile* sdc_open(FILE *fp);
int sdc_read(SDCfile *f, int t, float *data);

int sdc_close(SDCfile *f);

/* Useful macros */

#define WRITE_VAR(a, s, n, fd, ret)		  \
  if(fwrite(a, s, n, fd) != n) {                  \
    printf("\tError: Could not write to file\n"); \
    return(ret);                                  \
  }

#define WRITE_STRING(a, fd, ret)                \
  if(a == (char*) NULL) {                       \
    n = 0;                                      \
    WRITE_VAR(&n, sizeof(int), 1, fd, ret);     \
  }else {                                       \
    n = strlen(a);				\
    WRITE_VAR(&n, sizeof(int), 1, fd, ret);     \
    WRITE_VAR(a, 1, n, fd, ret);                \
  }

#define READ_VAR(a, s, n, fd, ret)              \
  if(fread(a, s, n, fd) != n) {                 \
    printf("\tError: Unexpected end of file\n");\
    return(ret);                                \
  }

#define READ_STRING(a, fd, ret)                 \
  READ_VAR(&n, sizeof(int), 1, fd, ret);        \
  if(n <= 0) {                                  \
    a = (char*) NULL;                           \
  }else {                                       \
    a = (char*) malloc(n+1);                    \
    READ_VAR(a, n, 1, fd, ret);                 \
    a[n] = 0;                                   \
  }

#endif /* __SDCLIB_H__ */
