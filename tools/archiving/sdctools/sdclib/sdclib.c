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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sdclib.h"

#define DEFAULT_IFRAME 10
#define DEFAULT_ORDER 4

#define DEFAULT_ABSTOL 1.0e-3
#define DEFAULT_RELTOL 1.0e-2
#define DEFAULT_ETA 1.0e-10

/* internal function prototypes */
int sdc_init(SDCfile *f);
void sdc_delete(SDCfile *f);
void sdc_fit(int *t, float *d, int n);
float sdc_interp(int *t, float *d, int n, int time);
void polint(int *xa, float *ya, int n, int x, float *y, float *dy);
void lagrange(int *xa, float *ya, int n, int x, float *y);
void sdc_add_iframe(SDCfile *f, int t, long pos);

/* External functions */

SDCfile* sdc_newfile(FILE *fp, int N, int order, int reset, int nregions)
{
  int n;
  char c;
  long head;
  SDCfile *file;
  
  /* Check parameters */
  
  if(nregions < 1)
    nregions = 256; /* Only need one byte to specify */

  if(reset < 1)
    reset = DEFAULT_IFRAME;

  if(order < 1)
    order = DEFAULT_ORDER;
  
  if(order > reset)
    order = reset;
  
  /* Write header into file */

  n = strlen(SDC_MAGIC);

  WRITE_VAR(SDC_MAGIC, 1, n, fp, SDC_NULL);
 
  /* Write size of variables */
  c = (char) sizeof(int);
  WRITE_VAR(&c, 1, 1, fp, SDC_NULL);
  c = (char) sizeof(float);
  WRITE_VAR(&c, 1, 1, fp, SDC_NULL);
  c = (char) sizeof(long);
  WRITE_VAR(&c, 1, 1, fp, SDC_NULL);
 
  head = ftell(fp);

  /* Description of data */
  WRITE_VAR(&N, sizeof(int), 1, fp, SDC_NULL);
  WRITE_VAR(&order, sizeof(int), 1, fp, SDC_NULL);
  WRITE_VAR(&reset, sizeof(int), 1, fp, SDC_NULL);
  n = 0;
  WRITE_VAR(&n, sizeof(int), 1, fp, SDC_NULL); /* number of time-points */
  WRITE_VAR(&n, sizeof(int), 1, fp, SDC_NULL); /* Number of i-frames */
  WRITE_VAR(&head, sizeof(long), 1, fp, SDC_NULL); /* Location of i-frame table */
  WRITE_VAR(&nregions, sizeof(int), 1, fp, SDC_NULL); /* Number of regions */

  /* Create structure */

  file = (SDCfile*) malloc(sizeof(SDCfile));
  
  /* Store settings */

  file->fp     = fp;
  file->header = head;
  file->N      = N;
  file->order  = order;
  file->reset  = reset;
  
  file->nt = 0;
  
  file->abstol = DEFAULT_ABSTOL;
  file->reltol = DEFAULT_RELTOL;
  file->eta    = DEFAULT_ETA;

  file->niframes = 0;

  file->nregions = nregions;

  file->writing = 1;

  /* Allocate memory */
  if(sdc_init(file)) {
    free(file);
    return(SDC_NULL);
  }

  return(file);
}

void sdc_set_tol(SDCfile *f, float abstol, float reltol, float eta)
{
  f->abstol = abstol;
  f->reltol = reltol;
  f->eta = eta;
}

SDCfile* sdc_open(FILE *fp)
{
  char buffer[256];
  int n;
  char c;
  SDCfile *file;
  long mark;
  
  /* Read the header */

  n = strlen(SDC_MAGIC);
  
  READ_VAR(buffer, 1, n, fp, SDC_NULL);
  
  if(strncmp(buffer, SDC_MAGIC, n) != 0) {
    printf("Error: This is not a valid SDC file\n");
    return(SDC_NULL);
  }
  
  /* Read sizes of variables */
  READ_VAR(&c, 1, 1, fp, SDC_NULL);
  if(c != sizeof(int)) {
    printf("Error: Size of integer does not match\n");
    return(SDC_NULL);
  }
  READ_VAR(&c, 1, 1, fp, SDC_NULL);
  if(c != sizeof(float)) {
    printf("Error: Size of float does not match\n");
    return(SDC_NULL);
  }
  READ_VAR(&c, 1, 1, fp, SDC_NULL);
  if(c != sizeof(long)) {
    printf("Error: Size of long does not match\n");
    return(SDC_NULL);
  }
  
  file = (SDCfile*) malloc(sizeof(SDCfile));

  file->fp = fp;
  file->header = ftell(fp);
  /* Read data description */
  READ_VAR(&(file->N)       , sizeof(int), 1, fp, SDC_NULL);
  READ_VAR(&(file->order)   , sizeof(int), 1, fp, SDC_NULL);
  READ_VAR(&(file->reset)   , sizeof(int), 1, fp, SDC_NULL);
  READ_VAR(&(file->nt)      , sizeof(int), 1, fp, SDC_NULL);
  READ_VAR(&(file->niframes), sizeof(int), 1, fp, SDC_NULL);
  READ_VAR(&(file->ifpos)   , sizeof(long), 1, fp, SDC_NULL);
  READ_VAR(&(file->nregions), sizeof(int), 1, fp, SDC_NULL);
  
  file->abstol = DEFAULT_ABSTOL;
  file->reltol = DEFAULT_RELTOL;
  file->eta    = DEFAULT_ETA;

  file->writing = 0;
  
  if(sdc_init(file)) {
    free(file);
    return(SDC_NULL);
  }

  mark = ftell(fp);

  fseek(fp, file->ifpos, SEEK_SET);
  file->time = (int*) malloc(sizeof(int)*file->niframes);
  file->iframe = (long*) malloc(sizeof(long)*file->niframes);
  READ_VAR(file->time, sizeof(int), file->niframes, fp, SDC_NULL);
  READ_VAR(file->iframe, sizeof(long), file->niframes, fp, SDC_NULL);
  fseek(fp, mark, SEEK_SET);

  printf("SDC: N = %d, nt = %d, niframes = %d, nregions = %d\n",
	 file->N, file->nt, file->niframes, file->nregions);

  for(n=0;n<file->niframes;n++)
    printf("%d: %d -> %ld\n", n, file->time[n], file->iframe[n]);

  return(file);
}

int sdc_close(SDCfile *f)
{
  int i, n;
  SDCregion *r;
  long pos, mark;
  
  if(f == SDC_NULL)
    return(1);

  if(f->writing) {
    if(( (f->nt-1) % f->reset) != 0) { /* If 0, last point written was an i-frame */
      /* Write an i-frame using the last set of data */
#ifdef DEBUG
      printf("Writing final i-frame\n");
#endif

      pos = ftell(f->fp);

      i = 0;
      WRITE_VAR(&i, sizeof(int), 1, f->fp, 2); /* Mark as an i-frame */

      for(i=0;i<f->nregions;i++) {
	r = &(f->region[i]);

	n = r->npoints + r->nextra - 1;

	WRITE_VAR(r->data[n], sizeof(float), r->n, f->fp, 2);
      }
      
      sdc_add_iframe(f, f->nt-1, pos);

      /* Update pointers */

      for(i=0;i<f->nregions;i++) {
	r = &(f->region[i]);
	
	mark = ftell(f->fp);

	if(r->last != 0) {
	  /* Need to store a pointer to the data */
	  fseek(f->fp, r->last, SEEK_SET);
	  WRITE_VAR(&pos, sizeof(long), 1, f->fp, 2);
	  
	  fseek(f->fp, mark, SEEK_SET);
	}
	r->last = mark;
	WRITE_VAR(&mark, sizeof(long), 1, f->fp, 2); /* Just a placeholder */
      }
    }
  

    /* Write the i-frame table */
    
    pos = ftell(f->fp);
    WRITE_VAR(f->time, sizeof(int), f->niframes, f->fp, 2);
    WRITE_VAR(f->iframe, sizeof(long), f->niframes, f->fp, 2);
    
    mark = ftell(f->fp);

    /* Update the header */
    
    fseek(f->fp, f->header + 3*sizeof(int), SEEK_SET);
    /* Write number of time-points, iframes and location of i-frame table */
    
    WRITE_VAR(&(f->nt), sizeof(int), 1, f->fp, 2);
    WRITE_VAR(&(f->niframes), sizeof(int), 1, f->fp, 2);
    WRITE_VAR(&pos, sizeof(long), 1, f->fp, 2);
    
    fseek(f->fp, mark, SEEK_SET);
  }

  /* Free memory */
  sdc_delete(f);

  return(0);
}

int sdc_add_data(SDCfile *f, int region, float *data)
{
  int output;
  int j, k, n;
  float val;
  SDCregion *r;
  long pos, mark;
  float *ptr;

  /* Arrays used for fitting polynomials */
  static int order = 0;
  static int *t;
  static float *d;

  if(order != f->order) {
    if(order == 0) {
      t = (int*) malloc(sizeof(int)*(f->order+1));
      d = (float*) malloc(sizeof(float)*(f->order+1));
    }else {
      t = (int*) realloc(t, sizeof(int)*(f->order+1));
      d = (float*) realloc(d, sizeof(float)*(f->order+1));
    }
    order = f->order;
  }

  r = &(f->region[region]);

  /* Copy the data into the array */

  n = r->npoints + r->nextra;
  r->time[n] = f->nt;
  for(j=0;j<r->n;j++) {
    r->data[n][j] = data[j];
  }
  
  output = 0;

  if(r->nextra > 0) {
    /* Check if interpolation still works */
    for(j=0;j<r->n;j++) {
      /* Fit polynomial */
      for(k=0;k<r->npoints;k++) {
	t[k] = r->time[k];
	d[k] = r->data[k][j];
      }
      t[r->npoints] = f->nt; /* Current time */
      d[r->npoints] = data[j];
      
      /*
      for(k=0;k<=r->npoints;k++)
	printf("Poly %d: %d, %f\n",k, t[k], d[k]); 
      */

      /* Check each of the data points */
      for(k=0;k<r->nextra;k++) {
	
	val = sdc_interp(t, d, r->npoints+1, r->time[r->npoints + k]);

	//printf("Extra: %d: %d, %f, %f\n", k, r->time[r->npoints+k], r->data[r->npoints+k][j], val);
	
	val = fabsf(val - r->data[r->npoints+k][j]);

	//printf("=> %f, %f || %f %f\n", val, f->abstol, 
	//     (val/(fabs(r->data[r->npoints+k][j]) + f->eta)), f->reltol);
	
	if((val > f->abstol) || ((val/(fabs(r->data[r->npoints+k][j]) + f->eta)) > f->reltol)) {
	  /* Exceeded error bounds */
	  output = 1;
	  break;
	}
      }
      
      if(output)
	break;
    }
  }

  if(output) {
    /* Exceeded error bounds.
       Need to store previous time-point (last known good) */

#ifdef DEBUG
    printf("Writing block %d (%d, %d)\n", region, r->npoints, r->nextra);
#endif

    pos = ftell(f->fp);
    
    WRITE_VAR(&(r->time[n-1]), sizeof(int), 1, f->fp, 2); /* Write the time */
    WRITE_VAR(r->data[n-1], sizeof(float), r->n, f->fp, 2); /* the data */
    
    mark = ftell(f->fp);
    
    if(r->last != 0) {
      /* Need to store pointer to this data block */
      fseek(f->fp, r->last, SEEK_SET);
      WRITE_VAR(&pos, sizeof(long), 1, f->fp, 2); 
      fseek(f->fp, mark, SEEK_SET);
    }
    r->last = mark;
    WRITE_VAR(&mark, sizeof(long), 1, f->fp, 2); /* Just a placeholder */
    
    /* re-arrange the data pointers */
    if(r->npoints == f->order) {
      /* Shuffle to keep npoints = r->order */
      ptr = r->data[0];
      for(j=1;j<r->npoints;j++) {
	r->data[j-1] = r->data[j];
	r->time[j-1] = r->time[j];
      }
      r->data[r->npoints-1] = r->data[n-1];
      r->time[r->npoints-1] = r->time[n-1];
      r->data[n-1] = ptr;
      
      /* Move data[n] to just after the stored points */
      ptr = r->data[r->npoints];
      r->data[r->npoints] = r->data[n];
      r->data[n] = ptr;
      r->time[r->npoints] = r->time[n];
	  
    }else {
      ptr = r->data[r->npoints];
      r->data[r->npoints] = r->data[n-1];
      r->data[n-1] = ptr;
      r->time[r->npoints] = r->time[n-1];
      r->npoints++;
	  
      ptr = r->data[r->npoints];
      r->data[r->npoints] = r->data[n];
      r->data[n] = ptr;
      r->time[r->npoints] = r->time[n];
	  
    }
    r->nextra = 1;
  }else {
#ifdef DEBUG
    printf("Skipping block %d (%d, %d)\n", region, r->npoints, r->nextra);
#endif
    /* Can still interpolate */
    r->nextra++;
  }
  return(0);
}

int sdc_write(SDCfile *f, float *data)
{
  long pos, mark;
  int i;
  SDCregion *r;
  float *ptr;  

  if(f->writing == 0) {
    printf("Error: Cannot append to files yet - sorry.\n");
    return(1);
  }
  
  /* Add data to each region. Only writes data when necessary */
  ptr = data;
  for(i=0;i<f->nregions;i++) {
    if(sdc_add_data(f, i, ptr))
      return(1);
    ptr += f->region[i].n;
  }

  if((f->nt % f->reset) == 0) {
    /* This is an i-frame. Resets polynomial fitting */
#ifdef DEBUG
    printf("Writing i-frame\n");
#endif
    
    /* Get the current position */
    pos = ftell(f->fp);
    
    i = 0;
    WRITE_VAR(&i, sizeof(int), 1, f->fp, 2); /* Mark as an i-frame */
    WRITE_VAR(data, sizeof(float), f->N, f->fp, 2);
    
    /* Record the location of this i-frame */
    sdc_add_iframe(f, f->nt, pos);
    
    /* Sort out the regions */
    for(i=0;i<f->nregions;i++) {
      r = &(f->region[i]);

      mark = ftell(f->fp);
      
      if(r->last != 0) {
	/* Need to store a pointer to the data */
	fseek(f->fp, r->last, SEEK_SET);
	WRITE_VAR(&pos, sizeof(long), 1, f->fp, 2);

	fseek(f->fp, mark, SEEK_SET);
      }
      r->last = mark;
      WRITE_VAR(&mark, sizeof(long), 1, f->fp, 2); /* Just a placeholder */

      r->time[0] = f->nt;
      ptr = r->data[0];
      r->data[0] = r->data[r->npoints + r->nextra-1];
      r->data[r->npoints + r->nextra-1] = ptr;
      
      r->npoints = 1;
      r->nextra = 0;
    }
  }

  f->nt++;

  return(0);
}

int sdc_read(SDCfile *f, int t, float *data)
{
  int i, j, k;
  int fnr;
  int n, p;
  int start, end;
  SDCregion *r;
  float *ptr;

  /* arrays for interpolation */
  static int order = 0;
  static float *da;

  if(f == SDC_NULL)
    return(1);

  if(order != f->order) {
    if(order == 0) {
      da = (float*) malloc(sizeof(float)*(f->order+1));
    }else {
      da = (float*) realloc(da, sizeof(float)*(f->order+1));
    }
    order = f->order;
  }

  if((t < 0) || (t >= f->nt)) {
    printf("Error: time index %d out of bounds (0 -> %d)\n", t, f->nt-1);
    return(1);
  }

#ifdef DEBUG
  printf("Reading t = %d\n", t);
#endif

  /* Check if this point's an i-frame */
  
  if(t == f->nt-1) {
    fnr = f->niframes-1;
  }else {
    fnr = -1;
    do {
      fnr++;
    }while(f->time[fnr] <= t);
    fnr--; /* This i-frame number is the one just before (or equal to) t */ 
  }

  if(f->time[fnr] == t) {
    /* This time is an i-frame */
#ifdef DEBUG
    printf(" => i-frame %d\n", fnr);
#endif

    fseek(f->fp, f->iframe[fnr], SEEK_SET);
    
    READ_VAR(&n, sizeof(int), 1, f->fp, 1);
    if(n != 0) {
      printf("Error: File corrupted: expected i-frame\n");
      return(1);
    }
    READ_VAR(data, sizeof(float), f->N, f->fp, 1);
    
    return(0);
  }

  /* Not an iframe. Check time ranges in regions */
  
  if(f->region[0].npoints < 2) {
    start = end = -1;
  }else {
    start = f->region[0].time[f->region[0].npoints-2];
    end = f->region[0].time[f->region[0].npoints-1];
    for(i=1;i<f->nregions;i++) {
      r = &(f->region[i]);
      /* Want to know the latest start and the earliest end */
      
      if(r->time[r->npoints-2] > start)
	start = r->time[r->npoints-2];
      
      if(r->time[r->npoints-1] < end)
	end = r->time[r->npoints-1];
    }
  }
  
  if( (t < start) || (f->time[fnr] >= end) ) {
    /* Need to reset to the new i-frame */
#ifdef DEBUG
    printf(" => Reset to %d (%d)\n", fnr, f->time[fnr]);
#endif

    fseek(f->fp, f->iframe[fnr], SEEK_SET);

    /* Read the marker */
    READ_VAR(&n, sizeof(int), 1, f->fp, 1);
    if(n != 0) {
      printf("Error: Expected i-frame\n");
      return(1);
    }

    /* Read the data */
    for(i=0;i<f->nregions;i++) {
      r = &(f->region[i]);
      r->time[0] = f->time[fnr];
      READ_VAR(r->data[0], sizeof(float), r->n, f->fp, 1);
    }
    
    /* Read pointers */
    for(i=0;i<f->nregions;i++) {
      r = &(f->region[i]);
      READ_VAR(&(r->last), sizeof(long), 1, f->fp, 1);
    }

    /* Each region needs at least two data points
       read the next data block */
    p = 0;
    for(i=0;i<f->nregions;i++) {
      r = &(f->region[i]);
      
      fseek(f->fp, r->last, SEEK_SET);
      
      READ_VAR(&(r->time[1]), sizeof(int), 1, f->fp, 1);

      if(r->time[1] == 0) {
	/* This is an iframe (fnr+1) */
	r->time[1] = f->time[fnr+1];

	fseek(f->fp, (long) p*sizeof(float), SEEK_CUR);
	READ_VAR(r->data[1], sizeof(float), r->n, f->fp, 1);
	r->last = -1; /* Should never be used */
      }else {
	/* A single data block */
	READ_VAR(r->data[1], sizeof(float), r->n, f->fp, 1);
	READ_VAR(&(r->last), sizeof(long), 1, f->fp, 1);
      }

      r->npoints = 2;

#ifdef DEBUG
      printf("    region %d: (%d, %f), (%d, %f)\n", i,
	     r->time[0], r->data[0][0],
	     r->time[1], r->data[1][0]);
#endif

      p += r->n;
    }
  }

  /* Now scan forwards until time t is in the range of each region */

  p = 0; /* Location in output array */
  for(i=0;i<f->nregions;i++) {
    r = &(f->region[i]);
    
#ifdef DEBUG
    printf(" => Scanning region %d\n", i);
#endif
    
    while((r->time[r->npoints-1] < t) && (r->time[r->npoints-1] != 0)) {
      /* Need to get next time-point */
      
#ifdef DEBUG
      printf("    Reading (%d ", r->time[r->npoints-1]);
#endif

      fseek(f->fp, r->last, SEEK_SET);
      READ_VAR(&(r->time[r->npoints]), sizeof(int), 1, f->fp, 1);
      
#ifdef DEBUG
      printf("-> %d", r->time[r->npoints]);
#endif
      
      if(r->time[r->npoints] == 0) {
	/* This is an iframe (fnr+1) */
	r->time[r->npoints] = f->time[fnr+1];
#ifdef DEBUG
	printf(" [iframe %d, t=%d])\n", fnr+1, f->time[fnr+1]);
#endif

	fseek(f->fp, (long) p*sizeof(float), SEEK_CUR);
	READ_VAR(r->data[r->npoints], sizeof(float), r->n, f->fp, 1);
	r->last = -1; /* Should never be used */
      }else {
#ifdef DEBUG
	printf(")\n");
#endif
	/* A single data block */
	READ_VAR(r->data[r->npoints], sizeof(float), r->n, f->fp, 1);
	READ_VAR(&(r->last), sizeof(long), 1, f->fp, 1);
      }
      

      if(r->npoints == (order + 1)) {
	/* Need to re-arrange to keep order+1 points in the list */
	ptr = r->data[0];
	for(j=0;j<r->npoints;j++) {
	  r->data[j] = r->data[j+1];
	  r->time[j] = r->time[j+1];
	}
	r->data[r->npoints] = ptr;
      }else {
	r->npoints++;
      }
    }
    
    if(r->time[r->npoints-1] == t) {
      /* Just copy the data across */
      for(j=0;j<r->n;j++) {
	data[p] = r->data[r->npoints-1][j];
	p++;
      }
    }else{
      /* Need to interpolate */
      for(j=0;j<r->n;j++) {
	for(k=0;k<r->npoints;k++) {
	  da[k] = r->data[k][j];
	}
	data[p] = sdc_interp(r->time, da, r->npoints, t);
	p++;
      }
    }
  }

  return(0);
}

/********************** INTERNAL FUNCTIONS ***************************/

int sdc_init(SDCfile *f)
{
  int n, i, j;
  //int p;
  SDCregion *r;

  if(f->nregions < 1)
    return(1);

  if(f->reset < 1)
    return(1);

  /* Region descriptions */
  f->region = (SDCregion*) malloc(sizeof(SDCregion)*f->nregions);

  if((f->N % f->nregions) == 0) {
    /* Divides equally into regions */
    n = f->N / f->nregions;
    for(i=0;i<f->nregions;i++)
      f->region[i].n = n;
  }else {
    /* First region different size to the others */
    n = f->N % (f->nregions-1);
    f->region[0].n = n; 
    n = (f->N - n) / (f->nregions-1);
    for(i=1;i<f->nregions;i++)
      f->region[i].n = n;
  }

  //p = 0;
  for(i=0;i<f->nregions;i++) {
    r = &(f->region[i]);
    
    /* Allocate data storage */
    r->data = (float**) malloc(sizeof(float*)*(f->reset+1));
    r->dptr = r->data[0] = (float*) malloc(sizeof(float)*r->n*(f->reset+1));
    for(j=1;j<=f->reset;j++) {
      r->data[j] = r->data[j-1] + r->n;
    }
    r->time = (int*) malloc(sizeof(int)*(f->reset+1));
    
    r->npoints = 0;
    r->nextra = 0;
    r->last = 0;

    //#ifdef DEBUG
    //printf("Region %d: %d (%d -> %d)\n", i, r->n, p, p+r->n - 1);
    //#endif
    //p += r->n;
  }

  return(0);
}

/* Free memory */
void sdc_delete(SDCfile *f)
{
  int i;

  if(f == SDC_NULL)
    return;

  if(f->nregions > 0) {
    for(i=0;i<f->nregions;i++) {
      free(f->region[i].dptr);
      free(f->region[i].data);
      free(f->region[i].time);
    }
    free(f->region);
  }
}

float sdc_interp(int *t, float *d, int n, int time)
{
  float y;
  
  //polint(t, d, n, time, &y, &dy);
  lagrange(t, d, n, time, &y);

  return(y);
}

void sdc_add_iframe(SDCfile *f, int t, long pos)
{
  if(f->niframes == 0) {
    f->time = (int*) malloc(sizeof(int));
    f->iframe = (long*) malloc(sizeof(long));
  }else {
    f->time = (int*) realloc(f->time, sizeof(int)*(f->niframes+1));
    f->iframe = (long*) realloc(f->iframe, sizeof(long)*(f->niframes+1));
  }
  
  f->time[f->niframes] = t;
  f->iframe[f->niframes] = pos;
  f->niframes++;
}


/* Polynomial interpolation/extrapolation: Input xa and ya, returns value y at x with error estimate dy */

/* SOMETHING WRONG WITH THIS ROUTINE! */
void polint(int *xa, float *ya, int n, int x, float *y, float *dy)
{
  int i, m, ns=0;
  float den, dif, dift, ho, hp, w;

  static int pol_n = 0;
  static float *c, *d;

  //printf("Interpolating %d\n", n);
  
  if(n > pol_n) {
    if(pol_n == 0) {
      c = (float*) malloc(sizeof(float)*n);
      d = (float*) malloc(sizeof(float)*n);
    }else {
      c = (float*) realloc(c, sizeof(float)*n);
      d = (float*) realloc(d, sizeof(float)*n);
    }
    pol_n = n;
  }

  dif = fabsf( (float) (x - xa[0]));
  for(i=0;i<n;i++) {
    /* Find closest index */
    if((dift = fabsf((float) (x - xa[i]))) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  
  *y = ya[ns--];
  for(m=1;m<n;m++) {
    for(i=1;i<=(n-m);i++) {
      ho = (float) xa[i-1] - x;
      hp = (float) xa[i+m-1] - x;
      w = c[i] - d[i-1];
      /* Two xa's within roundoff */
      if( (den = ho-hp) == 0.0) printf("Error in polint\n");
      den = w / den;
      d[i-1] = hp*den;
      c[i-1] = ho*den;
    }
    //printf("%d, %f, %d, %f\n", m, *y, ns, c[ns]);
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[(ns--)-1]));
  }
}


/* Just a naive implementation of lagrange's formula */
void lagrange(int *xa, float *ya, int n, int x, float *y)
{
  int i,j;
  float v;

  *y = 0.0;
  for(i=0;i<n;i++) {
    v = 1.0;
    for(j=0;j<n;j++) {
      if(j != i) {
	v *= ((float) (x-xa[j])) / ((float) (xa[i] - xa[j]));
      }
    }
    *y += v*ya[i];
  }
}
