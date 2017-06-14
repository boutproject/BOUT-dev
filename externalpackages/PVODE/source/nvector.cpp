/****************************************************************
 *                                                              *
 * File          : nvector.c                                    *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, and       *
 *                 Allan G. Taylor, LLNL                        *
 * Last Modified : 10 November 1998                             *
 *--------------------------------------------------------------*
 *                                                              *
 * This is the implementation file for an MPI VECTOR            *
 * package. It contains the implementation of the N_Vector      *
 * kernels listed in nvector.h.                                 *
 * Modified for MPI.                                            *
 *                                                              *
 ****************************************************************/

/*......................................................................

                            LEGAL NOTICES

This work was performed at the University of California, Lawrence
Livermore National Laboratory (UC LLNL) under contract no.
W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy
(DOE) and The Regents of the University of California (the University)
for the operation of UC LLNL.  The rights of the Federal Government are 
reserved under Contract 48 subject to the restrictions agreed upon by the 
DOE and University as allowed under DOE Acquisition Letter 97-1.

This work was prepared as an account of work sponsored by an agency of
the United States Government.  Neither the United States Government
nor the University of California nor any of their empolyees makes any
warranty, express or implied, or assumes any liability or
responsibility for the accuracy, completeness, or usefulness of any
information, apparatus, product, or process disclosed, or represents
that its use would not infringe privately owned rights.  Reference
herein to any specific commercial products, process, or service by
trade name, trademark, manufacturer, or otherwise, does not
necessarily constitute or imply its endorsement, recommendation, or
favoring by the United States Government or the University of
California.  The views and opinions of authors expressed herein do not
necessarily state or reflect those of the United States Government or
the University of California, and shall not be used for advertising or
product endorsement purposes.

......................................................................*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nvector.h"
#include "llnltyps.h"
#include "llnlmath.h" 

namespace pvode {

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Error message */

#define BAD_N1  "PVecInitMPI-- Sum of local vector lengths differs from "
#define BAD_N2  "input global length. \n\n"
#define BAD_N    BAD_N1 BAD_N2

/* Private Helper Prototypes */

static real PVecAllReduce(real d, int op, machEnvType machEnv);
/* Reduction operations add/max/min over the processor group */
static void VCopy(N_Vector x, N_Vector z); /* z=x */
static void VSum(N_Vector x, N_Vector y, N_Vector z); /* z=x+y */
static void VDiff(N_Vector x, N_Vector y, N_Vector z); /* z=x-y */
static void VNeg(N_Vector x, N_Vector z); /* z=-x */
/* z=c(x+y) */
static void VScaleSum(real c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff(real c, N_Vector x, N_Vector y, N_Vector z); 
static void VLin1(real a, N_Vector x, N_Vector y, N_Vector z); /* z=ax+y */
static void VLin2(real a, N_Vector x, N_Vector y, N_Vector z); /* z=ax-y */
static void Vaxpy(real a, N_Vector x, N_Vector y); /* y <- ax+y */
static void VScaleBy(real a, N_Vector x); /* x <- ax */


/********************* Exported Functions ************************/


void *PVecInitMPI(MPI_Comm comm,  integer local_vec_length, 
		  integer global_vec_length, int *argc, char ***argv)
{
  int initflag, initerr;
  integer n, Nsum;
  machEnvType env;

  /* Create structure env and begin loading it */
  env = (machEnvType) malloc (sizeof *env);
  if (env == NULL) return(NULL);

  env->local_vec_length = local_vec_length;
  env->global_vec_length = global_vec_length;

  MPI_Initialized(&initflag);
  if (!initflag) {
    initerr = MPI_Init(argc,argv);
    if (initerr != MPI_SUCCESS) return(NULL);
    }
  env->init_by_user = initflag;

  env->comm = (comm == MPI_COMM_NULL) ? MPI_COMM_WORLD : comm; 

  /* If this PE is inactive, return now */
  if (local_vec_length <= 0) return(env);

  /* Compute global length as sum of local lengths */
  n = local_vec_length;
  MPI_Allreduce(&n, &Nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
  env->global_vec_length = Nsum;

  /* Check input global length against computed value */
  if (Nsum != global_vec_length) {
    printf(BAD_N);
    PVecFreeMPI(env);
    return(NULL);
    } 

  /* Return the pointer env */
  return(env);
}


void PVecFreeMPI(void *machEnv)
{
  machEnvType env;

  env = (machEnvType) machEnv;
  if (env == NULL) return;

  if (!(env->init_by_user)) MPI_Finalize();

  free(machEnv);
}

 
N_Vector N_VNew(integer N, machEnvType machEnv)
{
  N_Vector v;
  int N_local, N_global;

  if (N <= 0) return(NULL);
  if (machEnv == NULL) return(NULL);

  N_local = machEnv->local_vec_length;
  N_global = machEnv->global_vec_length;

  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  v->data = (real *) malloc(N_local * sizeof(real));
  if (v->data == NULL) {
    free(v);
    return(NULL);
  }

  v->length = N_local;
  v->global_length = N_global;

  v->machEnv = (machEnvType) malloc(sizeof *machEnv);
  if (v->machEnv == NULL) {
    free(v->data);
    free(v);
    return(NULL);
  }
  memcpy(v->machEnv,machEnv,sizeof *machEnv);
  
  return(v);
}


void N_VFree(N_Vector x)
{
  free(x->data);
  free(x->machEnv);
  free(x);
}


void N_VLinearSum(real a, N_Vector x, real b, N_Vector y, N_Vector z)
{
  integer N;
  real c, *xd, *yd, *zd;
  N_Vector v, v1, v2;
  boole test;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy(a,x,y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */
  
  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

#ifndef _OPENMP
  for (integer i=0; i < N; i++) 
    *zd++ = a * (*xd++) + b * (*yd++);
#else
  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = a * xd[i] + b * yd[i];
#endif
}


void N_VConst(real c, N_Vector z)
{
  integer N;
  real *zd;

  N = z->length;
  zd = z->data;

#ifndef _OPENMP
  for (integer i=0; i < N; i++) 
    *zd++ = c;
#else
  #pragma omp parallel for
  for (integer i=0; i < N; i++) 
    zd[i] = c;
#endif
}


void N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
  integer i, N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = (*xd++) * (*yd++);
}


void N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
  integer N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

#ifndef _OPENMP
  for (integer i=0; i < N; i++)
    *zd++ = (*xd++) / (*yd++);
#else
  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = xd[i] / yd[i];
#endif
}


void N_VScale(real c, N_Vector x, N_Vector z)
{
  integer N;
  real *xd, *zd;

  if (z == x) {       /* BLAS usage: scale x <- cx */
    VScaleBy(c, x);
    return;
  }

  if (c == ONE) {
    VCopy(x, z);
  } else if (c == -ONE) {
    VNeg(x, z);
  } else {
    N = x->length;
    xd = x->data;
    zd = z->data;
#ifndef _OPENMP
    for (integer i=0; i < N; i++) *zd++ = c * (*xd++);
#else
    #pragma omp parallel for
    for (integer i=0; i < N; i++)
      zd[i] = c * xd[i];
#endif
  }
}


void N_VAbs(N_Vector x, N_Vector z)
{
  integer N;
  real *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

#ifndef _OPENMP
  for (integer i=0; i < N; i++, xd++, zd++)
    *zd = ABS(*xd);
#else
  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = ABS(xd[i]);
#endif
}


void N_VInv(N_Vector x, N_Vector z)
{
  integer N;
  real *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;
  
#ifndef _OPENMP
  for (integer i=0; i < N; i++)
    *zd++ = ONE / (*xd++);
#else
  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = ONE / xd[i];
#endif
}


void N_VAddConst(N_Vector x, real b, N_Vector z)
{
  integer N;
  real *xd, *zd;
  
  N = x->length;
  xd = x->data;
  zd = z->data;
  
#ifndef _OPENMP
  for (integer i=0; i < N; i++) *zd++ = (*xd++) + b; 
#else
  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = xd[i] + b;
#endif
}


real N_VDotProd(N_Vector x, N_Vector y)
{
  integer loclen;
  real sum = ZERO, *xd, *yd, gsum;
  machEnvType machenv;

  loclen = x->length;
  xd = x->data;
  yd = y->data;
  machenv = x->machEnv;
  
  #pragma omp parallel for reduction(+: sum)
  for (integer i=0; i < loclen; i++)
    sum += xd[i] * yd[i];
  
  gsum = PVecAllReduce(sum, 1, machenv);
  return(gsum);
}


real N_VMaxNorm(N_Vector x)
{
  integer N;
  real max = ZERO, *xd, gmax;
  machEnvType machenv;

  N = x->length;
  xd = x->data;
  machenv = x->machEnv;

#ifndef _OPENMP
  for (integer i=0; i < N; i++, xd++) {
    if (ABS(*xd) > max) max = ABS(*xd);
  }
#else
  #pragma omp parallel 
  {
    real max2 = ZERO;
    #pragma omp for nowait
    for (integer i=0; i < N; i++) {
      if (ABS(xd[i]) > max2)
        max2 = ABS(xd[i]);
    }
    #pragma omp critical
    {
      if(max2 > max)
        max = max2;
    }
  }
#endif  
  gmax = PVecAllReduce(max, 2, machenv);
  return(gmax);
}


real N_VWrmsNorm(N_Vector x, N_Vector w)
{
  integer N, N_global;
  real sum = ZERO, *xd, *wd, gsum;
  machEnvType machenv;

  N = x->length;
  N_global = x->global_length;
  xd = x->data;
  wd = w->data;
  machenv = x->machEnv;

  #pragma omp parallel for reduction(+: sum)
  for (integer i=0; i < N; i++) {
    real prodi = xd[i] * wd[i];
    sum += prodi * prodi;
  }

  gsum = PVecAllReduce(sum, 1, machenv);
  return(RSqrt(gsum / N_global));
}

real N_VMin(N_Vector x)
{
  integer N;
  real min, *xd, gmin;
  machEnvType machenv;

  N = x->length;
  xd = x->data;
  min = xd[0];
  machenv = x->machEnv;

#ifndef _OPENMP
  for (integer i=1; i < N; i++, xd++) {
    if ((*xd) < min) min = *xd;
  }
#else
  #pragma omp parallel 
  {
    real min2 = min;
    #pragma omp for nowait
    for (integer i=0; i < N; i++) {
      if (xd[i] < min) min = xd[i];
    }
    #pragma omp critical
    {
      if(min2 < min)
        min = min2;
    }
  }
#endif

  gmin = PVecAllReduce(min, 3, machenv);
  return(gmin);
}


real N_VWL2Norm(N_Vector x, N_Vector w)
{
  integer N, N_global;
  real sum = ZERO, *xd, *wd, gsum;
  machEnvType machenv;

  N = x->length;
  N_global = x->global_length;
  xd = x->data;
  wd = w->data;
  machenv = x->machEnv;
  
  #pragma omp parallel for reduction(+: sum)
  for (integer i=0; i < N; i++) {
    real prodi = xd[i] * wd[i];
    sum += prodi * prodi;
  }

  gsum = PVecAllReduce(sum, 1, machenv);
  return(RSqrt(gsum));
}


real N_VL1Norm(N_Vector x)
{
  integer N;
  real sum = ZERO, gsum, *xd;
  machEnvType machenv;


  N = x->length;
  xd = x->data;
  machenv = x->machEnv;

  #pragma omp parallel for reduction(+: sum)
  for (integer i=0; i<N; i++) sum += ABS(xd[i]);

  gsum = PVecAllReduce(sum, 1, machenv);

  return(gsum);
}


void N_VCompare(real c, N_Vector x, N_Vector z)
{
  integer i, N;
  real *xd, *zd;
  
  N = x->length;
  xd = x->data;
  zd = z->data;
  
  #pragma omp parallel for
  for (i=0; i < N; i++) {
    zd[i] = (ABS(xd[i]) >= c) ? ONE : ZERO;
  }
}


boole N_VInvTest(N_Vector x, N_Vector z)
{
  integer i, N;
  real *xd, *zd, val, gval;
  machEnvType machenv;

  N = x->length;
  xd = x->data;
  zd = z->data;
  machenv = x->machEnv;

  val = 1;
  for (i=0; i < N; i++) {
    if (*xd == ZERO) 
      val = 0;
    else
      *zd++ = ONE / (*xd++);
  }

  gval = PVecAllReduce(val, 3, machenv);
  if (gval == 0)
    return(FALSE);
  else
    return(TRUE);
}


boole N_VConstrProdPos(N_Vector c, N_Vector x)
{
  integer i, N;
  boole test;
  real  *xd, *cd;
  machEnvType machenv;

  N =  x->length;
  xd = x->data;
  cd = c->data;
  machenv = x->machEnv;

  test = TRUE;
  for (i=0; i < N; i++, xd++,cd++) {
    if (*cd !=ZERO){
      if(*xd * *cd <= ZERO) {
	test = FALSE;
        break;
      }
    }
  }
  return((boole)PVecAllReduce((int)test,3,machenv));
}

 
void N_VPrint(N_Vector x)
{
  integer i, N;
  real *xd;

  N = x->length;
  xd = x->data;

  for (i=0; i < N; i++) printf("%g\n", *xd++);

  printf("\n");
}


/***************** Private Helper Functions **********************/

static real PVecAllReduce(real d, int op, machEnvType machEnv)
{
  /* This function does a global reduction.  The operation is
       sum if op = 1,
       max if op = 2,
       min if op = 3.
     The operation is over all processors in the group defined by
     the parameters within machEnv. */

  MPI_Comm comm;
  real out;

  comm = machEnv->comm;

  switch (op) {
   case 1: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

  return(out);
}


static void VCopy(N_Vector x, N_Vector z)
{
  integer N;
  real *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = xd[i]; 
}


static void VSum(N_Vector x, N_Vector y, N_Vector z)
{
  integer N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = xd[i] + yd[i];
}


static void VDiff(N_Vector x, N_Vector y, N_Vector z)
{
  integer N;
  real *xd, *yd, *zd;
 
  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = xd[i] - yd[i];
}


static void VNeg(N_Vector x, N_Vector z)
{
  integer N;
  real *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = -xd[i];
}


static void VScaleSum(real c, N_Vector x, N_Vector y, N_Vector z)
{
  integer N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = c * (xd[i] + yd[i]);
}


static void VScaleDiff(real c, N_Vector x, N_Vector y, N_Vector z)
{
  integer N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = c * (xd[i] - yd[i]);
}


static void VLin1(real a, N_Vector x, N_Vector y, N_Vector z)
{
  integer N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = a * xd[i] + yd[i];
}


static void VLin2(real a, N_Vector x, N_Vector y, N_Vector z)
{
  integer N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    zd[i] = a * xd[i] - yd[i];
}

static void Vaxpy(real a, N_Vector x, N_Vector y)
{
  integer N;
  real *xd, *yd;

  N = x->length;
  xd = x->data;
  yd = y->data;

  if (a == ONE) {
    #pragma omp parallel for
    for (integer i=0; i < N; i++)
      yd[i] += xd[i];
    return;
  }

  if (a == -ONE) {
    #pragma omp parallel for
    for (integer i=0; i < N; i++)
      yd[i] -= xd[i];
    return;
  }    

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    yd[i] += a * xd[i];
}

static void VScaleBy(real a, N_Vector x)
{
  integer N;
  real *xd;

  N = x->length;
  xd = x->data;

  #pragma omp parallel for
  for (integer i=0; i < N; i++)
    xd[i] *= a;
}

} //namespace
