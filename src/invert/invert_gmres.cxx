/*******************************************************************************
 * Inversion routine using GMRES
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 *******************************************************************************/

#include <invert_gmres.hxx>
#include <globals.hxx>

typedef int (*opfunc1D) (BoutReal *b, BoutReal *x, void* data);
int gmres_solve(BoutReal *b, BoutReal *x, int n, int m, opfunc1D A, int itmax, BoutReal tol, void* data);

int operator_2d(BoutReal *b, BoutReal *x, void *data);
int operator_2d_bndry(BoutReal *b, BoutReal *x, void *data);

typedef struct {
  opfunc2D op;
  
  int y;

  BoutReal **b2d, **x2d;
  FieldPerp *bp, *xp;
  
  int flags;
  void *data;
}iter2d;

/*******************************************************************************
 *      ITERATIVE SOLVER FUNCTIONS
 *******************************************************************************/

#define ITER_RESTART 128
#define ITER_ITMAX 1000
#define ITER_TOL 1.0e-4

int iter_solve(FieldPerp &b, FieldPerp &x, opfunc2D A, void *extra)
{
  int n, ret;
  static int size = 0;
  static iter2d data;

  n = mesh->ngx*mesh->ngz;
  if(size == 0) {
    if(size == 0) {
      data.b2d = (BoutReal**) malloc(sizeof(BoutReal*)*mesh->ngx);
      data.x2d = (BoutReal**) malloc(sizeof(BoutReal*)*mesh->ngx);
    }else {
      data.b2d = (BoutReal**) BoutRealloc(data.b2d, sizeof(BoutReal*)*mesh->ngx);
      data.x2d = (BoutReal**) BoutRealloc(data.x2d, sizeof(BoutReal*)*mesh->ngx);
    }
    size = n;
  }

  data.op = A;
  data.y = jy;
  data.data = extra;
  
  // Solve with GMRES
  ret = gmres_solve(b[0], x[0], n, ITER_RESTART, operator_2d, ITER_ITMAX, ITER_TOL, (void*) &data);

  return(ret);
}

/* Same as iter_solve_2d above, except imposes boundary conditions
   NOTE: CHANGES b input variable */
int iter_solve_bndry(FieldPerp &b, FieldPerp &x, opfunc2D A, int flags, void *extra)
{
  int n, ret;
  int jx, jz;
  BoutReal dc;
  static int size = 0;
  static iter2d data;

  n = mesh->ngx*mesh->ngz;
  if(size == 0) {
    if(size == 0) {
      data.b2d = (BoutReal**) malloc(sizeof(BoutReal*)*mesh->ngx);
      data.x2d = (BoutReal**) malloc(sizeof(BoutReal*)*mesh->ngx);
    }else {
      data.b2d = (BoutReal**) BoutRealloc(data.b2d, sizeof(BoutReal*)*mesh->ngx);
      data.x2d = (BoutReal**) BoutRealloc(data.x2d, sizeof(BoutReal*)*mesh->ngx);
    }
    size = n;
  }

  // FFT FILTER INPUT
  fft_filter_2d(b);
  
  data.op = A;
  data.y = jy;
  data.flags = flags;
  data.data = extra;
  
  /* Set boundaries to zero */
  for(jx=0;jx<mesh->ngx;jx++) {
    b[jx][ncz] = 0.0;
  }
  
  for(jx=0;jx<MXG;jx++) {
    for(jz=0;jz<=ncz;jz++) {
      /* Inner boundary */
      b[jx][jz] = 0.0;
      
      /* Outer boundary */
      b[ncx-jx][jz] = 0.0;
    }
  }

  if(flags & INVERT_ZERO_DC) {
    /* Set the DC component to zero */
    for(jx=MXG;jx<(mesh->ngx-MXG);jx++) {
      dc = 0.0;
      for(jz=0;jz<ncz;jz++) {
	dc += b[jx][jz];
      }
      dc /= (BoutReal) ncz;
      for(jz=0;jz<ncz;jz++) {
	b[jx][jz] -= dc;
      }
    }
  }

  // Solve with GMRES
  ret = gmres_solve(b[0], x[0], n, ITER_RESTART, operator_2d_bndry, ITER_ITMAX, ITER_TOL, (void*) &data);

  return(ret);
}

int iter_solve(Field3D &b, Field3D &x, opfunc2D A, void *extra)
{
  // Solve each perpendicular slice separately
  
}

int iter_solve_bndry(Field3D &b, Field3D &x, opfunc2D A, int flags, void *extra)
{
  
}

/*******************************************************************************
 * Interface operator functions
 *******************************************************************************/

int operator_2d(BoutReal *b, BoutReal *x, void *extra)
{ 
  iter2d *data;
  int i;
  int ret;

  data = (iter2d*) extra;

  /* Set 2D array pointers */

  data->b2d[0] = b;
  data->x2d[0] = x;

  for(i=1;i<mesh->ngx;i++) {
    data->b2d[i] = data->b2d[i-1] + mesh->ngz;
    data->x2d[i] = data->x2d[i-1] + mesh->ngz;
  }

  data->bp->setData(data->b2d);
  data->xp->setData(data->x2d);

  /* Call 2D operator */

  ret = data->op(*(data->bp), *(data->xp), data->data);

  return(ret);
}

int operator_2d_bndry(BoutReal *b, BoutReal *x, void *extra)
{ 
  iter2d *data;
  int ret;
  int i, j;
  BoutReal dc1;

  data = (iter2d*) extra;

  /* Set 2D array pointers */

  data->b2d[0] = b;
  data->x2d[0] = x;

  for(i=1;i<mesh->ngx;i++) {
    data->b2d[i] = data->b2d[i-1] + mesh->ngz;
    data->x2d[i] = data->x2d[i-1] + mesh->ngz;
  }
  data->bp->setData(data->b2d);
  data->xp->setData(data->x2d);
  
  /* Call 2D operator */

  ret = data->op(*(data->bp), *(data->xp), data->data);

  /* FFT filter */
  fft_filter_2d(data->b2d);

  if(data->flags & INVERT_ZERO_DC) {
    /* Set the DC component to zero */
    for(i=MXG;i<(mesh->ngx-MXG);i++) {
      dc1 = 0.0;
      for(j=0;j<ncz;j++) {
	dc1 += data->b2d[i][j];
      }
      dc1 /= (BoutReal) ncz;
      for(j=0;j<ncz;j++) {
	data->b2d[i][j] -= dc1;
      }
    }
  }

  /* Set boundaries */

  /***** Toroidal *****/
  for(i=0;i<mesh->ngx;i++) {
    data->b2d[i][ncz] = data->x2d[i][ncz] - data->x2d[i][0]; /* Enforce periodicity */
  }

  /***** Inner Radial *****/

  for(i=0;i<MXG;i++) {
    if(data->flags & (INVERT_AC_IN_GRAD | INVERT_DC_IN_GRAD)) {
      for(j=0;j<ncz;j++) {
	data->b2d[i][j] = data->x2d[i][j] - data->x2d[i+1][j];
      }
    }else {
      for(j=0;j<ncz;j++) {
	data->b2d[i][j] = data->x2d[i][j];
      }
    }
  }
 
  /***** Outer Radial *****/

  for(i=ncx;i>(ncx-MXG);i--) {
    if(data->flags & (INVERT_AC_OUT_GRAD | INVERT_DC_OUT_GRAD)) {
      for(j=0;j<ncz;j++) {
	data->b2d[i][j] = data->x2d[i][j] - data->x2d[i-1][j];
      }
    }else {
      for(j=0;j<ncz;j++) {
	data->b2d[i][j] = data->x2d[i][j];
      }
    }
  }

  return(ret);
}

/*******************************************************************************
 * GMRES Iterative solver
 *******************************************************************************/

BoutReal norm_vector(BoutReal *b, int n)
{
  int i;
  BoutReal val = 0.0;

  for(i=0;i<n;i++)
    val += b[i]*b[i];
  
  return(sqrt(val));
}

BoutReal dot_product(BoutReal *a, BoutReal *b, int n)
{
  int i;
  BoutReal val = 0.0;

  for(i=0;i<n;i++)
    val += a[i]*b[i];

  return(val);
}

void Update(BoutReal *x, int it, BoutReal **h, BoutReal *s, BoutReal *y, BoutReal **v, int n)
{
  int i, j, p;
  
  /* y = s */
  for(i=0;i!=(it+1);i++) {
    y[i] = s[i];
  }
  
  /* backsolve */
  for(i = it; i >= 0; i--) {
    y[i] /= h[i][i];
    for(j=i-1; j >= 0; j--) {
      y[j] -= h[j][i] * y[i];
    }
  }

  for(p = 0; p != it+1; p++) {
    /* x += v_(p) * y[p] */
    for(i=0;i<n;i++)
      x[i] += v[p][i] * y[p];
  }
}

void GeneratePlaneRotation(BoutReal dx, BoutReal dy, BoutReal *cs, BoutReal *sn)
{
  BoutReal temp;
  if(dy == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if(fabs(dy) > fabs(dx)) {
    temp = dx / dy;
    *sn = 1.0 / sqrt(1.0 + temp*temp);
    *cs = temp * (*sn);
  } else {
    temp = dy / dx;
    *cs = 1.0 / sqrt(1.0 + temp*temp);
    *sn = temp * (*cs);
  }
}

void ApplyPlaneRotation(BoutReal *dx, BoutReal *dy, BoutReal cs, BoutReal sn)
{
  BoutReal temp;

  temp = *dx;
  *dx = cs * (*dx) + sn * (*dy);
  *dy = cs * (*dy) - sn * temp;
}

int gmres_solve(BoutReal *b, BoutReal *x, int n, int m, opfunc1D A, int itmax, BoutReal tol, void* data)
{
  int i;
  int it, itt, p;
  BoutReal normb, beta, resid;
  
  /* Problem array storage */
  static int size = 0, msize = 0;
  static BoutReal *y, *s, *cs, *sn;
  static BoutReal **H;
  static BoutReal *r, *w;
  static BoutReal **v;
  
  if((n < 1) || (m < 1))
    return(1);

  /************************************/
  
  /*.allocate memory if problem size increased */
  if((size < n) || (msize < m)) {
    if(size != 0) {
      free(y);
      free(s);
      free(cs);
      free(sn);
      
      free_rmatrix(H);
      
      free(r);
      free(w);
      
      free_rmatrix(v);
    }
    
    size = n;
    msize = m;

    y  = rarray(m+1);
    s  = rarray(m+1);
    cs = rarray(m+1);
    sn = rarray(m+1);
    
    H  = rmatrix(m+1, m+1);

    r = rarray(n);
    w = rarray(n);

    v = rmatrix(m+1, n);
  }

  /************************************/

  /* normb = |b| */
  normb = norm_vector(b, n);
  if(normb == 0.0)
    normb = 1.0;

  /* r = b - Ax */
  A(r, x, data);
  for(i=0;i<n;i++)
    r[i] = b[i] - r[i];

  /* beta = |r| */
  beta = norm_vector(r, n);
  
  if((resid = beta / normb) <= tol) {
    iterations = 0;
    residual = resid;
    return(0);
  }
  
  it = 1;

  while(it <= itmax) {
    /* v_(0) = r / beta */
    for(i=0;i<n;i++)
      v[0][i] = r[i] / beta;

    s[0] = beta;
    
    for(itt=0; (itt < m) && (it <= itmax); itt++, it++) {
      /* w = A*v_(itt) */
      A(w, v[itt], data);
      
      for(p=0;p<=itt;p++) {
	H[p][itt] = dot_product(w, v[p], n);
	/* w = w - H[p][itt] * v[p] */
	for(i=0;i<n;i++)
	  w[i] -= H[p][itt] * v[p][i];
      }
      
      /* H[itt+1][itt] = |w| */
      H[itt+1][itt] = norm_vector(w, n);
      
      /* v[itt+1] = w / |w| */
      for(i=0;i<n;i++)
	v[itt+1][i] = w[i] / H[itt+1][itt];
      
      for(p=0; p < itt; p++) {
	ApplyPlaneRotation(&(H[p][itt]), &(H[p+1][itt]), cs[p], sn[p]);
      }
      GeneratePlaneRotation(H[itt][itt], H[itt+1][itt], &(cs[itt]), &(sn[itt]));
      ApplyPlaneRotation(&(H[itt][itt]), &(H[itt+1][itt]), cs[itt], sn[itt]);
      s[itt+1] = 0.0;
      ApplyPlaneRotation(&(s[itt]), &(s[itt+1]), cs[itt], sn[itt]);

      if((resid = fabs(s[itt+1] / normb)) < tol) {
	Update(x, itt, H, s, y, v, n);
	iterations = it;
	residual = resid;
	return(0);
      }
    }
   
    Update(x, itt-1, H, s, y, v, n);
    /* r = b - Ax */
    A(r, x, data);
    for(i=0;i<n;i++)
      r[i] = b[i] - r[i];
    
    beta = norm_vector(r, n);
    if((resid = beta / normb) < tol) {
      iterations = it;
      residual = resid;
      return(0);
    }
  }
  iterations = it;
  residual = resid;
  return(-1);
}
