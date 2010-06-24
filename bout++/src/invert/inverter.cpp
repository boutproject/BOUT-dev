/**************************************************************************
 * Class for non-linear inversion problems
 *
 * Uses GMRES to solve problems of form F(x)=b, either on a single processor
 * or in parallel
 * 
 * Changelog: 
 *
 * 2007-10 Ben Dudson <bd512@york.ac.uk>
 *    * Initial version. Not working yet.
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
 **************************************************************************/

#include "globals.h"
#include "inverter.h"
#include "boundary.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/**************************************************************************
 * Constructor / Destructor
 **************************************************************************/

Inverter::Inverter()
{
  
}

Inverter::~Inverter()
{
  
}

/**************************************************************************
 * Solve function
 **************************************************************************/

int Inverter::solve(const FieldPerp &b, FieldPerp &x, int flags, int restart, int itmax, BoutReal tol)
{
  int iterations;
  BoutReal residual;
  
  int status = gmres_solve(*(b.getData()), *(x.getData()), mesh->ngx*mesh->ngz, 
			   restart, itmax, tol, iterations, residual);
  
  // Iterations and residual now set by GMRES method
  
  return status;
}

int Inverter::solve(const Field3D &b, Field3D &x, 
		    int flags, 
		    int restart, int itmax,
		    BoutReal tol)
{
  int ys = mesh->ystart, ye = mesh->yend;
 
  if(MYPE_IN_CORE == 0) {
    // NOTE: REFINE THIS TO ONLY SOLVE IN BOUNDARY Y CELLS
    ys = 0;
    ye = mesh->ngy-1;
  }
  
  FieldPerp xperp;
  for(int jy=ys; jy <= ye; jy++) {
    xperp = x.slice(jy); // For starting values
    
    int ret;
    if((ret = solve(b.slice(jy), xperp, flags, restart, itmax, tol)))
      return ret;
    x = xperp;
  }
  return 0;
}

void Inverter::A(BoutReal *b, BoutReal *x)
{
  FieldPerp Fb, Fx;
  
  Fb.setData(&b);
  Fx.setData(&x);

  Fb = function(Fx);
}

/**************************************************************************
 * GMRES iterative solver
 **************************************************************************/

BoutReal Inverter::norm_vector(BoutReal *b, int n)
{
  int i;
  BoutReal val = 0.0;

  for(i=0;i<n;i++)
    val += b[i]*b[i];
  
  return(sqrt(val));
}

BoutReal Inverter::dot_product(BoutReal *a, BoutReal *b, int n)
{
  int i;
  BoutReal val = 0.0;

  for(i=0;i<n;i++)
    val += a[i]*b[i];

  return(val);
}

void Inverter::Update(BoutReal *x, int it, BoutReal **h, BoutReal *s, BoutReal *y, BoutReal **v, int n)
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

void Inverter::GeneratePlaneRotation(BoutReal dx, BoutReal dy, BoutReal *cs, BoutReal *sn)
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

void Inverter::ApplyPlaneRotation(BoutReal *dx, BoutReal *dy, BoutReal cs, BoutReal sn)
{
  BoutReal temp;

  temp = *dx;
  *dx = cs * (*dx) + sn * (*dy);
  *dy = cs * (*dy) - sn * temp;
}

int Inverter::gmres_solve(BoutReal *b, BoutReal *x, int n, int m, int itmax, BoutReal tol, int &iterations, BoutReal &residual)
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

    y  = rvector(m+1);
    s  = rvector(m+1);
    cs = rvector(m+1);
    sn = rvector(m+1);
    
    H  = rmatrix(m+1, m+1);

    r = rvector(n);
    w = rvector(n);

    v = rmatrix(m+1, n);
  }

  /************************************/

  /* normb = |b| */
  normb = norm_vector(b, n);
  if(normb == 0.0)
    normb = 1.0;

  /* r = b - Ax */
  A(r, x);
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
      A(w, v[itt]);
      
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
    A(r, x);
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
