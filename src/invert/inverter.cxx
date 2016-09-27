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

#include <globals.hxx>
#include <inverter.hxx>
#include <utils.hxx>
#include <fft.hxx>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/**************************************************************************
 * Constructor / Destructor
 **************************************************************************/

Inverter::Inverter()
{
  nxgt1 = !(mesh->firstX() & (mesh->lastX()));
  parallel = nxgt1;
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
  
  int status = gmres_solve(*(b.getData()), *(x.getData()), mesh->LocalNx*mesh->LocalNz, 
			   restart, itmax, tol, iterations, residual);
  
  // Iterations and residual now set by GMRES method
  
  return status;
}

int Inverter::solve(const Field3D &b, Field3D &x, 
		    int flags, 
		    int restart, int itmax,
		    BoutReal tol)
{
  int ys = 0, ye = mesh->LocalNy-1;
  // NOTE: REFINE THIS TO ONLY SOLVE IN BOUNDARY Y CELLS
  
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
  
  //Fb.setData(&b);
  //Fx.setData(&x);
  
  if(parallel) {
    // Communicate Fx
    mesh->communicate(Fx);
  }
  
  // Need to set boundary conditions on x
  applyBoundary(Fx, bndry_flags);

  Fb = function(Fx);
}

/**************************************************************************
 * Protected functions
 **************************************************************************/

// #include <invert_laplace.hxx>
const int INVERT_DC_IN_GRAD = 1;
const int INVERT_AC_IN_GRAD = 2;
const int INVERT_DC_OUT_GRAD = 4;
const int INVERT_AC_OUT_GRAD = 8;
const int INVERT_ZERO_DC = 16;
const int INVERT_START_NEW = 32;
const int INVERT_BNDRY_ONE = 64; // Sets the width of the boundary to 1
const int INVERT_4TH_ORDER = 128; // Use band solver for 4th order in x

const int INVERT_AC_IN_LAP = 256;
const int INVERT_AC_OUT_LAP = 512;

const int INVERT_IN_SYM = 1024; // Use symmetry to enforce either zero-value or zero-gradient
const int INVERT_OUT_SYM = 2048; // Same for outer boundary
const int INVERT_IN_SET = 4096; // Set inner boundary
const int INVERT_OUT_SET = 8192; // Set outer boundary
const int INVERT_IN_RHS = 16384; // Use input value in RHS at inner boundary
const int INVERT_OUT_RHS = 32768; // Use input value in RHS at outer boundary
const int INVERT_KX_ZERO = 65536; // Zero the kx=0, n = 0 component

const int INVERT_DC_IN_LAP = 131072;

const int INVERT_BNDRY_IN_ONE = 262144;
const int INVERT_BNDRY_OUT_ONE = 524288;
const int INVERT_DC_IN_GRADPAR = 1048576;
const int INVERT_DC_IN_GRADPARINV = 2097152;

/// NOTE: This should be changed/merged with Field2D/3D boundary system
void Inverter::applyBoundary(FieldPerp &f, int flags)
{ 
  // Set boundaries in Fourier space (to be compatible with
  // invert_laplace)
  
  Coordinates *coord = mesh->coordinates();
  
  int nin = mesh->xstart; // Number of inner points
  int nout = mesh->LocalNx-mesh->xend-1; // Number of outer points
  
  int ncz = mesh->LocalNz;
  
  int jy = f.getIndex();

  // Allocate working memory
  static dcomplex **cdata = NULL;
  static BoutReal *h;
  if(cdata == NULL) {
    int size = BOUTMAX(nin, nout)+2;
    cdata = cmatrix(size, ncz/2 + 1);
    h = new BoutReal[size];
  }
  
  //////////////////////////////////////
  // Inner boundary
  
  ZFFT(f[nin+1], mesh->zShift(nin+1,jy), cdata[0]);
  ZFFT(f[nin], mesh->zShift(nin,jy), cdata[1]);
  for(int i=0;i<=nin+1;i++)
    h[i] = coord->dx(nin+1-i,jy);
  
  int mask = INVERT_DC_IN_GRAD | INVERT_AC_IN_GRAD | INVERT_AC_IN_LAP;
  calcBoundary(cdata, nin, h, flags & mask);
  
  for(int i=0;i<nin;i++)
    ZFFT_rev(cdata[2+i], mesh->zShift(nin-1-i,jy), f[nin-1-i]);
  
  //////////////////////////////////////
  // Outer boundary
  
  int xe = mesh->xend;
  ZFFT(f[xe-1], mesh->zShift(xe-1,jy), cdata[0]);
  ZFFT(f[xe], mesh->zShift(xe,jy), cdata[1]);
  for(int i=0;i<=nout+1;i++)
    h[i] = coord->dx(xe-1+i,jy);
  
  mask = INVERT_DC_OUT_GRAD | INVERT_AC_OUT_GRAD | INVERT_AC_OUT_LAP;
  calcBoundary(cdata, nout, h, flags & mask);
  
  for(int i=0;i<nout;i++)
    ZFFT_rev(cdata[2+i], mesh->zShift(xe+1+i,jy), f[xe+1+i]);
}

void Inverter::calcBoundary(dcomplex **cdata, int n, BoutReal *h, int flags)
{
  int ncz = mesh->LocalNz;
  
  // DC component
  if(flags & (INVERT_DC_IN_GRAD | INVERT_DC_OUT_GRAD)) {
    // Zero gradient
    for(int i=0;i<n;i++)
      cdata[2+i][0] = cdata[1][0];
  }else {
    // Zero value
    for(int i=0;i<n;i++)
      cdata[2+i][0] = 0.0;
  }
  
  // AC component
  if(flags & (INVERT_AC_IN_GRAD | INVERT_AC_OUT_GRAD)) {
    // Zero gradient
    for(int i=0;i<n;i++)
      for(int k=1;k<=ncz/2;k++)
	cdata[2+i][k] = cdata[1][k];
    
  }else if(flags & (INVERT_AC_IN_LAP | INVERT_AC_OUT_LAP)) {
    // Zero Laplacian
    
  }else {
    // Zero value
    for(int i=0;i<n;i++)
      for(int k=1;k<=ncz/2;k++)
	cdata[2+i][k] = 0.0;
  }
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
  
  if(parallel) {
    // Add together across processors in X
    
  }

  return(sqrt(val));
}

BoutReal Inverter::dot_product(BoutReal *a, BoutReal *b, int n)
{
  int i;
  BoutReal val = 0.0;

  for(i=0;i<n;i++)
    val += a[i]*b[i];

  if(parallel) {
    // Add together across processors in X
    
  }
  
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
