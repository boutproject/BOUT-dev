/*!
 * \file invert_laplace.cpp
 *
 * \brief Perpendicular Laplacian inversion using FFT and Tridiagonal solver
 *
 * Equation solved is \f$\nabla^2_\perp x + (1./c)\nabla_perp c\cdot\nabla_\perp x + a x = b \f$, where
 * \f$x\f$ and \f$x\f$ are perpendicular (X-Z) or 3D fields, 
 * and \f$a\f$ is a 2D field.
 * 
 * Flags control the boundary conditions (see header file)
 *
 * Parallel inversion done using two methods
 * - Either a simple parallelisation of the serial algorithm (same operations). Reasonably
 *   parallel as long as MYSUB > NXPE
 * - (EXPERIMENTAL) The Parallel Diagonally Dominant (PDD) algorithm. This doesn't seem
 *   to work properly for some simulations (works ok for some benchmarks).
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
 */

#include "mpi.h"

#include "invert_laplace.h"
#include "bout_types.h"
#include "globals.h"
#include "fft.h"
#include "utils.h"
#include "dcomplex.h"
#include "mesh_topology.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lapack_routines.h" // Tridiagonal & band inversion routines

// This was defined in nvector.h
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

/**********************************************************************************
 *                                 INITIALISATION
 **********************************************************************************/

int laplace_maxmode; ///< The maximum Z mode to solve for
bool invert_async_send; ///< If true, use asyncronous send in parallel algorithms
bool invert_use_pdd; ///< If true, use PDD algorithm
bool invert_low_mem;    ///< If true, reduce the amount of memory used
bool laplace_all_terms; // applies to Delp2 operator and laplacian inversion
bool laplace_nonuniform; // Non-uniform mesh correction

/// Laplacian inversion initialisation. Called once at the start to get settings
int invert_init()
{
  real filter; ///< Fraction of Z modes to filter out. Between 0 and 1

  output.write("Initialising Laplacian inversion routines\n");

  // Communication options
  options.setSection("comms");
  options.get("async", invert_async_send, true);
  
  // Inversion options
  options.setSection("laplace");
  options.get("filter", filter, 0.2);
  options.get("low_mem", invert_low_mem, false);
  options.get("use_pdd", invert_use_pdd, false);
  options.get("all_terms", laplace_all_terms, false); 
  OPTION(laplace_nonuniform, false);

  if(NXPE > 1) {
    if(invert_use_pdd) {
      output.write("\tUsing PDD algorithm\n");
    }else
      output.write("\tUsing parallel Thomas algorithm\n");
  }else
    output.write("\tUsing serial algorithm\n");

  // convert into an integer
  laplace_maxmode = ROUND((1.0 - filter) * ((double) (ncz / 2)));

  options.get("max_mode", laplace_maxmode, laplace_maxmode);
  
  if(laplace_maxmode < 0) laplace_maxmode = 0;
  if(laplace_maxmode > ncz/2) laplace_maxmode = ncz/2;
  
  // Broadcast this since rounding errors could cause mismatch across processors
  // THIS LINE CAUSES SEGFAULT ON LLNL GRENDEL
  //MPI_Bcast(&laplace_maxmode, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return 0;
}

/// Returns the coefficients for a tridiagonal matrix for laplace. Used by Delp2 too
void laplace_tridag_coefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b, dcomplex &c, const Field2D *ccoef)
{
  real coef1, coef2, coef3, coef4, coef5, kwave;
  
  kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
  
  coef1=mesh->g11[jx][jy]/(SQ(dx[jx][jy])); ///< X 2nd derivative
  coef2=mesh->g33[jx][jy];                  ///< Z 2nd derivative
  coef3=2.*mesh->g13[jx][jy]/(dx[jx][jy]);  ///< X-Z mixed derivative
  
  coef4 = 0.0;
  coef5 = 0.0;
  if(laplace_all_terms) {
    coef4 = mesh->G1[jx][jy] / (2.0*dx[jx][jy]); // X 1st derivative
    coef5 = mesh->G3[jx][jy];
  }

  if(laplace_nonuniform) {
    // non-uniform mesh correction
    if((jx != 0) && (jx != (ngx-1))) {
      //coef4 += mesh->g11[jx][jy]*0.25*( (1.0/dx[jx+1][jy]) - (1.0/dx[jx-1][jy]) )/dx[jx][jy]; // SHOULD BE THIS (?)
      coef4 -= 0.25*((dx[jx+1][jy] - dx[jx-1][jy])/dx[jx][jy])*coef1; // BOUT-06 term
    }
  }

  if(ccoef != NULL) {
    // A first order derivative term
    
    if((jx > 0) && (jx < (ngx-1)))
      coef4 += mesh->g11[jx][jy] * 0.25 * ((*ccoef)[jx+1][jy] - (*ccoef)[jx-1][jy]) / SQ(dx[jx][jy]*((*ccoef)[jx][jy]));
  }

  if(ShiftXderivs && IncIntShear) {
    // d2dz2 term
    coef2 += mesh->g11[jx][jy] * IntShiftTorsion[jx][jy] * IntShiftTorsion[jx][jy];
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }
  
  a = dcomplex(coef1 - coef4,-kwave*coef3);
  b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,kwave*coef5);
  c = dcomplex(coef1 + coef4,kwave*coef3);
}

/**********************************************************************************
 *                                 SERIAL CODE
 **********************************************************************************/

/// Perpendicular laplacian inversion (serial)
/*!
 * Inverts an X-Z slice (FieldPerp) using band-diagonal solvers
 * This code is only for serial i.e. NXPE == 1
 */
int invert_laplace_ser(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *ccoef=NULL)
{
  int ix, jy, iz;
  static dcomplex **bk = NULL, *bk1d;
  static dcomplex **xk, *xk1d;
  int xbndry; // Width of the x boundary
  
  real coef1=0.0, coef2=0.0, coef3=0.0, coef4=0.0, coef5=0.0, coef6=0.0, kwave, flt;

  if(NXPE != 1) {
    output.write("Error: invert_laplace only works for NXPE = 1\n");
    return 1;
  }
  
  x.Allocate();

  jy = b.getIndex();
  x.setIndex(jy);

  if(bk == NULL) {
    // Allocate memory

    bk = cmatrix(ngx, ncz/2 + 1);
    bk1d = new dcomplex[ngx];
    
    xk = cmatrix(ngx, ncz/2 + 1);
    xk1d = new dcomplex[ngx];
  }

  xbndry = MXG;
  if(flags & INVERT_BNDRY_ONE)
    xbndry = 1;

  for(ix=0;ix<=ncx;ix++) {
    // for fixed ix,jy set a complex vector rho(z)
    
    ZFFT(b[ix], zShift[ix][jy], bk[ix]);
  }
  
  if(flags & INVERT_IN_SET) {
    // Setting the inner boundary from x
    
    for(ix=0;ix<xbndry;ix++)
      ZFFT(x[ix], zShift[ix][jy], xk[ix]);
  }

  if(flags & INVERT_OUT_SET) {
    // Setting the outer boundary from x
    
    for(ix=0;ix<xbndry;ix++)
      ZFFT(x[ncx-ix], zShift[ncx-ix][jy], xk[ncx-ix]);
  }
  
  if(flags & INVERT_4TH_ORDER) { // Not implemented for parallel calculations
    // Use band solver - 4th order

    static dcomplex **A = (dcomplex**) NULL;
    int xstart, xend;

    if(A == (dcomplex**) NULL)
      A = cmatrix(ngx, 5);
    
    // Get range for 4th order: Need at least 2 each side
    if(xbndry > 1) {
      xstart = xbndry;
      xend = ngx-1-xbndry;
    }else {
      xstart = 2;
      xend = ngx-2;
    }

    for(iz=0;iz<=ncz/2;iz++) {
      // solve differential equation in x
    

      ///////// PERFORM INVERSION /////////
      
      // shift freqs according to FFT convention
      kwave=iz*2.0*PI/zlength; // wave number is 1/[rad]
      
      if (iz>laplace_maxmode) flt=0.0; else flt=1.0;

      // set bk1d
      for(ix=0;ix<=ncx;ix++)
	bk1d[ix] = bk[ix][iz]*flt;

      // Fill in interior points

      for(ix=xstart;ix<=xend;ix++) {

	// Set coefficients
	coef1 = mesh->g11[ix][jy];  // X 2nd derivative
	coef2 = mesh->g33[ix][jy];  // Z 2nd derivative
	coef3 = mesh->g13[ix][jy];  // X-Z mixed derivatives
	coef4 = 0.0;          // X 1st derivative
	coef5 = 0.0;          // Z 1st derivative
	coef6 = 0.0;          // Constant

	if(a != (Field2D*) NULL)
	  coef6 = (*a)[ix][jy];
	
	if(laplace_all_terms) {
	  coef4 = mesh->G1[ix][jy];
	  coef5 = mesh->G3[ix][jy];
	}

	if(laplace_nonuniform) {
	  // non-uniform mesh correction
	  if((ix != 0) && (ix != (ngx-1)))
	    coef4 += mesh->g11[ix][jy]*( (1.0/dx[ix+1][jy]) - (1.0/dx[ix-1][jy]) )/(2.0*dx[ix][jy]);
	}

	if(ccoef != NULL) {
	  // A first order derivative term (1/c)\nabla_perp c\cdot\nabla_\perp x
    
	  if((ix > 1) && (ix < (ngx-2)))
	    coef4 += mesh->g11[ix][jy] * ((*ccoef)[ix-2][jy] - 8.*(*ccoef)[ix-1][jy] + 8.*(*ccoef)[ix+1][jy] - (*ccoef)[ix+2][jy]) / (12.*dx[ix][jy]*((*ccoef)[ix][jy]));
	}

	// Put into matrix
	coef1 /= 12.* SQ(dx[ix][jy]);
	coef2 *= SQ(kwave);
	coef3 *= kwave / (12. * dx[ix][jy]);
	coef4 /= 12. * dx[ix][jy];
	coef5 *= kwave;

	A[ix][0] = dcomplex(    -coef1 +   coef4 ,     coef3 );
	A[ix][1] = dcomplex( 16.*coef1 - 8*coef4 , -8.*coef3 );
	A[ix][2] = dcomplex(-30.*coef1 - coef2 + coef6, coef5);
	A[ix][3] = dcomplex( 16.*coef1 + 8*coef4 ,  8.*coef3 );
	A[ix][4] = dcomplex(    -coef1 -   coef4 ,    -coef3 );
      }

      if(xbndry < 2) {
	// Use 2nd order near edges

	ix = 1;

	coef1=mesh->g11[ix][jy]/(SQ(dx[ix][jy]));
	coef2=mesh->g33[ix][jy];
	coef3= kwave * mesh->g13[ix][jy]/(2. * dx[ix][jy]);

	A[ix][0] = 0.0; // Should never be used
	A[ix][1] = dcomplex(coef1, -coef3);
	A[ix][2] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	A[ix][3] = dcomplex(coef1,  coef3);
	A[ix][4] = 0.0;

	ix = ncx-1;

	coef1=mesh->g11[ix][jy]/(SQ(dx[ix][jy]));
	coef2=mesh->g33[ix][jy];
	coef3= kwave * mesh->g13[ix][jy]/(2. * dx[ix][jy]);

	A[ix][0] = 0.0;
	A[ix][1] = dcomplex(coef1, -coef3);
	A[ix][2] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	A[ix][3] = dcomplex(coef1,  coef3);
	A[ix][4] = 0.0;  // Should never be used
      }

      // Boundary conditions

      for(ix=0;ix<xbndry;ix++) {
	// Set zero-value. Change to zero-gradient if needed
	bk1d[ix] = 0.0;
	bk1d[ncx-ix] = 0.0;

	A[ix][0] = A[ix][1] = A[ix][3] = A[ix][4] = 0.0;
	A[ix][2] = 1.0;

	A[ncx-ix][0] = A[ncx-ix][1] = A[ncx-ix][3] = A[ncx-ix][4] = 0.0;
	A[ncx-ix][2] = 1.0;
      }

      if(flags & INVERT_IN_SET) {
	// Set values of inner boundary from X
	for(ix=0;ix<xbndry;ix++)
	  bk1d[ix] = xk[ix][iz];
      }
      if(flags & INVERT_OUT_SET) {
	// Set values of outer boundary from X
	for(ix=0;ix<xbndry;ix++)
	  bk1d[ncx-ix] = xk[ncx-ix][iz];
      }

      if(iz == 0) {
	// DC
	
	// Inner boundary
	if(flags & INVERT_DC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  for (ix=0;ix<xbndry;ix++)
	    A[ix][3] = -1.0;
	}
	
	// Outer boundary
	if(flags & INVERT_DC_OUT_GRAD) {
	  // Zero gradient at outer boundary
	  for (ix=0;ix<xbndry;ix++)
	    A[ncx-ix][1] = -1.0;
	}
	
      }else {
	// AC
	
	// Inner boundarySQ(kwave)*coef2
	if(flags & INVERT_AC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  for (ix=0;ix<xbndry;ix++)
	    A[ix][3] = -1.0;
	}else if(flags & INVERT_AC_IN_LAP) {
	  // Enforce zero laplacian for 2nd and 4th-order
	  
	  ix = 1;
	  
	  coef1=mesh->g11[ix][jy]/(12.* SQ(dx[ix][jy]));
	
	  coef2=mesh->g33[ix][jy];
	
	  coef3= kwave * mesh->g13[ix][jy]/(2. * dx[ix][jy]);
	  
	  coef4 = 0.0;
	  if(a != (Field2D*) NULL)
	    coef4 = (*a)[ix][jy];
	  
	  // Combine 4th order at 1 with 2nd order at 0
	  A[1][0] = 0.0; // Not used
	  A[1][1] = dcomplex( (14. - SQ(dx[0][jy]*kwave)*mesh->g33[0][jy]/mesh->g11[0][jy])*coef1  ,  -coef3 );
	  A[1][2] = dcomplex(-29.*coef1 - SQ(kwave)*coef2 + coef4, 0.0);
	  A[1][3] = dcomplex( 16.*coef1  , coef3 );
	  A[1][4] = dcomplex(    -coef1  ,     0.0 );
	  
	  coef1=mesh->g11[ix][jy]/(SQ(dx[ix][jy]));
	  coef2=mesh->g33[ix][jy];
	  coef3= kwave * mesh->g13[ix][jy]/(2. * dx[ix][jy]);

	  // Use 2nd order at 1
	  A[0][0] = 0.0;  // Should never be used
	  A[0][1] = 0.0;
	  A[0][2] = dcomplex(coef1, -coef3);
	  A[0][3] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	  A[0][4] = dcomplex(coef1,  coef3);
	  
	}
	
	// Outer boundary
	if(flags & INVERT_AC_OUT_GRAD) {
	  // Zero gradient at outer boundary
	  for (ix=0;ix<xbndry;ix++)
	    A[ncx-ix][1] = -1.0;
	}else if(flags & INVERT_AC_OUT_LAP) {
	  // Enforce zero laplacian for 2nd and 4th-order
	  // NOTE: Currently ignoring XZ term and coef4 assumed zero on boundary
	  // FIX THIS IF IT WORKS

	  ix = ncx-1;
	  
	  coef1=mesh->g11[ix][jy]/(12.* SQ(dx[ix][jy]));
	
	  coef2=mesh->g33[ix][jy];
	
	  coef3= kwave * mesh->g13[ix][jy]/(2. * dx[ix][jy]);
	  
	  coef4 = 0.0;
	  if(a != (Field2D*) NULL)
	    coef4 = (*a)[ix][jy];
	  
	  // Combine 4th order at ncx-1 with 2nd order at ncx
	  A[ix][0] = dcomplex(    -coef1  ,     0.0 );
	  A[ix][1] = dcomplex( 16.*coef1  , -coef3 );
	  A[ix][2] = dcomplex(-29.*coef1 - SQ(kwave)*coef2 + coef4, 0.0);
	  A[ix][3] = dcomplex( (14. - SQ(dx[ncx][jy]*kwave)*mesh->g33[ncx][jy]/mesh->g11[ncx][jy])*coef1  ,  coef3 );
	  A[ix][4] = 0.0; // Not used
	  
	  coef1=mesh->g11[ix][jy]/(SQ(dx[ix][jy]));
	  coef2=mesh->g33[ix][jy];
	  coef3= kwave * mesh->g13[ix][jy]/(2. * dx[ix][jy]);

	  // Use 2nd order at ncx - 1
	  A[ncx][0] = dcomplex(coef1, -coef3);
	  A[ncx][1] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	  A[ncx][2] = dcomplex(coef1,  coef3);
	  A[ncx][3] = 0.0;  // Should never be used
	  A[ncx][4] = 0.0;
	}
      }

      // Perform inversion
      cband_solve(A, ngx, 2, 2, bk1d);
      // Fill xk
      
      for (ix=0; ix<=ncx; ix++)
	xk[ix][iz]=bk1d[ix];
      
    }
  }else {
    // Use tridiagonal system in x - 2nd order
    
    static dcomplex *avec = (dcomplex*) NULL, *bvec, *cvec;
    
    if(avec == (dcomplex*) NULL) {
      avec = new dcomplex[ngx];
      bvec = new dcomplex[ngx];
      cvec = new dcomplex[ngx];
    }

    for(iz=0;iz<=ncz/2;iz++) {
      // solve differential equation in x

      // set bk1d
      for(ix=0;ix<=ncx;ix++)
	bk1d[ix] = bk[ix][iz];

      ///////// PERFORM INVERSION /////////
      
      if (iz>laplace_maxmode) flt=0.0; else flt=1.0;
      
      for(ix=xbndry;ix<=ncx-xbndry;ix++) {
	bk1d[ix] *= flt;
	
	laplace_tridag_coefs(ix, jy, iz, avec[ix], bvec[ix], cvec[ix], ccoef);
	
	if(a != (Field2D*) NULL)
	  bvec[ix] += (*a)[ix][jy];
      }
      
      // Set boundary conditions
      
      if(iz == 0) {
	// DC
	
	// Inner boundary
	if(flags & INVERT_DC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  
	  if((flags & INVERT_IN_SYM) && (xbndry > 1) && BoundaryOnCell) {
	    // Use symmetric boundary to set zero-gradient
	    
	    for (ix=0;ix<xbndry-1;ix++) {
	      avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      bk1d[ix]=0.0;
	    }
	    // Symmetric on last point
	    avec[xbndry-1] = 1.0; bvec[xbndry-1] = 0.0; cvec[xbndry-1] = -1.0;
	    bk1d[xbndry-1] = 0.0;
	  }else {
	    for (ix=0;ix<xbndry;ix++){
	      avec[ix]=dcomplex(0.0,0.0);
	      bvec[ix]=dcomplex(1.,0.); cvec[ix]=dcomplex(-1.,0.);bk1d[ix]=dcomplex(0.0,0.0);
	    }
	  }
	}else if(flags & INVERT_IN_SET) {
	  for(ix=0;ix<xbndry;ix++) {
	    avec[ix] = 0.0;
	    bvec[ix] = 1.0;
	    cvec[ix] = 0.0;
	    bk1d[ix] = xk[ix][iz];
	  }
	}else {
	  // Zero value at inner boundary
	  if(flags & INVERT_IN_SYM) {
	    // Use anti-symmetric boundary to set zero-value
	    
	    // Zero-gradient for first point(s)
	    for(ix=0;ix<xbndry-1;ix++) {
	      avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      bk1d[ix]=0.0;
	    }

	    if(BoundaryOnCell) {
	      // Antisymmetric about boundary on cell
	      avec[xbndry-1]=1.0; bvec[xbndry-1]=0.0; cvec[xbndry-1]= 1.0;
	      bk1d[xbndry-1]=0.0;
	    }else { 
	      // Antisymmetric across boundary between cells
	      avec[xbndry-1]=0.0; bvec[xbndry-1]=1.0; cvec[xbndry-1]= 1.0;
	      bk1d[xbndry-1]=0.0;
	    }
	    
	  }else {
	    for (ix=0;ix<xbndry;ix++){
	      avec[ix]=dcomplex(0.,0.);
	      bvec[ix]=dcomplex(1.,0.);cvec[ix]=dcomplex(0.,0.);bk1d[ix]=dcomplex(0.0,0.0);
	    }
	  }
	}
	
	// Outer boundary
	if(flags & INVERT_DC_OUT_GRAD) {
	  // Zero gradient at outer boundary

	  if((flags & INVERT_OUT_SYM) && (xbndry > 1) && BoundaryOnCell) {
	    // Use symmetric boundary to set zero-gradient
	    
	    for (ix=0;ix<xbndry-1;ix++) {
	      avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      bk1d[ncx-ix]=0.0;
	    }
	    // Symmetric on last point
	    ix = xbndry-1;
	    avec[ncx-ix] = 1.0; bvec[ncx-ix] = 0.0; cvec[ncx-ix] = -1.0;
	    bk1d[ncx-ix] = 0.0;
	    
	  }else {
	    for (ix=0;ix<xbndry;ix++){
	      cvec[ncx-ix]=dcomplex(0.,0.);
	      bvec[ncx-ix]=dcomplex(1.,0.);avec[ncx-ix]=dcomplex(-1.,0.);bk1d[ncx-ix]=dcomplex(0.0,0.0);
	    }
	  }
	}else if(flags & INVERT_OUT_SET) {
	  // Setting the values in the outer boundary
	  for(ix=0;ix<xbndry;ix++) {
	    avec[ncx-ix] = 0.0;
	    bvec[ncx-ix] = 1.0;
	    cvec[ncx-ix] = 0.0;
	    bk1d[ncx-ix] = xk[ncx-ix][iz];
	  }
	}else {
	  // Zero value at outer boundary
	  if(flags & INVERT_OUT_SYM) {
	    // Use anti-symmetric boundary to set zero-value
	    
	    // Zero-gradient for first point(s)
	    for(ix=0;ix<xbndry-1;ix++) {
	      avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      bk1d[ncx-ix]=0.0;
	    }
	    ix = xbndry-1;
	    if(BoundaryOnCell) {
	      // Antisymmetric about boundary on cell
	      avec[ncx-ix]=1.0; bvec[ncx-ix]=0.0; cvec[ncx-ix]= 1.0;
	      bk1d[ncx-ix]=0.0;
	    }else { 
	      // Antisymmetric across boundary between cells
	      avec[ncx-ix]=1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      bk1d[ncx-ix]=0.0;
	    }
	  }else {
	    for (ix=0;ix<xbndry;ix++){
	      cvec[ncx-ix]=dcomplex(0.,0.);
	      bvec[ncx-ix]=dcomplex(1.,0.);avec[ncx-ix]=dcomplex(0.,0.);bk1d[ncx-ix]=dcomplex(0.0,0.0);
	    }
	  }
	}
      }else {
	// AC
	
	// Inner boundary
	if(flags & INVERT_AC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  
	  if((flags & INVERT_IN_SYM) && (xbndry > 1) && BoundaryOnCell) {
	    // Use symmetric boundary to set zero-gradient
	    
	    for (ix=0;ix<xbndry-1;ix++) {
	      avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      bk1d[ix]=0.0;
	    }
	    // Symmetric on last point
	    avec[xbndry-1] = 1.0; bvec[xbndry-1] = 0.0; cvec[xbndry-1] = -1.0;
	    bk1d[xbndry-1] = 0.0;
	  }else {
	    for (ix=0;ix<xbndry;ix++){
	      avec[ix]=dcomplex(0.,0.);
	      bvec[ix]=dcomplex(1.,0.);cvec[ix]=dcomplex(-1.,0.);bk1d[ix]=dcomplex(0.0,0.0);
	    }
	  }
	}else if(flags & INVERT_IN_SET) {
	  // Setting the values in the boundary
	  for(ix=0;ix<xbndry;ix++) {
	    avec[ix] = 0.0;
	    bvec[ix] = 1.0;
	    cvec[ix] = 0.0;
	    bk1d[ix] = xk[ix][iz];
	  }
	}else if(flags & INVERT_AC_IN_LAP) {
	  // Use decaying zero-Laplacian solution in the boundary
	  real kwave=iz*2.0*PI/zlength; // wave number is 1/[rad]
	  for (ix=0;ix<xbndry;ix++) {
	    avec[ix] = 0.0;
	    bvec[ix] = -1.0;
	    cvec[ix] = exp(-1.0*sqrt(mesh->g33[ix][jy]/mesh->g11[ix][jy])*kwave*dx[ix][jy]);
	    bk1d[ix] = 0.0;
	  }
	}else {
	  // Zero value at inner boundary

	  if(flags & INVERT_IN_SYM) {
	    // Use anti-symmetric boundary to set zero-value
	    
	    // Zero-gradient for first point(s)
	    for(ix=0;ix<xbndry-1;ix++) {
	      avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      bk1d[ix]=0.0;
	    }

	    if(BoundaryOnCell) {
	      // Antisymmetric about boundary on cell
	      avec[xbndry-1]=1.0; bvec[xbndry-1]=0.0; cvec[xbndry-1]= 1.0;
	      bk1d[xbndry-1]=0.0;
	    }else { 
	      // Antisymmetric across boundary between cells
	      avec[xbndry-1]=0.0; bvec[xbndry-1]=1.0; cvec[xbndry-1]= 1.0;
	      bk1d[xbndry-1]=0.0;
	    }
	    
	  }else {
	    for (ix=0;ix<xbndry;ix++){
	      avec[ix]=dcomplex(0.,0.);
	      bvec[ix]=dcomplex(1.,0.);cvec[ix]=dcomplex(0.,0.);bk1d[ix]=dcomplex(0.0,0.0);
	    }
	  }
	}
	
	// Outer boundary
	if(flags & INVERT_AC_OUT_GRAD) {
	  // Zero gradient at outer boundary
	  
	  if((flags & INVERT_OUT_SYM) && (xbndry > 1) && BoundaryOnCell) {
	    // Use symmetric boundary to set zero-gradient
	    
	    for (ix=0;ix<xbndry-1;ix++) {
	      avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      bk1d[ncx-ix]=0.0;
	    }
	    // Symmetric on last point
	    ix = xbndry-1;
	    avec[ncx-ix] = 1.0; bvec[ncx-ix] = 0.0; cvec[ncx-ix] = -1.0;
	    bk1d[ncx-ix] = 0.0;
	    
	  }else {
	    for (ix=0;ix<xbndry;ix++){
	      cvec[ncx-ix]=dcomplex(0.,0.);
	      bvec[ncx-ix]=dcomplex(1.,0.);avec[ncx-ix]=dcomplex(-1.,0.);bk1d[ncx-ix]=dcomplex(0.0,0.0);
	    }
	  }
	}else if(flags & INVERT_AC_OUT_LAP) {
	  // Use decaying zero-Laplacian solution in the boundary
	  real kwave=iz*2.0*PI/zlength; // wave number is 1/[rad]
	  for (ix=0;ix<xbndry;ix++) {
	    avec[ncx-ix] = exp(-1.0*sqrt(mesh->g33[ncx-ix][jy]/mesh->g11[ncx-ix][jy])*kwave*dx[ncx-ix][jy]);;
	    bvec[ncx-ix] = -1.0;
	    cvec[ncx-ix] = 0.0;
	    bk1d[ncx-ix] = 0.0;
	  }
	}else if(flags & INVERT_OUT_SET) {
	  // Setting the values in the outer boundary
	  for(ix=0;ix<xbndry;ix++) {
	    avec[ncx-ix] = 0.0;
	    bvec[ncx-ix] = 1.0;
	    cvec[ncx-ix] = 0.0;
	    bk1d[ncx-ix] = xk[ncx-ix][iz];
	  }
	}else {
	  // Zero value at outer boundary

	  if(flags & INVERT_OUT_SYM) {
	    // Use anti-symmetric boundary to set zero-value
	    
	    // Zero-gradient for first point(s)
	    for(ix=0;ix<xbndry-1;ix++) {
	      avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      bk1d[ncx-ix]=0.0;
	    }
	    ix = xbndry-1;
	    if(BoundaryOnCell) {
	      // Antisymmetric about boundary on cell
	      avec[ncx-ix]=1.0; bvec[ncx-ix]=0.0; cvec[ncx-ix]= 1.0;
	      bk1d[ncx-ix]=0.0;
	    }else { 
	      // Antisymmetric across boundary between cells
	      avec[ncx-ix]=1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      bk1d[ncx-ix]=0.0;
	    }
	  }else {
	    for (ix=0;ix<xbndry;ix++){
	      cvec[ncx-ix]=dcomplex(0.,0.);
	      bvec[ncx-ix]=dcomplex(1.,0.);avec[ncx-ix]=dcomplex(0.,0.);bk1d[ncx-ix]=dcomplex(0.0,0.0);
	    }
	  }
	}
      }
      
      // Call tridiagonal solver
      tridag(avec, bvec, cvec, bk1d, xk1d, ngx);

      if((flags & INVERT_IN_SYM) && (xbndry > 1)) {
	// (Anti-)symmetry on inner boundary. Nothing to do if only one boundary cell
	int xloc = 2*xbndry;
	if(!BoundaryOnCell)
	  xloc--;
	
	if( ((iz == 0) && (flags & INVERT_DC_IN_GRAD)) || ((iz != 0) && (flags & INVERT_AC_IN_GRAD)) ) {
	  // Inner gradient zero - symmetric
	  for(ix=0;ix<xbndry-1;ix++)
	    xk1d[ix] = xk1d[xloc-ix];
	}else {
	  // Inner value zero - antisymmetric
	  for(ix=0;ix<xbndry-1;ix++)
	    xk1d[ix] = -xk1d[xloc-ix];
	}
      }
      if((flags & INVERT_OUT_SYM) && (xbndry > 1)) {
	// (Anti-)symmetry on outer boundary. Nothing to do if only one boundary cell
	
	int xloc =  ngx - 2*xbndry;
	if(BoundaryOnCell)
	  xloc--;
	
	if( ((iz == 0) && (flags & INVERT_DC_IN_GRAD)) || ((iz != 0) && (flags & INVERT_AC_IN_GRAD)) ) {
	  // Outer gradient zero - symmetric
	  for(ix=0;ix<xbndry-1;ix++)
	    xk1d[ncx-ix] = xk1d[xloc + ix];
	}else {
	  // Outer value zero - antisymmetric
	  for(ix=0;ix<xbndry-1;ix++)
	    xk1d[ncx-ix] = -xk1d[xloc + ix];
	}
      }
      
      // Fill xk
      
      for (ix=0; ix<=ncx; ix++){
	xk[ix][iz]=xk1d[ix];
      }
    }
  }

  // Done inversion, transform back

  for(ix=0; ix<=ncx; ix++){
    
    if(flags & INVERT_ZERO_DC)
      xk[ix][0] = 0.0;

    ZFFT_rev(xk[ix], zShift[ix][jy], x[ix]);
    
    x[ix][ncz] = x[ix][0]; // enforce periodicity
  }

  return 0;
}

/**********************************************************************************
 *                           PARALLEL CODE - COMMON
 **********************************************************************************/

/// Sets the coefficients for parallel tridiagonal matrix inversion
/*!
 * Uses the laplace_tridag_coefs routine above to fill a matrix [kz][ix] of coefficients
 */
void par_tridag_matrix(dcomplex **avec, dcomplex **bvec, dcomplex **cvec,
		       dcomplex **bk, int jy, int flags, const Field2D *a = NULL, const Field2D *ccoef=NULL)
{
  int ix, kz;
  
  int xbndry = MXG;
  if(flags & INVERT_BNDRY_ONE)
    xbndry = 1;
  
  for(kz = 0; kz <= laplace_maxmode; kz++) {
    
    // Entire domain. Change to boundaries later

    for(ix=0;ix<=ncx;ix++) {

      laplace_tridag_coefs(ix, jy, kz, avec[kz][ix], bvec[kz][ix], cvec[kz][ix], ccoef);
	
      if(a != (Field2D*) NULL)
	bvec[kz][ix] += (*a)[ix][jy];
    }

    // Boundary conditions

    if(PE_XIND == 0) {
      // INNER BOUNDARY ON THIS PROCESSOR
      
      if(kz == 0) {
	// DC
	
	if(flags & INVERT_DC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  for (ix=0;ix<xbndry;ix++){
	    avec[kz][ix]=dcomplex(0.0,0.0);
	    bvec[kz][ix]=dcomplex(1.,0.);
	    cvec[kz][ix]=dcomplex(-1.,0.);
	    bk[kz][ix]=dcomplex(0.0,0.0);
	  }
	}else {
	  // Zero value at inner boundary
	  for (ix=0;ix<xbndry;ix++){
	    avec[kz][ix]=dcomplex(0.,0.);
	    bvec[kz][ix]=dcomplex(1.,0.);
	    cvec[kz][ix]=dcomplex(0.,0.);
	    bk[kz][ix]=dcomplex(0.0,0.0);
	  }
	}
	
      }else {
	// AC
	
	if(flags & INVERT_AC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  for (ix=0;ix<xbndry;ix++){
	    avec[kz][ix]=dcomplex(0.,0.);
	    bvec[kz][ix]=dcomplex(1.,0.);
	    cvec[kz][ix]=dcomplex(-1.,0.);
	    bk[kz][ix]=dcomplex(0.0,0.0);
	  }
	}else if(flags & INVERT_AC_IN_LAP) {
	  // Use decaying zero-Laplacian solution in the boundary
	  real kwave=kz*2.0*PI/zlength; // wave number is 1/[rad]
	  for (ix=0;ix<xbndry;ix++) {
	    avec[kz][ix] = 0.0;
	    bvec[kz][ix] = 1.0;
	    cvec[kz][ix] = -exp(-1.0*sqrt(mesh->g33[ix][jy]/mesh->g11[ix][jy])*kwave*dx[ix][jy]);
	    bk[kz][ix] = 0.0;
	  }
	}else {
	  // Zero value at inner boundary
	  for (ix=0;ix<xbndry;ix++){
	    avec[kz][ix]=dcomplex(0.,0.);
	    bvec[kz][ix]=dcomplex(1.,0.);
	    cvec[kz][ix]=dcomplex(0.,0.);
	    bk[kz][ix]=dcomplex(0.0,0.0);
	  }
	}
      }
    }else if(PE_XIND == (NXPE - 1)) {
      // OUTER BOUNDARY
      
      if(kz == 0) {
	// DC
	
	if(flags & INVERT_DC_OUT_GRAD) {
	  // Zero gradient at outer boundary
	  for (ix=0;ix<xbndry;ix++){
	    cvec[kz][ncx-ix]=dcomplex(0.,0.);
	    bvec[kz][ncx-ix]=dcomplex(1.,0.);
	    avec[kz][ncx-ix]=dcomplex(-1.,0.);
	    bk[kz][ncx-ix]=dcomplex(0.0,0.0);
	  }
	}else {
	  // Zero value at outer boundary
	  for (ix=0;ix<xbndry;ix++){
	    cvec[kz][ncx-ix]=dcomplex(0.,0.);
	    bvec[kz][ncx-ix]=dcomplex(1.,0.);
	    avec[kz][ncx-ix]=dcomplex(0.,0.);
	    bk[kz][ncx-ix]=dcomplex(0.0,0.0);
	  }
	}
      }else {
	// AC
	
	if(flags & INVERT_AC_OUT_GRAD) {
	  // Zero gradient at outer boundary
	  for (ix=0;ix<xbndry;ix++){
	    cvec[kz][ncx-ix]=dcomplex(0.,0.);
	    bvec[kz][ncx-ix]=dcomplex(1.,0.);
	    avec[kz][ncx-ix]=dcomplex(-1.,0.);
	    bk[kz][ncx-ix]=dcomplex(0.0,0.0);
	  }
	}else if(flags & INVERT_AC_OUT_LAP) {
	  // Use decaying zero-Laplacian solution in the boundary
	  real kwave=kz*2.0*PI/zlength; // wave number is 1/[rad]
	  for (ix=0;ix<xbndry;ix++) {
	    avec[kz][ncx-ix] = -exp(-1.0*sqrt(mesh->g33[ncx-ix][jy]/mesh->g11[ncx-ix][jy])*kwave*dx[ncx-ix][jy]);;
	    bvec[kz][ncx-ix] = 1.0;
	    cvec[kz][ncx-ix] = 0.0;
	    bk[kz][ncx-ix] = 0.0;
	  }
	}else {
	  // Zero value at outer boundary
	  for (ix=0;ix<xbndry;ix++){
	    cvec[kz][ncx-ix]=dcomplex(0.,0.);
	    bvec[kz][ncx-ix]=dcomplex(1.,0.);
	    avec[kz][ncx-ix]=dcomplex(0.,0.);
	    bk[kz][ncx-ix]=dcomplex(0.0,0.0);
	  }
	}
      }
    }
  }
}

/**********************************************************************************
 *                           PARALLEL CODE - SIMPLE ALGORITHM
 * 
 * I'm just calling this Simple Parallel Tridag. Naive parallelisation of
 * the serial code. For use as a reference case.
 * 
 * Overlap calculation / communication of poloidal slices to achieve some
 * parallelism.
 **********************************************************************************/

/// Data structure for SPT algorithm
typedef struct {
  int jy; ///< Y index

  dcomplex **bk;  ///< b vector in Fourier space
  dcomplex **xk;

  dcomplex **gam;
  
  dcomplex **avec, **bvec, **cvec; ///< Diagonal bands of matrix

  int proc; // Which processor has this reached?
  int dir;  // Which direction is it going?
  
  MPI_Request send_req, recv_req;

  real *buffer;
}SPT_data;


/// This is the first half of the Thomas algorithm for parallel calculations
/*!
 * Two complex quantities have to be propagated between processors: bet and u[-1].
 * This routine takes bet and um from the last processor (if start == false),
 * and returns the values to be passed to the next processor in the same variables.
 *
 * @param[in]  a    Vector of matrix coefficients (Left of diagonal)
 * @param[in]  b    Vector of matrix coefficients (Diagonal)
 * @param[in]  c    Vector of matrix coefficients (Right of diagonal)
 * @param[in]  r    RHS vector
 * @param[in]  u    Result vector (Au = r)
 * @param[in]  n    Size of the matrix
 * @param[out] gam  Intermediate values used for backsolve stage
 * @param[inout] bet
 * @param[inout] um
 * @param[in] start
 */
void spt_tridag_forward(dcomplex *a, dcomplex *b, dcomplex *c,
			dcomplex *r, dcomplex *u, int n,
			dcomplex *gam,
			dcomplex &bet, dcomplex &um, bool start=false)
{
  int j;
  
  if(start) {
    bet = b[0];
    u[0] = r[0] / bet;
  }else {
    //output.write("um = %e,%e\n", um.Real(), um.Imag());
    gam[0] = c[-1] / bet; // NOTE: ASSUMES C NOT CHANGING
    bet = b[0] - a[0]*gam[0];
    u[0] = (r[0]-a[0]*um)/bet;
  }
  
  for(j=1;j<n;j++) {
    gam[j] = c[j-1]/bet;
    bet = b[j]-a[j]*gam[j];
    if(bet == 0.0) {
      printf("Tridag: Zero pivot\n");
      exit(1);
    }
    u[j] = (r[j]-a[j]*u[j-1])/bet;
  }

  um = u[n-1];
}

/// Second (backsolve) part of the Thomas algorithm
/*!
 * @param[inout] u    Result to be solved (Au = r)
 * @param[in]    n    Size of the problem
 * @param[in]    gam  Intermediate values produced by the forward part
 * @param[inout] gp   gam from the processor PE_XIND + 1, and returned to PE_XIND - 1
 * @param[inout] up   u from processor PE_XIND + 1, and returned to PE_XIND - 1
 */
void spt_tridag_back(dcomplex *u, int n,
		     dcomplex *gam, dcomplex &gp, dcomplex &up)
{
  int j;

  u[n-1] = u[n-1] - gp*up;

  for(j=n-2;j>=0;j--) {
    u[j] = u[j]-gam[j+1]*u[j+1];
  }
  gp = gam[0];
  up = u[0];
}

const int SPT_DATA = 1123; ///< 'magic' number for SPT MPI messages

/// Simple parallelisation of the Thomas tridiagonal solver algorithm (serial code)
/*!
 * This is a reference code which performs the same operations as the serial code.
 * To invert a single XZ slice (FieldPerp object), data must pass from the innermost
 * processor (PE_XIND = 0) to the outermost (PE_XIND = NXPE-1) and back again.
 *
 * Some parallelism is achieved by running several inversions simultaneously, so while
 * processor #1 is inverting Y=0, processor #0 is starting on Y=1. This works ok as long
 * as the number of slices to be inverted is greater than the number of X processors (MYSUB > NXPE).
 * If MYSUB < NXPE then not all processors can be busy at once, and so efficiency will fall sharply.
 *
 * @param[in]    b      RHS values (Ax = b)
 * @param[in]    flags  Inversion settings (see boundary.h for values)
 * @param[in]    a      This is a 2D matrix which allows solution of A = Delp2 + a
 * @param[out]   data   Structure containing data needed for second half of inversion
 * @param[in]    ccoef  Optional coefficient for first-order derivative
 */
int invert_spt_start(const FieldPerp &b, int flags, const Field2D *a, SPT_data &data, const Field2D *ccoef = NULL)
{
  if(NXPE == 1) {
    output.write("Error: SPT method only works for NXPE > 1\n");
    return 1;
  }

  data.send_req = data.recv_req = MPI_REQUEST_NULL;

  data.jy = b.getIndex();

  if(data.bk == NULL) {
    /// Allocate memory
    
    // RHS vector
    data.bk = cmatrix(laplace_maxmode + 1, ngx);
    data.xk = cmatrix(laplace_maxmode + 1, ngx);
    
    data.gam = cmatrix(laplace_maxmode + 1, ngx);

    // Matrix to be solved
    data.avec = cmatrix(laplace_maxmode + 1, ngx);
    data.bvec = cmatrix(laplace_maxmode + 1, ngx);
    data.cvec = cmatrix(laplace_maxmode + 1, ngx);
    
    data.buffer  = new real[4*(laplace_maxmode + 1)];
  }

  /// Take FFTs of data
  static dcomplex *bk1d = NULL; ///< 1D in Z for taking FFTs
  int ix, kz;

  if(bk1d == NULL)
    bk1d = new dcomplex[ncz/2 + 1];

  for(ix=0; ix < ngx; ix++) {
    ZFFT(b[ix], zShift[ix][data.jy], bk1d);
    for(kz = 0; kz <= laplace_maxmode; kz++)
      data.bk[kz][ix] = bk1d[kz];
  }
  
  /// Set matrix elements
  par_tridag_matrix(data.avec, data.bvec, data.cvec,
		    data.bk, data.jy, flags, a, ccoef);

  data.proc = 0; //< Starts at processor 0
  data.dir = 1;
  
  if(PE_XIND == 0) {
    dcomplex bet, u0;
    for(kz = 0; kz <= laplace_maxmode; kz++) {
      // Start tridiagonal solve
      spt_tridag_forward(data.avec[kz], data.bvec[kz], data.cvec[kz],
			 data.bk[kz], data.xk[kz], MXG+MXSUB,
			 data.gam[kz],
			 bet, u0, true);
      // Load intermediate values into buffers
      data.buffer[4*kz]     = bet.Real();
      data.buffer[4*kz + 1] = bet.Imag();
      data.buffer[4*kz + 2] = u0.Real();
      data.buffer[4*kz + 3] = u0.Imag();
    }
    
    // Send data

    if(invert_async_send) {
      MPI_Isend(data.buffer, 
		4*(laplace_maxmode+1),
		PVEC_REAL_MPI_TYPE,
		PROC_NUM(1, PE_YIND),
		SPT_DATA,
		MPI_COMM_WORLD,
		&data.send_req);
    }else
      MPI_Send(data.buffer, 
	       4*(laplace_maxmode+1),
	       PVEC_REAL_MPI_TYPE,
	       PROC_NUM(1, PE_YIND),
	       SPT_DATA,
	       MPI_COMM_WORLD);
    
  }else if(PE_XIND == 1) {
    // Post a receive
    
    MPI_Irecv(data.buffer,
	      4*(laplace_maxmode+1),
	      PVEC_REAL_MPI_TYPE,
	      PROC_NUM(0, PE_YIND),
	      SPT_DATA,
	      MPI_COMM_WORLD,
	      &data.recv_req);
  }
  
  data.proc++; // Now moved onto the next processor
  if(NXPE == 2)	
    data.dir = -1; // Special case. Otherwise reversal handled in spt_continue

  return 0;
}

/// Shifts the parallelised Thomas algorithm along one processor.
/*!
  Returns non-zero when the calculation is complete.

  @param[inout] data  Structure which keeps track of the calculation
*/
int invert_spt_continue(SPT_data &data)
{
  MPI_Status status;

  if(data.proc < 0) // Already finished
    return 1;
  
  if(PE_XIND == data.proc) {
    /// This processor's turn to do inversion

    // Wait for data to arrive
    MPI_Wait(&data.recv_req, &status);

    if(PE_XIND == (NXPE - 1)) {
      // Last processor, turn-around
      
      dcomplex bet, u0;
      dcomplex gp, up;
      for(int kz = 0; kz <= laplace_maxmode; kz++) {
	bet = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	u0 = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);
	spt_tridag_forward(data.avec[kz]+MXG, data.bvec[kz]+MXG, data.cvec[kz]+MXG,
			   data.bk[kz]+MXG, data.xk[kz]+MXG, MXG+MXSUB,
			   data.gam[kz]+MXG,
			   bet, u0);
	
	// Back-substitute
	gp = 0.0;
	up = 0.0;
	spt_tridag_back(data.xk[kz]+MXG, MXG+MXSUB, data.gam[kz]+MXG, gp, up);
	data.buffer[4*kz]     = gp.Real();
	data.buffer[4*kz + 1] = gp.Imag();
	data.buffer[4*kz + 2] = up.Real();
	data.buffer[4*kz + 3] = up.Imag();
      }

    }else if(data.dir > 0) {
      // In the middle of X, forward direction

      dcomplex bet, u0;
      for(int kz = 0; kz <= laplace_maxmode; kz++) {
	
	bet = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	u0 = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);
	spt_tridag_forward(data.avec[kz]+MXG, data.bvec[kz]+MXG, data.cvec[kz]+MXG,
			   data.bk[kz]+MXG, data.xk[kz]+MXG, MXSUB,
			   data.gam[kz]+MXG,
			   bet, u0);
	// Load intermediate values into buffers
	data.buffer[4*kz]     = bet.Real();
	data.buffer[4*kz + 1] = bet.Imag();
	data.buffer[4*kz + 2] = u0.Real();
	data.buffer[4*kz + 3] = u0.Imag();
      }
      
    }else if(PE_XIND == 0) {
      // Back to the start
      
      dcomplex gp, up;
      for(int kz = 0; kz <= laplace_maxmode; kz++) {
	gp = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	up = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);

	spt_tridag_back(data.xk[kz], MXG+MXSUB, data.gam[kz], gp, up);
      }

    }else {
      // Middle of X, back-substitution stage

      dcomplex gp, up;
      for(int kz = 0; kz <= laplace_maxmode; kz++) {
	gp = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	up = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);

	spt_tridag_back(data.xk[kz]+MXG, MXSUB, data.gam[kz]+MXG, gp, up);
	
	data.buffer[4*kz]     = gp.Real();
	data.buffer[4*kz + 1] = gp.Imag();
	data.buffer[4*kz + 2] = up.Real();
	data.buffer[4*kz + 3] = up.Imag();
      }
    }

    if(PE_XIND != 0) { // If not finished yet
       /// Send data

      if(invert_async_send) {
	if(data.send_req != MPI_REQUEST_NULL) // Check the last send finished
	  MPI_Wait(&data.send_req, &status);
	
	MPI_Isend(data.buffer, 
		  4*(laplace_maxmode+1),
		  PVEC_REAL_MPI_TYPE,
		  PROC_NUM(data.proc + data.dir, PE_YIND),
		  SPT_DATA,
		  MPI_COMM_WORLD,
		  &data.send_req);
      }else
	MPI_Send(data.buffer, 
		 4*(laplace_maxmode+1),
		 PVEC_REAL_MPI_TYPE,
		 PROC_NUM(data.proc + data.dir, PE_YIND),
		 SPT_DATA,
		 MPI_COMM_WORLD);
    }

  }else if(PE_XIND == data.proc + data.dir) {
    // This processor is next, post receive
    
    if(invert_async_send && (data.send_req != MPI_REQUEST_NULL)) { 
      // Check the last send finished
      MPI_Wait(&data.send_req, &status);
      data.send_req = MPI_REQUEST_NULL;
    }
    
    MPI_Irecv(data.buffer,
	      4*(laplace_maxmode+1),
	      PVEC_REAL_MPI_TYPE,
	      PROC_NUM(data.proc, PE_YIND),
	      SPT_DATA,
	      MPI_COMM_WORLD,
	      &data.recv_req);
  }
  
  data.proc += data.dir;
  
  if(data.proc == NXPE-1)
    data.dir = -1; // Reverses direction at the end

  return 0;
}

/// Finishes the parallelised Thomas algorithm
/*!
  @param[inout] data   Structure keeping track of calculation
  @param[in]    flags  Inversion flags (same as passed to invert_spt_start)
  @param[out]   x      The result
*/
void invert_spt_finish(SPT_data &data, int flags, FieldPerp &x)
{
  int ix, kz;
  MPI_Status status;
  
  x.Allocate();
  x.setIndex(data.jy);

  // Make sure calculation has finished
  while(invert_spt_continue(data) == 0) {}

   // Have result in Fourier space. Convert back to real space

  static dcomplex *xk1d = NULL; ///< 1D in Z for taking FFTs

  if(xk1d == NULL) {
    xk1d = new dcomplex[ncz/2 + 1];
    for(kz=0;kz<=ncz/2;kz++)
      xk1d[kz] = 0.0;
  }

  // Make sure all comms finished (necessary to free memory)
  if(data.send_req != MPI_REQUEST_NULL) {
    MPI_Wait(&data.send_req, &status);
    data.send_req = MPI_REQUEST_NULL;
  }

  //output.write("xk[100][%d][1] = %e,%e\n",data.jy,data.xk[1][100].Real(), data.xk[1][100].Imag());

  for(ix=0; ix<=ncx; ix++){
    
    for(kz = 0; kz<= laplace_maxmode; kz++) {
      xk1d[kz] = data.xk[kz][ix];
    }

    if(flags & INVERT_ZERO_DC)
      xk1d[0] = 0.0;

    ZFFT_rev(xk1d, zShift[ix][data.jy], x[ix]);
    
    x[ix][ncz] = x[ix][0]; // enforce periodicity
  }

  if(PE_XIND != 0) {
    // Set left boundary to zero (Prevent unassigned values in corners)
    for(ix=0; ix<MXG; ix++){
      for(kz=0;kz<ngz;kz++)
	x[ix][kz] = 0.0;
    }
  }
  if(PE_XIND != (NXPE-1)) {
    // Same for right boundary
    for(ix=ngx-MXG; ix<ngx; ix++){
      for(kz=0;kz<ngz;kz++)
	x[ix][kz] = 0.0;
    }
  }
}

/**********************************************************************************
 *                           PARALLEL CODE - PDD ALGORITHM
 * 
 * This code uses the Parallel Diagonally Dominant algorithm. This is very efficient
 * (constant number of communications), but achieves this by neglecting "small"
 * corrections. For ELM simulations these seem to be non-negligable, hence:
 *
 * CHECK IF THIS ALGORITHM PRODUCES REASONABLE RESULTS FOR YOUR PROBLEM
 **********************************************************************************/

const int PDD_COMM_XV = 123; // First message tag
const int PDD_COMM_Y = 456;  // Second tag

/// Data structure for PDD algorithm
typedef struct {
  dcomplex **bk;  ///< b vector in Fourier space

  dcomplex **avec, **bvec, **cvec; ///< Diagonal bands of matrix
  
  int jy; ///< Y index
  
  dcomplex **xk;
  dcomplex **v, **w;

  real *snd; // send buffer
  real *rcv; // receive buffer
  
  MPI_Request snd_req, rcv_req; // Send and receive requests

  dcomplex *y2i;
}PDD_data;

/// Laplacian inversion using Parallel Diagonal Dominant (PDD) method
/*!
 *
 * July 2008: Adapted from serial version to run in parallel (split in X) for tridiagonal system
 * i.e. no 4th order inversion yet.
 *
 * \note This code stores intermediate results and takes significantly more memory than
 * the serial version. This can be balanced against communication time i.e. faster communications
 * can allow less memory use.
 *
 * @param[in] data  Internal data used for multiple calls in parallel mode
 * @param[in] stage Which stage of the inversion, used to overlap calculation and communications.
 */
int invert_pdd_start(const FieldPerp &b, int flags, const Field2D *a, PDD_data &data, const Field2D *ccoef = NULL)
{
  int ix, kz;
  
  data.jy = b.getIndex();

  if(NXPE == 1) {
    output.write("Error: PDD method only works for NXPE > 1\n");
    return 1;
  }

  if(data.bk == NULL) {
    // Need to allocate working memory
    
    // RHS vector
    data.bk = cmatrix(laplace_maxmode + 1, ngx);
    
    // Matrix to be solved
    data.avec = cmatrix(laplace_maxmode + 1, ngx);
    data.bvec = cmatrix(laplace_maxmode + 1, ngx);
    data.cvec = cmatrix(laplace_maxmode + 1, ngx);
    
    // Working vectors
    data.v = cmatrix(laplace_maxmode + 1, ngx);
    data.w = cmatrix(laplace_maxmode + 1, ngx);

    // Result
    data.xk = cmatrix(laplace_maxmode + 1, ngx);

    // Communication buffers. Space for 2 complex values for each kz
    data.snd = new real[4*(laplace_maxmode+1)];
    data.rcv = new real[4*(laplace_maxmode+1)];

    data.y2i = new dcomplex[laplace_maxmode + 1];
  }

  /// Take FFTs of data
  static dcomplex *bk1d = NULL; ///< 1D in Z for taking FFTs

  if(bk1d == NULL)
    bk1d = new dcomplex[ncz/2 + 1];

  for(ix=0; ix < ngx; ix++) {
    ZFFT(b[ix], zShift[ix][data.jy], bk1d);
    for(kz = 0; kz <= laplace_maxmode; kz++)
      data.bk[kz][ix] = bk1d[kz];
  }

  /// Create the matrices to be inverted (one for each z point)

  /// Set matrix elements
  par_tridag_matrix(data.avec, data.bvec, data.cvec,
		    data.bk, data.jy, flags, a, ccoef);

  for(kz = 0; kz <= laplace_maxmode; kz++) {
    // Start PDD algorithm

    // Solve for xtilde, v and w (step 2)

    static dcomplex *e = NULL;
    if(e == NULL) {
      e = new dcomplex[ngx];
      for(ix=0;ix<ngx;ix++)
	e[ix] = 0.0;
    }

    dcomplex v0, x0; // Values to be sent to processor i-1

    if(PE_XIND == 0) {
      // Domain includes inner boundary
      tridag(data.avec[kz], data.bvec[kz], data.cvec[kz], 
	     data.bk[kz], data.xk[kz], MXG+MXSUB);
      
      // Add C (row m-1) from next processor
      
      e[MXG+MXSUB-1] = data.cvec[kz][MXG+MXSUB-1];
      tridag(data.avec[kz], data.bvec[kz], data.cvec[kz], 
	     e, data.w[kz], MXG+MXSUB);

    }else if(PE_XIND == (NXPE - 1)) {
      // Domain includes outer boundary
      tridag(data.avec[kz]+MXG, data.bvec[kz]+MXG, data.cvec[kz]+MXG, 
	     data.bk[kz]+MXG, data.xk[kz]+MXG, MXSUB+MXG);
      
      // Add A (row 0) from previous processor
      e[0] = data.avec[kz][MXG];
      tridag(data.avec[kz]+MXG, data.bvec[kz]+MXG, data.cvec[kz]+MXG, 
	     e, data.v[kz]+MXG, MXSUB+MXG);
      
      x0 = data.xk[kz][MXG];
      v0 = data.v[kz][MXG];

    }else {
      // No boundaries
      tridag(data.avec[kz]+MXG, data.bvec[kz]+MXG, data.cvec[kz]+MXG, 
	     data.bk[kz]+MXG, data.xk[kz]+MXG, MXSUB);

      // Add A (row 0) from previous processor
      e[0] = data.avec[kz][MXG];
      tridag(data.avec[kz]+MXG, data.bvec[kz]+MXG, data.cvec[kz]+MXG, 
	     e+MXG, data.v[kz]+MXG, MXSUB);
      e[0] = 0.0;
      
      // Add C (row m-1) from next processor
      e[MXG+MXSUB-1] = data.cvec[kz][MXG+MXSUB-1];
      tridag(data.avec[kz]+MXG, data.bvec[kz]+MXG, data.cvec[kz]+MXG, 
	     e+MXG, data.v[kz]+MXG, MXSUB);
      e[MXG+MXSUB-1] = 0.0;
    }
    
    // Put values into communication buffers
    data.snd[4*kz]   = x0.Real();
    data.snd[4*kz+1] = x0.Imag();
    data.snd[4*kz+2] = v0.Real();
    data.snd[4*kz+3] = v0.Imag();
  }
  
  // Stage 3: Communicate x0, v0 from node i to i-1
  
  if(PE_XIND != (NXPE-1)) {
    // All except the last processor expect to receive data
    // Post async receive
    MPI_Irecv(data.rcv,
	      4*(laplace_maxmode+1),
	      PVEC_REAL_MPI_TYPE,
	      PROC_NUM(PE_XIND+1, PE_YIND), // from processor + 1
	      PDD_COMM_XV,
	      MPI_COMM_WORLD,
	      &data.rcv_req);
  }

  if(PE_XIND != 0) {
    // Send the data
    
    if(invert_async_send) {
      MPI_Isend(data.snd, 
		4*(laplace_maxmode+1),
		PVEC_REAL_MPI_TYPE,
		PROC_NUM(PE_XIND-1, PE_YIND),
		PDD_COMM_XV,
		MPI_COMM_WORLD,
		&data.snd_req);
    }else
      MPI_Send(data.snd, 
	       4*(laplace_maxmode+1),
	       PVEC_REAL_MPI_TYPE,
	       PROC_NUM(PE_XIND-1, PE_YIND),
	       PDD_COMM_XV,
	       MPI_COMM_WORLD);
  }

  return 0;
}

/// Middle part of the PDD algorithm
int invert_pdd_continue(PDD_data &data)
{
  // Wait for x0 and v0 to arrive from processor i+1
  
  MPI_Status status;
  
  if(PE_XIND != (NXPE-1)) {
    MPI_Wait(&data.rcv_req, &status);

    /*! Now solving on all except the last processor
     * 
     * |    1       w^(i)_(m-1) | | y_{2i}   | = | x^(i)_{m-1} |
     * | v^(i+1)_0       1      | | y_{2i+1} |   | x^(i+1)_0   |
     *
     * Only interested in the value of y_2i however
     */
    
    for(int kz = 0; kz <= laplace_maxmode; kz++) {
      dcomplex v0, x0;
      
      // Get x and v0 from processor
      x0 = dcomplex(data.rcv[4*kz], data.rcv[4*kz+1]);
      v0 = dcomplex(data.rcv[4*kz+2], data.rcv[4*kz+3]);
      
      data.y2i[kz] = (data.xk[kz][MXG+MXSUB-1] - data.w[kz][MXG+MXSUB-1]*x0) / (1. - data.w[kz][MXG+MXSUB-1]*v0);
      
    }
  }
  
  if(PE_XIND != 0) {
    // All except pe=0 receive values from i-1. Posting async receive
    MPI_Irecv(data.rcv,
	      2*(laplace_maxmode+1),
	      PVEC_REAL_MPI_TYPE,
	      PROC_NUM(PE_XIND-1, PE_YIND), // from processor - 1
	      PDD_COMM_Y,
	      MPI_COMM_WORLD,
	      &data.rcv_req);
  }


  if(PE_XIND != (NXPE-1)) {
    // Send value to the (i+1)th processor

    if(invert_async_send && (PE_XIND != 0)) // Wait for the previous send to finish before changing the send buffer
      MPI_Wait(&data.snd_req, &status);
    
    for(int kz = 0; kz <= laplace_maxmode; kz++) {
      data.snd[2*kz]   = data.y2i[kz].Real();
      data.snd[2*kz+1] = data.y2i[kz].Imag();
    }
    
    if(invert_async_send) {
      MPI_Isend(data.snd, 
		2*(laplace_maxmode+1),
		PVEC_REAL_MPI_TYPE,
		PROC_NUM(PE_XIND+1, PE_YIND),
		PDD_COMM_Y,
		MPI_COMM_WORLD,
		&data.snd_req);
    }else
      MPI_Send(data.snd, 
	       2*(laplace_maxmode+1),
	       PVEC_REAL_MPI_TYPE,
	       PROC_NUM(PE_XIND+1, PE_YIND),
	       PDD_COMM_Y,
	       MPI_COMM_WORLD);
  }
  
  return 0;
}

/// Last part of the PDD algorithm
int invert_pdd_finish(PDD_data &data, int flags, FieldPerp &x)
{
  int ix, kz;
  MPI_Status status;

  x.Allocate();
  x.setIndex(data.jy);
  
  if(PE_XIND != (NXPE-1)) {
    for(kz = 0; kz <= laplace_maxmode; kz++) {
      for(ix=0; ix < ngx; ix++)
	data.xk[kz][ix] -= data.w[kz][ix] * data.y2i[kz];
    }
  }

  if(PE_XIND != 0) {
    MPI_Wait(&data.rcv_req, &status);
  
    for(kz = 0; kz <= laplace_maxmode; kz++) {
      dcomplex y2m = dcomplex(data.rcv[2*kz], data.rcv[2*kz+1]);
      
      for(ix=0; ix < ngx; ix++)
	data.xk[kz][ix] -= data.v[kz][ix] * y2m;
    }
  }
  
  // Have result in Fourier space. Convert back to real space

  static dcomplex *xk1d = NULL; ///< 1D in Z for taking FFTs

  if(xk1d == NULL) {
    xk1d = new dcomplex[ncz/2 + 1];
    for(kz=0;kz<=ncz/2;kz++)
      xk1d[kz] = 0.0;
  }

  for(ix=0; ix<=ncx; ix++){
    
    for(kz = 0; kz<= laplace_maxmode; kz++) {
      xk1d[kz] = data.xk[kz][ix];
    }

    if(flags & INVERT_ZERO_DC)
      xk1d[0] = 0.0;

    ZFFT_rev(xk1d, zShift[ix][data.jy], x[ix]);
    
    x[ix][ncz] = x[ix][0]; // enforce periodicity
  }

  // Make sure all communication has completed
  if(invert_async_send)
    MPI_Wait(&data.snd_req, &status);

  return 0;
}

/**********************************************************************************
 *                              EXTERNAL INTERFACE
 **********************************************************************************/

/// Invert FieldPerp 
int invert_laplace(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *c)
{
  if(NXPE == 1) {
    // Just use the serial code
    return invert_laplace_ser(b, x, flags, a, c);
  }else {
    // Parallel inversion using PDD

    if(invert_use_pdd) {
      static PDD_data data;
      static bool allocated = false;
      if(!allocated) {
	data.bk = NULL;
	allocated = true;
      }
      
      invert_pdd_start(b, flags, a, data);
      invert_pdd_continue(data);
      invert_pdd_finish(data, flags, x);
    }else {
      static SPT_data data;
      static bool allocated = false;
      if(!allocated) {
	data.bk = NULL;
	allocated = true;
      }
      
      invert_spt_start(b, flags, a, data, c);
      invert_spt_finish(data, flags, x);
    }
  }

  return 0;
}

/// Extracts perpendicular slices from 3D fields and inverts separately
/*!
 * In parallel (NXPE > 1) this tries to overlap computation and communication.
 * This is done at the expense of more memory useage. Setting low_mem
 * in the config file uses less memory, and less communication overlap
 */
int invert_laplace(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c)
{
  int jy, jy2;
  FieldPerp xperp;
  int ret;
  real t;
  
  t = MPI_Wtime();
  
  x.Allocate();

  int ys = jstart, ye = jend;
 
  if(MYPE_IN_CORE == 0) {
    // NOTE: REFINE THIS TO ONLY SOLVE IN BOUNDARY Y CELLS
    ys = 0;
    ye = ngy-1;
  }
  
  if((NXPE == 1) || invert_low_mem) {
    
    for(jy=ys; jy <= ye; jy++) {
      if((flags & INVERT_IN_SET) || (flags & INVERT_OUT_SET))
	xperp = x.Slice(jy); // Using boundary values
      
      if((ret = invert_laplace(b.Slice(jy), xperp, flags, a, c)))
	return(ret);
      x = xperp;
    }
  }else {
    // Use more memory to overlap calculation and communication
    
    if(invert_use_pdd) {
      
      static PDD_data *data = NULL;
    
      if(data == NULL) {
	data = new PDD_data[ye - ys + 1];
	data -= ys; // Re-number indices to start at jstart
	for(jy=ys;jy<=ye;jy++)
	  data[jy].bk = NULL; // Mark as unallocated for PDD routine
      }

      /// PDD algorithm communicates twice, so done in 3 stages
      
      for(jy=ys; jy <= ye; jy++)
	invert_pdd_start(b.Slice(jy), flags, a, data[jy]);
      
      for(jy=ys; jy <= ye; jy++)
	invert_pdd_continue(data[jy]);
      
      for(jy=ys; jy <= ye; jy++) {
	invert_pdd_finish(data[jy], flags, xperp);
	x = xperp;
      }
      
    }else {
      static SPT_data *data = NULL;
      if(data == NULL) {
	data = new SPT_data[ye - ys + 1];
	data -= ys; // Re-number indices to start at ys
	for(jy=ys;jy<=ye;jy++)
	  data[jy].bk = NULL; // Mark as unallocated for PDD routine
      }
      
      
      for(jy=ys; jy <= ye; jy++) {	
	// And start another one going
	invert_spt_start(b.Slice(jy), flags, a, data[jy]);

	// Move each calculation along one processor
	for(jy2=ys; jy2 < jy; jy2++) 
	  invert_spt_continue(data[jy2]);
      }
      
      bool running = true;
      do {
	// Move each calculation along until the last one is finished
	for(jy=ys; jy <= ye; jy++)
	  running = invert_spt_continue(data[jy]) == 0;
      }while(running);
      
      // All calculations finished. Get result
      for(jy=ys; jy <= ye; jy++) {
	invert_spt_finish(data[jy], flags, xperp);
	x = xperp;
      }
      
    }
  }

  wtime_invert += MPI_Wtime() - t;

  x.setLocation(b.getLocation());

  return 0;
}
const Field3D invert_laplace(const Field3D &b, int flags, const Field2D *a, const Field2D *c)
{
  Field3D x;
  
  invert_laplace(b, x, flags, a, c);
  return x;
}

