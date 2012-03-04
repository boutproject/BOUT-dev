/*!
 * \file invert_laplace.cpp
 *
 * \brief Perpendicular Laplacian inversion using FFT and Tridiagonal solver
 *
 * Equation solved is \f$d*\nabla^2_\perp x + (1./c)\nabla_perp c\cdot\nabla_\perp x + a x = b \f$, where
 * \f$x\f$ and \f$x\f$ are perpendicular (X-Z) or 3D fields, 
 * and \f$a\f$ and d are 2D fields. If d is not supplied then it is 1
 * 
 * Flags control the boundary conditions (see header file)
 *
 * Parallel inversion done using two methods
 * - Either a simple parallelisation of the serial algorithm (same operations). Reasonably
 *   parallel as long as MYSUB > mesh->NXPE
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

#include <invert_laplace.hxx>
#include <bout_types.hxx>
#include <globals.hxx>
#include <options.hxx>
#include <fft.hxx>
#include <utils.hxx>
#include <dcomplex.hxx>
#include <cmath>
#include <output.hxx>

#include <lapack_routines.hxx> // Tridiagonal & band inversion routines
#include <boutexception.hxx>
#include <bout/sys/timer.hxx>

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
  BoutReal filter; ///< Fraction of Z modes to filter out. Between 0 and 1

  output.write("Initialising Laplacian inversion routines\n");

  Options *options = Options::getRoot();
  
  // Communication options
  Options *commOpts = options->getSection("comms");
  commOpts->get("async", invert_async_send, true);
  
  // Inversion options
  Options *lapOpts = options->getSection("laplace");
  OPTION(lapOpts, filter, 0.2);
  lapOpts->get("low_mem", invert_low_mem, false);
  lapOpts->get("use_pdd", invert_use_pdd, false);
  lapOpts->get("all_terms", laplace_all_terms, false);
  OPTION(lapOpts, laplace_nonuniform, false);

  if(mesh->firstX() && mesh->lastX()) {
    // This processor is both the first and the last in X
    // -> No parallelisation needed
    output.write("\tUsing serial algorithm\n");
    
  }else {
    // Need to use a parallel algorithm
    if(invert_use_pdd) {
      output.write("\tUsing PDD algorithm\n");
    }else
      output.write("\tUsing parallel Thomas algorithm\n");
  }
  
  int ncz = mesh->ngz-1;

  // convert filtering into an integer number of modes
  laplace_maxmode = ROUND((1.0 - filter) * ((double) (ncz / 2)));

  // Can be overriden by max_mode option
  lapOpts->get("max_mode", laplace_maxmode, laplace_maxmode);
  
  if(laplace_maxmode < 0) laplace_maxmode = 0;
  if(laplace_maxmode > ncz/2) laplace_maxmode = ncz/2;
  
  // Broadcast this since rounding errors could cause mismatch across processors
  // THIS LINE CAUSES SEGFAULT ON LLNL GRENDEL
  //MPI_Bcast(&laplace_maxmode, 1, MPI_INT, 0, BoutComm::get());

  return 0;
}

/// Returns the coefficients for a tridiagonal matrix for laplace. Used by Delp2 too
void laplace_tridag_coefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b, dcomplex &c, 
                          const Field2D *ccoef, const Field2D *d)
{
  BoutReal coef1, coef2, coef3, coef4, coef5, kwave;
  
  kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
  
  coef1=mesh->g11[jx][jy];     ///< X 2nd derivative coefficient
  coef2=mesh->g33[jx][jy];     ///< Z 2nd derivative coefficient
  coef3=2.*mesh->g13[jx][jy];  ///< X-Z mixed derivative coefficient

  coef4 = 0.0;
  coef5 = 0.0;
  if(laplace_all_terms) {
    coef4 = mesh->G1[jx][jy]; // X 1st derivative
    coef5 = mesh->G3[jx][jy]; // Z 1st derivative
  }

  if(d != (Field2D*) NULL) {
    // Multiply Delp2 component by a factor
    coef1 *= (*d)[jx][jy];
    coef2 *= (*d)[jx][jy];
    coef3 *= (*d)[jx][jy];
    coef4 *= (*d)[jx][jy];
    coef5 *= (*d)[jx][jy];
  }

  if(laplace_nonuniform) {
    // non-uniform mesh correction
    if((jx != 0) && (jx != (mesh->ngx-1))) {
      //coef4 += mesh->g11[jx][jy]*0.25*( (1.0/dx[jx+1][jy]) - (1.0/dx[jx-1][jy]) )/dx[jx][jy]; // SHOULD BE THIS (?)
      coef4 -= 0.5*((mesh->dx[jx+1][jy] - mesh->dx[jx-1][jy])/SQ(mesh->dx[jx][jy]))*coef1; // BOUT-06 term
    }
  }

  if(ccoef != NULL) {
    // A first order derivative term
    
    if((jx > 0) && (jx < (mesh->ngx-1)))
      coef4 += mesh->g11[jx][jy] * ((*ccoef)[jx+1][jy] - (*ccoef)[jx-1][jy]) / (2.*mesh->dx[jx][jy]*((*ccoef)[jx][jy]));
  }
  
  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    // d2dz2 term
    coef2 += mesh->g11[jx][jy] * mesh->IntShiftTorsion[jx][jy] * mesh->IntShiftTorsion[jx][jy];
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }
  
  coef1 /= SQ(mesh->dx[jx][jy]);
  coef3 /= 2.*mesh->dx[jx][jy];
  coef4 /= 2.*mesh->dx[jx][jy];

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
 * This code is only for serial i.e. mesh->NXPE == 1
 */
int invert_laplace_ser(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a,
                       const Field2D *ccoef=NULL, const Field2D *d=NULL) {
  int ncx = mesh->ngx-1;
  int ncz = mesh->ngz-1;
  
  static dcomplex **bk = NULL, *bk1d;
  static dcomplex **xk, *xk1d;
  int xbndry; // Width of the x boundary

  if(!mesh->firstX() || !mesh->lastX()) {
    output.write("Error: invert_laplace only works for mesh->NXPE = 1\n");
    return 1;
  }
  
  x.allocate();

  int jy = b.getIndex();
  x.setIndex(jy);

  if(bk == NULL) {
    // Allocate memory

    bk = cmatrix(mesh->ngx, ncz/2 + 1);
    bk1d = new dcomplex[mesh->ngx];
    
    xk = cmatrix(mesh->ngx, ncz/2 + 1);
    xk1d = new dcomplex[mesh->ngx];
  }

  xbndry = 2;
  if(flags & INVERT_BNDRY_ONE)
    xbndry = 1;
  
  #pragma omp parallel for
  for(int ix=0;ix<mesh->ngx;ix++) {
    // for fixed ix,jy set a complex vector rho(z)
    
    ZFFT(b[ix], mesh->zShift[ix][jy], bk[ix]);
  }
  
  if(!mesh->periodicX) {
    if(flags & INVERT_IN_SET) {
      // Setting the inner boundary from x
      #pragma omp parallel for
      for(int ix=0;ix<xbndry;ix++)
	ZFFT(x[ix], mesh->zShift[ix][jy], xk[ix]);
    }
    
    if(flags & INVERT_OUT_SET) {
      // Setting the outer boundary from x
      #pragma omp parallel for
      for(int ix=0;ix<xbndry;ix++)
	ZFFT(x[ncx-ix], mesh->zShift[ncx-ix][jy], xk[ncx-ix]);
    }
  }
  
  if((flags & INVERT_4TH_ORDER) && (!mesh->periodicX)) { // Not implemented for parallel calculations or periodic X
    // Use band solver - 4th order

    static dcomplex **A = (dcomplex**) NULL;
    int xstart, xend;

    if(A == (dcomplex**) NULL)
      A = cmatrix(mesh->ngx, 5);
    
    // Get range for 4th order: Need at least 2 each side
    if(xbndry > 1) {
      xstart = xbndry;
      xend = ncx-xbndry;
    }else {
      xstart = 2;
      xend = mesh->ngx-2;
    }
    
    for(int iz=0;iz<=ncz/2;iz++) {
      // solve differential equation in x
    
      BoutReal coef1=0.0, coef2=0.0, coef3=0.0, coef4=0.0, 
        coef5=0.0, coef6=0.0, kwave, flt;
      ///////// PERFORM INVERSION /////////
      
      // shift freqs according to FFT convention
      kwave=iz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
      
      if (iz>laplace_maxmode) flt=0.0; else flt=1.0;

      // set bk1d
      for(int ix=0;ix<mesh->ngx;ix++)
	bk1d[ix] = bk[ix][iz]*flt;

      // Fill in interior points

      for(int ix=xstart;ix<=xend;ix++) {

	// Set coefficients
	coef1 = mesh->g11[ix][jy];  // X 2nd derivative
	coef2 = mesh->g33[ix][jy];  // Z 2nd derivative
	coef3 = mesh->g13[ix][jy];  // X-Z mixed derivatives
	coef4 = 0.0;          // X 1st derivative
	coef5 = 0.0;          // Z 1st derivative
	coef6 = 0.0;          // Constant

        if(d != (Field2D*) NULL) {
          // Multiply Delp2 component by a factor
          coef1 *= (*d)[ix][jy];
          coef2 *= (*d)[ix][jy];
          coef3 *= (*d)[ix][jy];
        }

	if(a != (Field2D*) NULL)
	  coef6 = (*a)[ix][jy];
	
	if(laplace_all_terms) {
	  coef4 = mesh->G1[ix][jy];
	  coef5 = mesh->G3[ix][jy];
	}

	if(laplace_nonuniform) {
	  // non-uniform mesh correction
	  if((ix != 0) && (ix != ncx))
	    coef4 += mesh->g11[ix][jy]*( (1.0/mesh->dx[ix+1][jy]) - (1.0/mesh->dx[ix-1][jy]) )/(2.0*mesh->dx[ix][jy]);
	}

	if(ccoef != NULL) {
	  // A first order derivative term (1/c)\nabla_perp c\cdot\nabla_\perp x
    
	  if((ix > 1) && (ix < (mesh->ngx-2)))
	    coef4 += mesh->g11[ix][jy] * ((*ccoef)[ix-2][jy] - 8.*(*ccoef)[ix-1][jy] + 8.*(*ccoef)[ix+1][jy] - (*ccoef)[ix+2][jy]) / (12.*mesh->dx[ix][jy]*((*ccoef)[ix][jy]));
	}

	// Put into matrix
	coef1 /= 12.* SQ(mesh->dx[ix][jy]);
	coef2 *= SQ(kwave);
	coef3 *= kwave / (12. * mesh->dx[ix][jy]);
	coef4 /= 12. * mesh->dx[ix][jy];
	coef5 *= kwave;

	A[ix][0] = dcomplex(    -coef1 +   coef4 ,     coef3 );
	A[ix][1] = dcomplex( 16.*coef1 - 8*coef4 , -8.*coef3 );
	A[ix][2] = dcomplex(-30.*coef1 - coef2 + coef6, coef5);
	A[ix][3] = dcomplex( 16.*coef1 + 8*coef4 ,  8.*coef3 );
	A[ix][4] = dcomplex(    -coef1 -   coef4 ,    -coef3 );
      }

      if(xbndry < 2) {
	// Use 2nd order near edges

	int ix = 1;

	coef1=mesh->g11[ix][jy]/(SQ(mesh->dx[ix][jy]));
	coef2=mesh->g33[ix][jy];
	coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);
        
        if(d != (Field2D*) NULL) {
          // Multiply Delp2 component by a factor
          coef1 *= (*d)[ix][jy];
          coef2 *= (*d)[ix][jy];
          coef3 *= (*d)[ix][jy];
        }
        
	A[ix][0] = 0.0; // Should never be used
	A[ix][1] = dcomplex(coef1, -coef3);
	A[ix][2] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	A[ix][3] = dcomplex(coef1,  coef3);
	A[ix][4] = 0.0;

	ix = ncx-1;

	coef1=mesh->g11[ix][jy]/(SQ(mesh->dx[ix][jy]));
	coef2=mesh->g33[ix][jy];
	coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);

	A[ix][0] = 0.0;
	A[ix][1] = dcomplex(coef1, -coef3);
	A[ix][2] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	A[ix][3] = dcomplex(coef1,  coef3);
	A[ix][4] = 0.0;  // Should never be used
      }

      // Boundary conditions

      for(int ix=0;ix<xbndry;ix++) {
	// Set zero-value. Change to zero-gradient if needed

        if(!(flags & INVERT_IN_RHS))
          bk1d[ix] = 0.0;
        if(!(flags & INVERT_OUT_RHS))
          bk1d[ncx-ix] = 0.0;

	A[ix][0] = A[ix][1] = A[ix][3] = A[ix][4] = 0.0;
	A[ix][2] = 1.0;

	A[ncx-ix][0] = A[ncx-ix][1] = A[ncx-ix][3] = A[ncx-ix][4] = 0.0;
	A[ncx-ix][2] = 1.0;
      }

      if(flags & INVERT_IN_SET) {
	// Set values of inner boundary from X
	for(int ix=0;ix<xbndry;ix++)
	  bk1d[ix] = xk[ix][iz];
      }
      
      if(flags & INVERT_OUT_SET) {
	// Set values of outer boundary from X
	for(int ix=0;ix<xbndry;ix++)
	  bk1d[ncx-ix] = xk[ncx-ix][iz];
      }

      if(iz == 0) {
	// DC
	
	// Inner boundary
	if(flags & INVERT_DC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  for (int ix=0;ix<xbndry;ix++)
	    A[ix][3] = -1.0;
	}
	
	// Outer boundary
	if(flags & INVERT_DC_OUT_GRAD) {
	  // Zero gradient at outer boundary
	  for (int ix=0;ix<xbndry;ix++)
	    A[ncx-ix][1] = -1.0;
	}
	
      }else {
	// AC
	
	// Inner boundarySQ(kwave)*coef2
	if(flags & INVERT_AC_IN_GRAD) {
	  // Zero gradient at inner boundary
	  for (int ix=0;ix<xbndry;ix++)
	    A[ix][3] = -1.0;
	}else if(flags & INVERT_AC_IN_LAP) {
	  // Enforce zero laplacian for 2nd and 4th-order
	  
	  int ix = 1;
	  
	  coef1=mesh->g11[ix][jy]/(12.* SQ(mesh->dx[ix][jy]));
	
	  coef2=mesh->g33[ix][jy];
	
	  coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);
	  
	  coef4 = 0.0;
	  if(a != (Field2D*) NULL)
	    coef4 = (*a)[ix][jy];
	  
	  // Combine 4th order at 1 with 2nd order at 0
	  A[1][0] = 0.0; // Not used
	  A[1][1] = dcomplex( (14. - SQ(mesh->dx[0][jy]*kwave)*mesh->g33[0][jy]/mesh->g11[0][jy])*coef1  ,  -coef3 );
	  A[1][2] = dcomplex(-29.*coef1 - SQ(kwave)*coef2 + coef4, 0.0);
	  A[1][3] = dcomplex( 16.*coef1  , coef3 );
	  A[1][4] = dcomplex(    -coef1  ,     0.0 );
	  
	  coef1=mesh->g11[ix][jy]/(SQ(mesh->dx[ix][jy]));
	  coef2=mesh->g33[ix][jy];
	  coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);

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
	  for (int ix=0;ix<xbndry;ix++)
	    A[ncx-ix][1] = -1.0;
	}else if(flags & INVERT_AC_OUT_LAP) {
	  // Enforce zero laplacian for 2nd and 4th-order
	  // NOTE: Currently ignoring XZ term and coef4 assumed zero on boundary
	  // FIX THIS IF IT WORKS

	  int ix = ncx-1;
	  
	  coef1=mesh->g11[ix][jy]/(12.* SQ(mesh->dx[ix][jy]));
	
	  coef2=mesh->g33[ix][jy];
	
	  coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);
	  
	  coef4 = 0.0;
	  if(a != (Field2D*) NULL)
	    coef4 = (*a)[ix][jy];
	  
	  // Combine 4th order at ncx-1 with 2nd order at ncx
	  A[ix][0] = dcomplex(    -coef1  ,     0.0 );
	  A[ix][1] = dcomplex( 16.*coef1  , -coef3 );
	  A[ix][2] = dcomplex(-29.*coef1 - SQ(kwave)*coef2 + coef4, 0.0);
	  A[ix][3] = dcomplex( (14. - SQ(mesh->dx[ncx][jy]*kwave)*mesh->g33[ncx][jy]/mesh->g11[ncx][jy])*coef1  ,  coef3 );
	  A[ix][4] = 0.0; // Not used
	  
	  coef1=mesh->g11[ix][jy]/(SQ(mesh->dx[ix][jy]));
	  coef2=mesh->g33[ix][jy];
	  coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);

	  // Use 2nd order at ncx - 1
	  A[ncx][0] = dcomplex(coef1, -coef3);
	  A[ncx][1] = dcomplex(-2.0*coef1 - SQ(kwave)*coef2 + coef4,0.0);
	  A[ncx][2] = dcomplex(coef1,  coef3);
	  A[ncx][3] = 0.0;  // Should never be used
	  A[ncx][4] = 0.0;
	}
      }
      
      // Perform inversion
      cband_solve(A, mesh->ngx, 2, 2, bk1d);
      
      if((flags & INVERT_KX_ZERO) && (iz == 0)) {
        // Set the Kx = 0, n = 0 component to zero. For now just subtract
        // Should do in the inversion e.g. Sherman-Morrison formula
        
        dcomplex offset(0.0);
        for(int ix=0;ix<=ncx;ix++)
          offset += bk1d[ix];
        offset /= (BoutReal) (ncx+1);
        for(int ix=0;ix<=ncx;ix++)
          bk1d[ix] -= offset;
      }
      
      // Fill xk
      for (int ix=0; ix<=ncx; ix++)
	xk[ix][iz]=bk1d[ix];
      
    }
  }else {
    // Use tridiagonal system in x - 2nd order
    
    static dcomplex *avec = (dcomplex*) NULL, *bvec, *cvec;
    
    if(avec == (dcomplex*) NULL) {
      avec = new dcomplex[mesh->ngx];
      bvec = new dcomplex[mesh->ngx];
      cvec = new dcomplex[mesh->ngx];
    }

    for(int iz=0;iz<=ncz/2;iz++) {
      // solve differential equation in x

      // set bk1d
      BoutReal flt;
      if (iz>laplace_maxmode) flt=0.0; else flt=1.0;
      
      for(int ix=0;ix<=ncx;ix++)
	bk1d[ix] = bk[ix][iz] * flt;

      ///////// PERFORM INVERSION /////////
      
      for(int ix=xbndry;ix<=ncx-xbndry;ix++) {
	laplace_tridag_coefs(ix, jy, iz, avec[ix], bvec[ix], cvec[ix], ccoef, d);
	
	if(a != (Field2D*) NULL)
	  bvec[ix] += (*a)[ix][jy];
      }

      if(!mesh->periodicX) {
	// Need boundary conditions
	/// By default, set RHS to zero, unless INVERT_*_RHS set
	if(!(flags & INVERT_IN_RHS)) {
	  for(int ix=0;ix<xbndry;ix++)
	    bk1d[ix] = 0.;
	}
	if(!(flags & INVERT_OUT_RHS)) {
	  for(int ix=mesh->ngx-xbndry;ix<mesh->ngx;ix++)
	    bk1d[ix] = 0.;
	}
	
	// Set boundary conditions
	
	if(iz == 0) {
	  // DC
	  
	  // Inner boundary
	  if(flags & INVERT_DC_IN_GRAD) {
	    // Zero gradient at inner boundary
	    
	    if((flags & INVERT_IN_SYM) && (xbndry > 1) && mesh->BoundaryOnCell) {
	      // Use symmetric boundary to set zero-gradient
	      
	      for (int ix=0;ix<xbndry-1;ix++) {
		avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      }
	      // Symmetric on last point
	      avec[xbndry-1] = 1.0; bvec[xbndry-1] = 0.0; cvec[xbndry-1] = -1.0;
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		avec[ix]=dcomplex(0.0,0.0);
		bvec[ix]=dcomplex(1.,0.);
		cvec[ix]=dcomplex(-1.,0.);
	      }
	    }
	  }else if(flags & INVERT_IN_SET) {
	    for(int ix=0;ix<xbndry;ix++) {
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
	      for(int ix=0;ix<xbndry-1;ix++) {
		avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      }

	      if(mesh->BoundaryOnCell) {
		// Antisymmetric about boundary on cell
		avec[xbndry-1]=1.0; bvec[xbndry-1]=0.0; cvec[xbndry-1]= 1.0;
	      }else { 
		// Antisymmetric across boundary between cells
		avec[xbndry-1]=0.0; bvec[xbndry-1]=1.0; cvec[xbndry-1]= 1.0;
	      }
	    
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		avec[ix]=dcomplex(0.,0.);
		bvec[ix]=dcomplex(1.,0.);
		cvec[ix]=dcomplex(0.,0.);
	      }
	    }
	  }
	
	  // Outer boundary
	  if(flags & INVERT_DC_OUT_GRAD) {
	    // Zero gradient at outer boundary

	    if((flags & INVERT_OUT_SYM) && (xbndry > 1) && mesh->BoundaryOnCell) {
	      // Use symmetric boundary to set zero-gradient
	    
	      for (int ix=0;ix<xbndry-1;ix++) {
		avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	      // Symmetric on last point
	      int ix = xbndry-1;
	      avec[ncx-ix] = 1.0; bvec[ncx-ix] = 0.0; cvec[ncx-ix] = -1.0;
	    
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		cvec[ncx-ix]=dcomplex(0.,0.);
		bvec[ncx-ix]=dcomplex(1.,0.);
		avec[ncx-ix]=dcomplex(-1.,0.);
	      }
	    }
	  }else if(flags & INVERT_OUT_SET) {
	    // Setting the values in the outer boundary
	    for(int ix=0;ix<xbndry;ix++) {
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
	      for(int ix=0;ix<xbndry-1;ix++) {
		avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	      int ix = xbndry-1;
	      if(mesh->BoundaryOnCell) {
		// Antisymmetric about boundary on cell
		avec[ncx-ix]=1.0; bvec[ncx-ix]=0.0; cvec[ncx-ix]= 1.0;
	      }else { 
		// Antisymmetric across boundary between cells
		avec[ncx-ix]=1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		cvec[ncx-ix]=dcomplex(0.,0.);
		bvec[ncx-ix]=dcomplex(1.,0.);
		avec[ncx-ix]=dcomplex(0.,0.);
	      }
	    }
	  }
	}else {
	  // AC
	
	  // Inner boundary
	  if(flags & INVERT_AC_IN_GRAD) {
	    // Zero gradient at inner boundary
	  
	    if((flags & INVERT_IN_SYM) && (xbndry > 1) && mesh->BoundaryOnCell) {
	      // Use symmetric boundary to set zero-gradient
	    
	      for (int ix=0;ix<xbndry-1;ix++) {
		avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      }
	      // Symmetric on last point
	      avec[xbndry-1] = 1.0; bvec[xbndry-1] = 0.0; cvec[xbndry-1] = -1.0;
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		avec[ix]=dcomplex(0.,0.);
		bvec[ix]=dcomplex(1.,0.);
		cvec[ix]=dcomplex(-1.,0.);
	      }
	    }
	  }else if(flags & INVERT_IN_SET) {
	    // Setting the values in the boundary
	    for(int ix=0;ix<xbndry;ix++) {
	      avec[ix] = 0.0;
	      bvec[ix] = 1.0;
	      cvec[ix] = 0.0;
	      bk1d[ix] = xk[ix][iz];
	    }
	  }else if(flags & INVERT_AC_IN_LAP) {
	    // Use decaying zero-Laplacian solution in the boundary
	    BoutReal kwave=iz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	    for (int ix=0;ix<xbndry;ix++) {
	      avec[ix] = 0.0;
	      bvec[ix] = -1.0;
	      cvec[ix] = exp(-1.0*sqrt(mesh->g33[ix][jy]/mesh->g11[ix][jy])*kwave*mesh->dx[ix][jy]);
	    }
	  }else {
	    // Zero value at inner boundary

	    if(flags & INVERT_IN_SYM) {
	      // Use anti-symmetric boundary to set zero-value
	    
	      // Zero-gradient for first point(s)
	      for(int ix=0;ix<xbndry-1;ix++) {
		avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
	      }

	      if(mesh->BoundaryOnCell) {
		// Antisymmetric about boundary on cell
		avec[xbndry-1]=1.0; bvec[xbndry-1]=0.0; cvec[xbndry-1]= 1.0;
	      }else { 
		// Antisymmetric across boundary between cells
		avec[xbndry-1]=0.0; bvec[xbndry-1]=1.0; cvec[xbndry-1]= 1.0;
	      }
	    
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		avec[ix]=dcomplex(0.,0.);
		bvec[ix]=dcomplex(1.,0.);
		cvec[ix]=dcomplex(0.,0.);
	      }
	    }
	  }
	
	  // Outer boundary
	  if(flags & INVERT_AC_OUT_GRAD) {
	    // Zero gradient at outer boundary
	  
	    if((flags & INVERT_OUT_SYM) && (xbndry > 1) && mesh->BoundaryOnCell) {
	      // Use symmetric boundary to set zero-gradient
	    
	      for (int ix=0;ix<xbndry-1;ix++) {
		avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	      // Symmetric on last point
	      int ix = xbndry-1;
	      avec[ncx-ix] = 1.0; bvec[ncx-ix] = 0.0; cvec[ncx-ix] = -1.0;
	    
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		cvec[ncx-ix]=dcomplex(0.,0.);
		bvec[ncx-ix]=dcomplex(1.,0.);
		avec[ncx-ix]=dcomplex(-1.,0.);
	      }
	    }
	  }else if(flags & INVERT_AC_OUT_LAP) {
	    // Use decaying zero-Laplacian solution in the boundary
	    BoutReal kwave=iz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	    for (int ix=0;ix<xbndry;ix++) {
	      avec[ncx-ix] = exp(-1.0*sqrt(mesh->g33[ncx-ix][jy]/mesh->g11[ncx-ix][jy])*kwave*mesh->dx[ncx-ix][jy]);;
	      bvec[ncx-ix] = -1.0;
	      cvec[ncx-ix] = 0.0;
	    }
	  }else if(flags & INVERT_OUT_SET) {
	    // Setting the values in the outer boundary
	    for(int ix=0;ix<xbndry;ix++) {
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
	      for(int ix=0;ix<xbndry-1;ix++) {
		avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	      int ix = xbndry-1;
	      if(mesh->BoundaryOnCell) {
		// Antisymmetric about boundary on cell
		avec[ncx-ix]=1.0; bvec[ncx-ix]=0.0; cvec[ncx-ix]= 1.0;
	      }else {
		// Antisymmetric across boundary between cells
		avec[ncx-ix]=1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
	      }
	    }else {
	      for (int ix=0;ix<xbndry;ix++){
		cvec[ncx-ix]=dcomplex(0.,0.);
		bvec[ncx-ix]=dcomplex(1.,0.);
		avec[ncx-ix]=dcomplex(0.,0.);
	      }
	    }
	  }
	}
        
	// Call tridiagonal solver
	tridag(avec, bvec, cvec, bk1d, xk1d, mesh->ngx);

	if((flags & INVERT_IN_SYM) && (xbndry > 1)) {
	  // (Anti-)symmetry on inner boundary. Nothing to do if only one boundary cell
	  int xloc = 2*xbndry;
	  if(!mesh->BoundaryOnCell)
	    xloc--;
	  
	  if( ((iz == 0) && (flags & INVERT_DC_IN_GRAD)) || ((iz != 0) && (flags & INVERT_AC_IN_GRAD)) ) {
	    // Inner gradient zero - symmetric
	    for(int ix=0;ix<xbndry-1;ix++)
	      xk1d[ix] = xk1d[xloc-ix];
	  }else {
	    // Inner value zero - antisymmetric
	    for(int ix=0;ix<xbndry-1;ix++)
	      xk1d[ix] = -xk1d[xloc-ix];
	  }
	}
	if((flags & INVERT_OUT_SYM) && (xbndry > 1)) {
	  // (Anti-)symmetry on outer boundary. Nothing to do if only one boundary cell
	  
	  int xloc =  mesh->ngx - 2*xbndry;
	  if(mesh->BoundaryOnCell)
	    xloc--;
	
	  if( ((iz == 0) && (flags & INVERT_DC_IN_GRAD)) || ((iz != 0) && (flags & INVERT_AC_IN_GRAD)) ) {
	    // Outer gradient zero - symmetric
	    for(int ix=0;ix<xbndry-1;ix++)
	      xk1d[ncx-ix] = xk1d[xloc + ix];
	  }else {
	    // Outer value zero - antisymmetric
	    for(int ix=0;ix<xbndry-1;ix++)
	      xk1d[ncx-ix] = -xk1d[xloc + ix];
	  }
	}
      } else {
	// Periodic in X, so no boundaries
	cyclic_tridag(avec+2, bvec+2, cvec+2, bk1d+2, xk1d+2, mesh->ngx-4);
	
	// Copy boundary regions
	for(int ix=0;ix<2;ix++) {
	  xk1d[ix] = xk1d[mesh->ngx-4+ix];
	  xk1d[mesh->ngx-2+ix] = xk1d[2+ix];
	}
      }
      
      if((flags & INVERT_KX_ZERO) && (iz == 0)) {
        dcomplex offset(0.0);
        for(int ix=0;ix<=ncx;ix++)
          offset += bk1d[ix];
        offset /= (BoutReal) (ncx+1);
        for(int ix=0;ix<=ncx;ix++)
          bk1d[ix] -= offset;
      }
      
      // Fill xk
      
      for (int ix=0; ix<=ncx; ix++){
	xk[ix][iz]=xk1d[ix];
      }
    }
  }

  // Done inversion, transform back

  for(int ix=0; ix<=ncx; ix++){
    
    if(flags & INVERT_ZERO_DC)
      xk[ix][0] = 0.0;

    ZFFT_rev(xk[ix], mesh->zShift[ix][jy], x[ix]);
    
    x[ix][mesh->ngz-1] = x[ix][0]; // enforce periodicity
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
		       dcomplex **bk, int jy, int flags, 
                       const Field2D *a = NULL, const Field2D *ccoef=NULL, const Field2D *d=NULL)
{
  int ix, kz;
  
  int ncx = mesh->ngx-1;

  int xbndry = 2;
  if(flags & INVERT_BNDRY_ONE)
    xbndry = 1;
  
  for(kz = 0; kz <= laplace_maxmode; kz++) {
    
    // Entire domain. Change to boundaries later

    for(ix=0;ix<=ncx;ix++) {

      laplace_tridag_coefs(ix, jy, kz, avec[kz][ix], bvec[kz][ix], cvec[kz][ix], ccoef, d);
	
      if(a != (Field2D*) NULL)
	bvec[kz][ix] += (*a)[ix][jy];
    }

    if(!mesh->periodicX) {
      // Boundary conditions
      
      if(mesh->firstX()) {
	// INNER BOUNDARY ON THIS PROCESSOR
	
	if(!(flags & INVERT_IN_RHS)) {
	  for(ix=0;ix<xbndry;ix++)
	    bk[kz][ix] = 0.;
	}
	
	if(kz == 0) {
	  // DC
	  
	  if(flags & INVERT_DC_IN_GRAD) {
	    // Zero gradient at inner boundary
	    for (ix=0;ix<xbndry;ix++){
	      avec[kz][ix] =  0.;
	      bvec[kz][ix] =  1.;
	      cvec[kz][ix] = -1.;
	    }
	  }else {
	    // Zero value at inner boundary
	    for (ix=0;ix<xbndry;ix++){
	      avec[kz][ix] = 0.;
	      bvec[kz][ix] = 1.;
	      cvec[kz][ix] = 0.;
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
	    }
	  }else if(flags & INVERT_AC_IN_LAP) {
	    // Use decaying zero-Laplacian solution in the boundary
	    BoutReal kwave=kz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	    for (ix=0;ix<xbndry;ix++) {
	      avec[kz][ix] = 0.0;
	      bvec[kz][ix] = 1.0;
	      cvec[kz][ix] = -exp(-1.0*sqrt(mesh->g33[ix][jy]/mesh->g11[ix][jy])*kwave*mesh->dx[ix][jy]);
	    }
	  }else {
	    // Zero value at inner boundary
	    for (ix=0;ix<xbndry;ix++){
	      avec[kz][ix]=dcomplex(0.,0.);
	      bvec[kz][ix]=dcomplex(1.,0.);
	      cvec[kz][ix]=dcomplex(0.,0.);
	    }
	  }
	}
      }else if(mesh->lastX()) {
	// OUTER BOUNDARY
      
	if(!(flags & INVERT_OUT_RHS)) {
	  for (ix=0;ix<xbndry;ix++)
	    bk[kz][ncx-ix] = 0.;
	}

	if(kz == 0) {
	  // DC
	
	  if(flags & INVERT_DC_OUT_GRAD) {
	    // Zero gradient at outer boundary
	    for (ix=0;ix<xbndry;ix++){
	      cvec[kz][ncx-ix]=dcomplex(0.,0.);
	      bvec[kz][ncx-ix]=dcomplex(1.,0.);
	      avec[kz][ncx-ix]=dcomplex(-1.,0.);
	    }
	  }else {
	    // Zero value at outer boundary
	    for (ix=0;ix<xbndry;ix++){
	      cvec[kz][ncx-ix]=dcomplex(0.,0.);
	      bvec[kz][ncx-ix]=dcomplex(1.,0.);
	      avec[kz][ncx-ix]=dcomplex(0.,0.);
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
	    }
	  }else if(flags & INVERT_AC_OUT_LAP) {
	    // Use decaying zero-Laplacian solution in the boundary
	    BoutReal kwave=kz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	    for (ix=0;ix<xbndry;ix++) {
	      avec[kz][ncx-ix] = -exp(-1.0*sqrt(mesh->g33[ncx-ix][jy]/mesh->g11[ncx-ix][jy])*kwave*mesh->dx[ncx-ix][jy]);;
	      bvec[kz][ncx-ix] = 1.0;
	      cvec[kz][ncx-ix] = 0.0;
	    }
	  }else {
	    // Zero value at outer boundary
	    for (ix=0;ix<xbndry;ix++){
	      cvec[kz][ncx-ix]=dcomplex(0.,0.);
	      bvec[kz][ncx-ix]=dcomplex(1.,0.);
	      avec[kz][ncx-ix]=dcomplex(0.,0.);
	    }
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
  
  comm_handle recv_handle; // Handle for receives
  
  int comm_tag; // Tag for communication
  
  BoutReal *buffer;
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
    if(bet == 0.0)
      throw BoutException("Tridag: Zero pivot\n");
    
    u[j] = (r[j]-a[j]*u[j-1])/bet;
  }

  um = u[n-1];
}

/// Second (backsolve) part of the Thomas algorithm
/*!
 * @param[inout] u    Result to be solved (Au = r)
 * @param[in]    n    Size of the problem
 * @param[in]    gam  Intermediate values produced by the forward part
 * @param[inout] gp   gam from the processor mesh->PE_XIND + 1, and returned to mesh->PE_XIND - 1
 * @param[inout] up   u from processor mesh->PE_XIND + 1, and returned to mesh->PE_XIND - 1
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
 * processor (mesh->PE_XIND = 0) to the outermost (mesh->PE_XIND = mesh->NXPE-1) and back again.
 *
 * Some parallelism is achieved by running several inversions simultaneously, so while
 * processor #1 is inverting Y=0, processor #0 is starting on Y=1. This works ok as long
 * as the number of slices to be inverted is greater than the number of X processors (MYSUB > mesh->NXPE).
 * If MYSUB < mesh->NXPE then not all processors can be busy at once, and so efficiency will fall sharply.
 *
 * @param[in]    b      RHS values (Ax = b)
 * @param[in]    flags  Inversion settings (see boundary.h for values)
 * @param[in]    a      This is a 2D matrix which allows solution of A = Delp2 + a
 * @param[out]   data   Structure containing data needed for second half of inversion
 * @param[in]    ccoef  Optional coefficient for first-order derivative
 * @param[in]    d      Optional factor to multiply the Delp2 operator
 */
int invert_spt_start(const FieldPerp &b, int flags, const Field2D *a, SPT_data &data, 
                     const Field2D *ccoef = NULL, const Field2D *d = NULL)
{
  if(mesh->NXPE == 1)
    throw BoutException("Error: SPT method only works for mesh->NXPE > 1\n");

  data.jy = b.getIndex();

  if(data.bk == NULL) {
    /// Allocate memory
    
    // RHS vector
    data.bk = cmatrix(laplace_maxmode + 1, mesh->ngx);
    data.xk = cmatrix(laplace_maxmode + 1, mesh->ngx);
    
    data.gam = cmatrix(laplace_maxmode + 1, mesh->ngx);

    // Matrix to be solved
    data.avec = cmatrix(laplace_maxmode + 1, mesh->ngx);
    data.bvec = cmatrix(laplace_maxmode + 1, mesh->ngx);
    data.cvec = cmatrix(laplace_maxmode + 1, mesh->ngx);
    
    data.buffer  = new BoutReal[4*(laplace_maxmode + 1)];
  }

  /// Take FFTs of data
  static dcomplex *bk1d = NULL; ///< 1D in Z for taking FFTs

  int ncz = mesh->ngz-1;

  if(bk1d == NULL)
    bk1d = new dcomplex[ncz/2 + 1];
  
  for(int ix=0; ix < mesh->ngx; ix++) {
    ZFFT(b[ix], mesh->zShift[ix][data.jy], bk1d);
    for(int kz = 0; kz <= laplace_maxmode; kz++)
      data.bk[kz][ix] = bk1d[kz];
  }
  
  /// Set matrix elements
  par_tridag_matrix(data.avec, data.bvec, data.cvec,
		    data.bk, data.jy, flags, a, ccoef, d);

  data.proc = 0; //< Starts at processor 0
  data.dir = 1;
  
  if(mesh->firstX()) {
    dcomplex bet, u0;
    #pragma omp parallel for
    for(int kz = 0; kz <= laplace_maxmode; kz++) {
      // Start tridiagonal solve
      spt_tridag_forward(data.avec[kz], data.bvec[kz], data.cvec[kz],
			 data.bk[kz], data.xk[kz], mesh->xend+1,
			 data.gam[kz],
			 bet, u0, true);
      // Load intermediate values into buffers
      data.buffer[4*kz]     = bet.Real();
      data.buffer[4*kz + 1] = bet.Imag();
      data.buffer[4*kz + 2] = u0.Real();
      data.buffer[4*kz + 3] = u0.Imag();
    }
    
    // Send data
    mesh->sendXOut(data.buffer, 4*(laplace_maxmode+1), data.comm_tag);
    
  }else if(mesh->PE_XIND == 1) {
    // Post a receive
    data.recv_handle = mesh->irecvXIn(data.buffer, 4*(laplace_maxmode+1), data.comm_tag);
  }
  
  data.proc++; // Now moved onto the next processor
  if(mesh->NXPE == 2)	
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
  if(data.proc < 0) // Already finished
    return 1;
  
  if(mesh->PE_XIND == data.proc) {
    /// This processor's turn to do inversion

    // Wait for data to arrive
    mesh->wait(data.recv_handle);

    if(mesh->lastX()) {
      // Last processor, turn-around
      
      #pragma omp parallel for
      for(int kz = 0; kz <= laplace_maxmode; kz++) {
        dcomplex bet, u0;
        dcomplex gp, up;
	bet = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	u0 = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);
	spt_tridag_forward(data.avec[kz]+mesh->xstart,
			   data.bvec[kz]+mesh->xstart, 
			   data.cvec[kz]+mesh->xstart,
			   data.bk[kz]+mesh->xstart, 
			   data.xk[kz]+mesh->xstart, mesh->xend+1,
			   data.gam[kz]+mesh->xstart,
			   bet, u0);
	
	// Back-substitute
	gp = 0.0;
	up = 0.0;
	spt_tridag_back(data.xk[kz]+mesh->xstart, mesh->ngx-mesh->xstart, 
			data.gam[kz]+mesh->xstart, gp, up);
	data.buffer[4*kz]     = gp.Real();
	data.buffer[4*kz + 1] = gp.Imag();
	data.buffer[4*kz + 2] = up.Real();
	data.buffer[4*kz + 3] = up.Imag();
      }

    }else if(data.dir > 0) {
      // In the middle of X, forward direction

      #pragma omp parallel for
      for(int kz = 0; kz <= laplace_maxmode; kz++) {
	dcomplex bet, u0;
	bet = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	u0 = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);
	spt_tridag_forward(data.avec[kz]+mesh->xstart, 
			   data.bvec[kz]+mesh->xstart, 
			   data.cvec[kz]+mesh->xstart,
			   data.bk[kz]+mesh->xstart, 
			   data.xk[kz]+mesh->xstart, 
			   mesh->xend - mesh->xstart+1,
			   data.gam[kz]+mesh->xstart,
			   bet, u0);
	// Load intermediate values into buffers
	data.buffer[4*kz]     = bet.Real();
	data.buffer[4*kz + 1] = bet.Imag();
	data.buffer[4*kz + 2] = u0.Real();
	data.buffer[4*kz + 3] = u0.Imag();
      }
      
    }else if(mesh->firstX()) {
      // Back to the start
      
      dcomplex gp, up;
      for(int kz = 0; kz <= laplace_maxmode; kz++) {
	gp = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	up = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);

	spt_tridag_back(data.xk[kz], mesh->xend+1, data.gam[kz], gp, up);
      }

    }else {
      // Middle of X, back-substitution stage

      #pragma omp parallel for
      for(int kz = 0; kz <= laplace_maxmode; kz++) {
	dcomplex gp = dcomplex(data.buffer[4*kz], data.buffer[4*kz + 1]);
	dcomplex up = dcomplex(data.buffer[4*kz + 2], data.buffer[4*kz + 3]);

	spt_tridag_back(data.xk[kz]+mesh->xstart, 
			mesh->xend-mesh->xstart+1, 
			data.gam[kz]+mesh->xstart, gp, up);
	
	data.buffer[4*kz]     = gp.Real();
	data.buffer[4*kz + 1] = gp.Imag();
	data.buffer[4*kz + 2] = up.Real();
	data.buffer[4*kz + 3] = up.Imag();
      }
    }

    if(mesh->PE_XIND != 0) { // If not finished yet
      /// Send data
      
      if(data.dir > 0) {
	mesh->sendXOut(data.buffer, 4*(laplace_maxmode+1), data.comm_tag);
      }else
	mesh->sendXIn(data.buffer, 4*(laplace_maxmode+1), data.comm_tag);
    }

  }else if(mesh->PE_XIND == data.proc + data.dir) {
    // This processor is next, post receive
    
    if(data.dir > 0) {
      data.recv_handle = mesh->irecvXIn(data.buffer, 4*(laplace_maxmode+1), data.comm_tag);
    }else
      data.recv_handle = mesh->irecvXOut(data.buffer, 4*(laplace_maxmode+1), data.comm_tag);
  }
  
  data.proc += data.dir;
  
  if(data.proc == mesh->NXPE-1)
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
  int ncx = mesh->ngx-1;
  int ncz = mesh->ngz-1;

  x.allocate();
  x.setIndex(data.jy);
  BoutReal **xdata = x.getData();

  // Make sure calculation has finished
  while(invert_spt_continue(data) == 0) {}

  // Have result in Fourier space. Convert back to real space

  static dcomplex *xk1d = NULL; ///< 1D in Z for taking FFTs

  if(xk1d == NULL) {
    xk1d = new dcomplex[ncz/2 + 1];
    for(int kz=0;kz<=ncz/2;kz++)
      xk1d[kz] = 0.0;
  }
  
  for(int ix=0; ix<=ncx; ix++){
    
    for(int kz = 0; kz<= laplace_maxmode; kz++) {
      xk1d[kz] = data.xk[kz][ix];
    }

    if(flags & INVERT_ZERO_DC)
      xk1d[0] = 0.0;

    ZFFT_rev(xk1d, mesh->zShift[ix][data.jy], xdata[ix]);
    
    xdata[ix][ncz] = xdata[ix][0]; // enforce periodicity
  }

  if(!mesh->firstX()) {
    // Set left boundary to zero (Prevent unassigned values in corners)
    for(int ix=0; ix<mesh->xstart; ix++){
      for(int kz=0;kz<mesh->ngz;kz++)
	xdata[ix][kz] = 0.0;
    }
  }
  if(!mesh->lastX()) {
    // Same for right boundary
    for(int ix=mesh->xend+1; ix<mesh->ngx; ix++){
      for(int kz=0;kz<mesh->ngz;kz++)
	xdata[ix][kz] = 0.0;
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

  BoutReal *snd; // send buffer
  BoutReal *rcv; // receive buffer
  
  comm_handle recv_handle;

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
int invert_pdd_start(const FieldPerp &b, int flags, const Field2D *a, PDD_data &data, 
                     const Field2D *ccoef = NULL, const Field2D *d=NULL)
{
  int ix, kz;
  
  int ncz = mesh->ngz-1;

  data.jy = b.getIndex();

  if(mesh->firstX() && mesh->lastX()) {
    output.write("Error: PDD method only works for NXPE > 1\n");
    return 1;
  }

  if(data.bk == NULL) {
    // Need to allocate working memory
    
    // RHS vector
    data.bk = cmatrix(laplace_maxmode + 1, mesh->ngx);
    
    // Matrix to be solved
    data.avec = cmatrix(laplace_maxmode + 1, mesh->ngx);
    data.bvec = cmatrix(laplace_maxmode + 1, mesh->ngx);
    data.cvec = cmatrix(laplace_maxmode + 1, mesh->ngx);
    
    // Working vectors
    data.v = cmatrix(laplace_maxmode + 1, mesh->ngx);
    data.w = cmatrix(laplace_maxmode + 1, mesh->ngx);

    // Result
    data.xk = cmatrix(laplace_maxmode + 1, mesh->ngx);

    // Communication buffers. Space for 2 complex values for each kz
    data.snd = new BoutReal[4*(laplace_maxmode+1)];
    data.rcv = new BoutReal[4*(laplace_maxmode+1)];

    data.y2i = new dcomplex[laplace_maxmode + 1];
  }

  /// Take FFTs of data
  static dcomplex *bk1d = NULL; ///< 1D in Z for taking FFTs

  if(bk1d == NULL)
    bk1d = new dcomplex[ncz/2 + 1];

  for(ix=0; ix < mesh->ngx; ix++) {
    ZFFT(b[ix], mesh->zShift[ix][data.jy], bk1d);
    for(kz = 0; kz <= laplace_maxmode; kz++)
      data.bk[kz][ix] = bk1d[kz];
  }

  /// Create the matrices to be inverted (one for each z point)

  /// Set matrix elements
  par_tridag_matrix(data.avec, data.bvec, data.cvec,
		    data.bk, data.jy, flags, a, ccoef, d);

  for(kz = 0; kz <= laplace_maxmode; kz++) {
    // Start PDD algorithm

    // Solve for xtilde, v and w (step 2)

    static dcomplex *e = NULL;
    if(e == NULL) {
      e = new dcomplex[mesh->ngx];
      for(ix=0;ix<mesh->ngx;ix++)
	e[ix] = 0.0;
    }

    dcomplex v0, x0; // Values to be sent to processor i-1

    if(mesh->firstX()) {
      // Domain includes inner boundary
      tridag(data.avec[kz], data.bvec[kz], data.cvec[kz], 
	     data.bk[kz], data.xk[kz], mesh->xend+1);
      
      // Add C (row m-1) from next processor
      
      e[mesh->xend] = data.cvec[kz][mesh->xend];
      tridag(data.avec[kz], data.bvec[kz], data.cvec[kz], 
	     e, data.w[kz], mesh->xend+1);

    }else if(mesh->lastX()) {
      // Domain includes outer boundary
      tridag(data.avec[kz]+mesh->xstart, 
	     data.bvec[kz]+mesh->xstart, 
	     data.cvec[kz]+mesh->xstart, 
	     data.bk[kz]+mesh->xstart, 
	     data.xk[kz]+mesh->xstart, 
	     mesh->xend - mesh->xend + 1);
      
      // Add A (row 0) from previous processor
      e[0] = data.avec[kz][mesh->xstart];
      tridag(data.avec[kz]+mesh->xstart, 
	     data.bvec[kz]+mesh->xstart, 
	     data.cvec[kz]+mesh->xstart, 
	     e, data.v[kz]+mesh->xstart,
	     mesh->xend+1);
      
      x0 = data.xk[kz][mesh->xstart];
      v0 = data.v[kz][mesh->xstart];

    }else {
      // No boundaries
      tridag(data.avec[kz]+mesh->xstart,
	     data.bvec[kz]+mesh->xstart,
	     data.cvec[kz]+mesh->xstart, 
	     data.bk[kz]+mesh->xstart, 
	     data.xk[kz]+mesh->xstart, 
	     mesh->xend - mesh->xstart + 1);

      // Add A (row 0) from previous processor
      e[0] = data.avec[kz][mesh->xstart];
      tridag(data.avec[kz]+mesh->xstart,
	     data.bvec[kz]+mesh->xstart,
	     data.cvec[kz]+mesh->xstart, 
	     e+mesh->xstart,
	     data.v[kz]+mesh->xstart,
	     mesh->xend - mesh->xstart + 1);
      e[0] = 0.0;
      
      // Add C (row m-1) from next processor
      e[mesh->xend] = data.cvec[kz][mesh->xend];
      tridag(data.avec[kz]+mesh->xstart,
	     data.bvec[kz]+mesh->xstart,
	     data.cvec[kz]+mesh->xstart, 
	     e+mesh->xstart,
	     data.v[kz]+mesh->xstart, 
	     mesh->xend - mesh->xstart + 1);
      e[mesh->xend] = 0.0;
    }
    
    // Put values into communication buffers
    data.snd[4*kz]   = x0.Real();
    data.snd[4*kz+1] = x0.Imag();
    data.snd[4*kz+2] = v0.Real();
    data.snd[4*kz+3] = v0.Imag();
  }
  
  // Stage 3: Communicate x0, v0 from node i to i-1
  
  if(!mesh->lastX()) {
    // All except the last processor expect to receive data
    // Post async receive
    data.recv_handle = mesh->irecvXOut(data.rcv, 4*(laplace_maxmode+1), PDD_COMM_XV);
  }

  if(!mesh->firstX()) {
    // Send the data
    
    mesh->sendXIn(data.snd, 4*(laplace_maxmode+1), PDD_COMM_XV);
  }

  return 0;
}

/// Middle part of the PDD algorithm
int invert_pdd_continue(PDD_data &data)
{
  // Wait for x0 and v0 to arrive from processor i+1
  
  if(!mesh->lastX()) {
    mesh->wait(data.recv_handle);

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
      
      data.y2i[kz] = (data.xk[kz][mesh->xend] - data.w[kz][mesh->xend]*x0) / (1. - data.w[kz][mesh->xend]*v0);
      
    }
  }
  
  if(!mesh->firstX()) {
    // All except pe=0 receive values from i-1. Posting async receive
    data.recv_handle = mesh->irecvXIn(data.rcv, 2*(laplace_maxmode+1), PDD_COMM_Y);
  }
  
  if(mesh->PE_XIND != (mesh->NXPE-1)) {
    // Send value to the (i+1)th processor
    
    for(int kz = 0; kz <= laplace_maxmode; kz++) {
      data.snd[2*kz]   = data.y2i[kz].Real();
      data.snd[2*kz+1] = data.y2i[kz].Imag();
    }
    
    mesh->sendXOut(data.snd, 2*(laplace_maxmode+1), PDD_COMM_Y);
  }
  
  return 0;
}

/// Last part of the PDD algorithm
int invert_pdd_finish(PDD_data &data, int flags, FieldPerp &x) {
  int ix, kz;

  x.allocate();
  x.setIndex(data.jy);
  
  if(!mesh->lastX()) {
    for(kz = 0; kz <= laplace_maxmode; kz++) {
      for(ix=0; ix < mesh->ngx; ix++)
	data.xk[kz][ix] -= data.w[kz][ix] * data.y2i[kz];
    }
  }

  if(!mesh->firstX()) {
    mesh->wait(data.recv_handle);
  
    for(kz = 0; kz <= laplace_maxmode; kz++) {
      dcomplex y2m = dcomplex(data.rcv[2*kz], data.rcv[2*kz+1]);
      
      for(ix=0; ix < mesh->ngx; ix++)
	data.xk[kz][ix] -= data.v[kz][ix] * y2m;
    }
  }
  
  // Have result in Fourier space. Convert back to BoutReal space

  static dcomplex *xk1d = NULL; ///< 1D in Z for taking FFTs

  int ncz = mesh->ngz-1;

  if(xk1d == NULL) {
    xk1d = new dcomplex[ncz/2 + 1];
    for(kz=0;kz<=ncz/2;kz++)
      xk1d[kz] = 0.0;
  }

  for(ix=0; ix<mesh->ngx; ix++){
    
    for(kz = 0; kz<= laplace_maxmode; kz++) {
      xk1d[kz] = data.xk[kz][ix];
    }

    if(flags & INVERT_ZERO_DC)
      xk1d[0] = 0.0;

    ZFFT_rev(xk1d, mesh->zShift[ix][data.jy], x[ix]);
    
    x[ix][ncz] = x[ix][0]; // enforce periodicity
  }

  return 0;
}

/**********************************************************************************
 *                              EXTERNAL INTERFACE
 **********************************************************************************/

/// Invert FieldPerp 
int invert_laplace(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d)
{
  if(mesh->NXPE == 1) {
    // Just use the serial code
    return invert_laplace_ser(b, x, flags, a, c, d);
  }else {
    // Parallel inversion using PDD

    if(invert_use_pdd) {
      static PDD_data data;
      static bool allocated = false;
      if(!allocated) {
	data.bk = NULL;
	allocated = true;
      }
      
      invert_pdd_start(b, flags, a, data, c, d);
      invert_pdd_continue(data);
      invert_pdd_finish(data, flags, x);
    }else {
      static SPT_data data;
      static bool allocated = false;
      if(!allocated) {
	data.bk = NULL;
        data.comm_tag = SPT_DATA;
	allocated = true;
      }
      
      invert_spt_start(b, flags, a, data, c, d);
      invert_spt_finish(data, flags, x);
    }
  }

  return 0;
}

/// Extracts perpendicular slices from 3D fields and inverts separately
/*!
 * In parallel (mesh->NXPE > 1) this tries to overlap computation and communication.
 * This is done at the expense of more memory useage. Setting low_mem
 * in the config file uses less memory, and less communication overlap
 */
int invert_laplace(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {
  int jy, jy2;
  FieldPerp xperp;
  int ret;
  
  Timer timer("invert");
  
  x.allocate();

  int ys = mesh->ystart, ye = mesh->yend;
  
  if(MYPE_IN_CORE == 0) {
    // NOTE: REFINE THIS TO ONLY SOLVE IN BOUNDARY Y CELLS
    ys = 0;
    ye = mesh->ngy-1;
  }
  
  if((mesh->NXPE == 1) || invert_low_mem) {
    
    for(jy=ys; jy <= ye; jy++) {
      if((flags & INVERT_IN_SET) || (flags & INVERT_OUT_SET))
	xperp = x.slice(jy); // Using boundary values
      
      if((ret = invert_laplace(b.slice(jy), xperp, flags, a, c, d)))
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
	invert_pdd_start(b.slice(jy), flags, a, data[jy], c, d);
      
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
	for(jy=ys;jy<=ye;jy++) {
	  data[jy].bk = NULL; // Mark as unallocated for PDD routine
          data[jy].comm_tag = SPT_DATA + jy; // Give each one a different tag
        }
      }
      
      for(jy=ys; jy <= ye; jy++) {	
	// And start another one going
	invert_spt_start(b.slice(jy), flags, a, data[jy], c, d);

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

  x.setLocation(b.getLocation());

  return 0;
}

const Field3D invert_laplace(const Field3D &b, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {
  Field3D x;
  
  invert_laplace(b, x, flags, a, c, d);
  return x;
}

