/**************************************************************************
 * Perpendicular Laplacian inversion. Serial code using FFT
 * and band-diagonal solver
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
#include "serial_band.hxx"

#include <fft.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <lapack_routines.hxx>
#include <bout/constants.hxx>

#include <output.hxx>

//#define SECONDORDER // Define to use 2nd order differencing

LaplaceSerialBand::LaplaceSerialBand(Options *opt) : Laplacian(opt), Acoef(0.0), Ccoef(1.0), Dcoef(1.0) {
  if(!mesh->firstX() || !mesh->lastX())
    throw BoutException("LaplaceSerialBand only works for mesh->NXPE = 1");
  if(mesh->periodicX)
    throw BoutException("LaplaceSerialBand doesn't work with periodic X");
  
  // Allocate memory

  int ncz = mesh->ngz-1;
  bk = cmatrix(mesh->ngx, ncz/2 + 1);
  bk1d = new dcomplex[mesh->ngx];
  
  xk = cmatrix(mesh->ngx, ncz/2 + 1);
  xk1d = new dcomplex[mesh->ngx];
  
  A = cmatrix(mesh->ngx, 5);
}

LaplaceSerialBand::~LaplaceSerialBand() {
  free_cmatrix(bk);
  delete[] bk1d;
  free_cmatrix(xk);
  delete[] xk1d;
  
  free_cmatrix(A);
}

const FieldPerp LaplaceSerialBand::solve(const FieldPerp &b) {
  return solve(b,b);
}

const FieldPerp LaplaceSerialBand::solve(const FieldPerp &b, const FieldPerp &x0) {
  FieldPerp x;
  x.allocate();

  int jy = b.getIndex();
  x.setIndex(jy);
  
  int ncz = mesh->ngz-1;
  int ncx = mesh->ngx-1;

  int xbndry = 2; // Width of the x boundary
  if(flags & INVERT_BNDRY_ONE)
    xbndry = 1;

  #pragma omp parallel for
  for(int ix=0;ix<mesh->ngx;ix++) {
    // for fixed ix,jy set a complex vector rho(z)
    
    if(((ix < xbndry) && (flags & INVERT_IN_SET)) ||
       ((ncx-ix < xbndry) && (flags & INVERT_OUT_SET))) {
      // Use the values in x0 in the boundary
      ZFFT(x0[ix], mesh->zShift[ix][jy], bk[ix]);
    }else
      ZFFT(b[ix], mesh->zShift[ix][jy], bk[ix]);
  }
  
  int xstart, xend;
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
      
    if (iz>maxmode) flt=0.0; else flt=1.0;

    // set bk1d
    for(int ix=0;ix<mesh->ngx;ix++)
      bk1d[ix] = bk[ix][iz]*flt;

    // Fill in interior points

    for(int ix=xstart;ix<=xend;ix++) {
#ifdef SECONDORDER 
      // Use second-order differencing. Useful for testing the tridiagonal solver
      // with different boundary conditions
      dcomplex a,b,c;
      tridagCoefs(ix, jy, iz, a, b, c, &Ccoef, &Dcoef);
      
      A[ix][0] = 0.;
      A[ix][1] = a;
      A[ix][2] = b + Acoef[ix][jy];
      A[ix][3] = c;
      A[ix][4] = 0.;
#else
      // Set coefficients
      coef1 = mesh->g11[ix][jy];  // X 2nd derivative
      coef2 = mesh->g33[ix][jy];  // Z 2nd derivative
      coef3 = mesh->g13[ix][jy];  // X-Z mixed derivatives
      coef4 = 0.0;          // X 1st derivative
      coef5 = 0.0;          // Z 1st derivative
      coef6 = 0.0;          // Constant

      // Multiply Delp2 component by a factor
      coef1 *= Dcoef[ix][jy];
      coef2 *= Dcoef[ix][jy];
      coef3 *= Dcoef[ix][jy];

      coef6 = Acoef[ix][jy];
	
      if(all_terms) {
        coef4 = mesh->G1[ix][jy];
        coef5 = mesh->G3[ix][jy];
      }

      if(nonuniform) {
        // non-uniform mesh correction
        if((ix != 0) && (ix != ncx))
          coef4 += mesh->g11[ix][jy]*( (1.0/mesh->dx[ix+1][jy]) - (1.0/mesh->dx[ix-1][jy]) )/(2.0*mesh->dx[ix][jy]);
      }

      // A first order derivative term (1/c)\nabla_perp c\cdot\nabla_\perp x
    
      if((ix > 1) && (ix < (mesh->ngx-2)))
        coef4 += mesh->g11[ix][jy] * (Ccoef[ix-2][jy] - 8.*Ccoef[ix-1][jy] + 8.*Ccoef[ix+1][jy] - Ccoef[ix+2][jy]) / (12.*mesh->dx[ix][jy]*(Ccoef[ix][jy]));

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
#endif
    }

    if(xbndry < 2) {
      // Use 2nd order near edges

      int ix = 1;

      coef1=mesh->g11[ix][jy]/(SQ(mesh->dx[ix][jy]));
      coef2=mesh->g33[ix][jy];
      coef3= kwave * mesh->g13[ix][jy]/(2. * mesh->dx[ix][jy]);
        
      // Multiply Delp2 component by a factor
      coef1 *= Dcoef[ix][jy];
      coef2 *= Dcoef[ix][jy];
      coef3 *= Dcoef[ix][jy];
        
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

    if(iz == 0) {
      // DC
	
      // Inner boundary
      if(flags & INVERT_DC_IN_GRAD) {
        // Zero gradient at inner boundary. 2nd-order accurate
        // One-sided differences
        for (int ix=0;ix<xbndry;ix++) {
          A[ix][0] =  0.;
          A[ix][1] =  0.;
          A[ix][2] =  -3.;
          A[ix][3] =  4.;
          A[ix][4] =  -1.;
        }
        
      }else if(flags & INVERT_DC_IN_GRADPAR) {
        for (int ix=0;ix<xbndry;ix++) {
          A[ix][0] =  0.;
          A[ix][1] =  0.;
          A[ix][2] =  -3./sqrt(mesh->g_22(ix,jy));
          A[ix][3] =  4./sqrt(mesh->g_22(ix+1,jy));
          A[ix][4] =  -1./sqrt(mesh->g_22(ix+2,jy));
        }
      }else if(flags & INVERT_DC_IN_GRADPARINV) {
        for (int ix=0;ix<xbndry;ix++) {
          A[ix][0] =  0.;
          A[ix][1] =  0.;
          A[ix][2] =  -3.*sqrt(mesh->g_22(ix,jy));
          A[ix][3] =  4.*sqrt(mesh->g_22(ix+1,jy));
          A[ix][4] =  -sqrt(mesh->g_22(ix+2,jy));
        }
      }else if (flags & INVERT_DC_IN_LAP) {
        for (int ix=0;ix<xbndry;ix++) {
          A[ix][0] = 0.;
          A[ix][1] = 0.;
          A[ix][2] = 1.;
          A[ix][3] = -2;
          A[ix][4] = 1.;
        }
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
        
        coef4 = Acoef[ix][jy];
	  
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
        
        coef4 = Acoef[ix][jy];
	  
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
  
  // Done inversion, transform back

  for(int ix=0; ix<=ncx; ix++){
    if(flags & INVERT_ZERO_DC)
      xk[ix][0] = 0.0;

    ZFFT_rev(xk[ix], mesh->zShift[ix][jy], x[ix]);
    
    x[ix][mesh->ngz-1] = x[ix][0]; // enforce periodicity
  }

  return x;
}
