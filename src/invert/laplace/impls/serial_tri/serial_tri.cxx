/**************************************************************************
 * Perpendicular Laplacian inversion. Serial code using FFT
 * and tridiagonal solver.
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

#include "globals.hxx"
#include "serial_tri.hxx"

#include <boutexception.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <lapack_routines.hxx>
#include <bout/constants.hxx>
#include <cmath>

#include <output.hxx>

LaplaceSerialTri::LaplaceSerialTri(Options *opt) : Laplacian(opt), A(0.0), C(1.0), D(1.0) {
  
  if(!mesh->firstX() || !mesh->lastX()) {
    throw BoutException("LaplaceSerialTri only works for mesh->NXPE = 1");
  }
  
  // Allocate memory

  int ncz = mesh->ngz-1;

  bk = cmatrix(mesh->ngx, ncz/2 + 1);
  bk1d = new dcomplex[mesh->ngx];
  
  xk = cmatrix(mesh->ngx, ncz/2 + 1);
  xk1d = new dcomplex[mesh->ngx];
  
  avec = new dcomplex[mesh->ngx];
  bvec = new dcomplex[mesh->ngx];
  cvec = new dcomplex[mesh->ngx];
}

LaplaceSerialTri::~LaplaceSerialTri() {
  free_cmatrix(bk);
  delete[] bk1d;
  free_cmatrix(xk);
  delete[] xk1d;
  
  delete[] avec;
  delete[] bvec;
  delete[] cvec;
}

const FieldPerp LaplaceSerialTri::solve(const FieldPerp &b) {
  return solve(b,b);
}

const FieldPerp LaplaceSerialTri::solve(const FieldPerp &b, const FieldPerp &x0) {
  FieldPerp x;
  x.allocate();

  int jy = b.getIndex();
  x.setIndex(jy);
  
  int ncz = mesh->ngz-1;
  int ncx = mesh->ngx-1;
  
  int inbndry = 2, outbndry=2;
  
  if(flags & INVERT_BNDRY_ONE) {
    inbndry = outbndry = 1;
  }
  if(flags & INVERT_BNDRY_IN_ONE)
    inbndry = 1;
  if(flags & INVERT_BNDRY_OUT_ONE)
    outbndry = 1;
  
  #pragma omp parallel for
  for(int ix=0;ix<mesh->ngx;ix++) {
    // for fixed ix,jy set a complex vector rho(z)
    
    if(((ix < inbndry) && (flags & INVERT_IN_SET)) ||
       ((ncx-ix < outbndry) && (flags & INVERT_OUT_SET))) {
      // Use the values in x0 in the boundary
      ZFFT(x0[ix], mesh->zShift[ix][jy], bk[ix]);
      
    }else
      ZFFT(b[ix], mesh->zShift[ix][jy], bk[ix]);
  }

  for(int iz=0;iz<=ncz/2;iz++) {
    // solve differential equation in x

    // set bk1d
    BoutReal flt;
    if (iz>maxmode) flt=0.0; else flt=1.0;
      
    for(int ix=0;ix<=ncx;ix++)
      bk1d[ix] = bk[ix][iz] * flt;
    
    ///////// PERFORM INVERSION /////////
      
    for(int ix=inbndry;ix<=ncx-outbndry;ix++) {
      tridagCoefs(ix, jy, iz, avec[ix], bvec[ix], cvec[ix], &C, &D);
      
      bvec[ix] += A[ix][jy];
    }

    if(!mesh->periodicX) {
      // Need boundary conditions
      /// By default, set RHS to zero, unless INVERT_*_RHS set
      if(!(flags & (INVERT_IN_RHS | INVERT_IN_SET))) {
        for(int ix=0;ix<inbndry;ix++)
          bk1d[ix] = 0.;
      }
      if(!(flags & (INVERT_OUT_RHS | INVERT_OUT_SET))) {
        for(int ix=mesh->ngx-outbndry;ix<mesh->ngx;ix++)
          bk1d[ix] = 0.;
      }
	
      // Set boundary conditions
	
      if(iz == 0) {
        // DC
	  
        // Inner boundary
        if(flags & INVERT_DC_IN_GRAD) {
          // Zero gradient at inner boundary
	    
          if((flags & INVERT_IN_SYM) && (inbndry > 1) && mesh->BoundaryOnCell) {
            // Use symmetric boundary to set zero-gradient
	      
            for (int ix=0;ix<inbndry-1;ix++) {
              avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
            }
            // Symmetric on last point
            avec[inbndry-1] = 1.0; bvec[inbndry-1] = 0.0; cvec[inbndry-1] = -1.0;
          }else {
            for (int ix=0;ix<inbndry;ix++){
              avec[ix]=dcomplex(0.0,0.0);
              bvec[ix]=dcomplex(1.,0.);
              cvec[ix]=dcomplex(-1.,0.);
            }
          }
        }else if(flags & INVERT_DC_IN_GRADPAR) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] =  0.0;
            bvec[ix] =  1.0/sqrt(mesh->g_22(ix,jy));
            cvec[ix] = -1.0/sqrt(mesh->g_22(ix+1,jy));
          }
        }else if(flags & INVERT_DC_IN_GRADPARINV) {
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] =  0.0;
            bvec[ix] =  sqrt(mesh->g_22(ix,jy));
            cvec[ix] = -sqrt(mesh->g_22(ix+1,jy));
          }
        }else if(flags & INVERT_IN_SET) {
          for(int ix=0;ix<inbndry;ix++) {
            avec[ix] = 0.0;
            bvec[ix] = 1.0;
            cvec[ix] = 0.0;
          }
        }else if (flags & INVERT_DC_IN_LAP) {
          // Decaying boundary conditions
          BoutReal ksq = -(A(inbndry, jy));
          if(ksq < 0.0)
            throw BoutException("ksq must be positive");
          BoutReal k = sqrt(ksq);
          //output << "k = " << k << endl;
          for(int ix=0;ix<inbndry;ix++){
            avec[ix] =  0.;
            bvec[ix] =  1.;
            //output << "factor = " << exp(-k*mesh->dx(ix,jy)/sqrt(mesh->g11(ix,jy))) << endl;
            cvec[ix] = -exp(-k*mesh->dx(ix,jy)/sqrt(mesh->g11(ix,jy)));
          }
          
        }else {
          // Zero value at inner boundary
          if(flags & INVERT_IN_SYM) {
            // Use anti-symmetric boundary to set zero-value
	    
            // Zero-gradient for first point(s)
            for(int ix=0;ix<inbndry-1;ix++) {
              avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
            }

            if(mesh->BoundaryOnCell) {
              // Antisymmetric about boundary on cell
              avec[inbndry-1]=1.0; bvec[inbndry-1]=0.0; cvec[inbndry-1]= 1.0;
            }else { 
              // Antisymmetric across boundary between cells
              avec[inbndry-1]=0.0; bvec[inbndry-1]=1.0; cvec[inbndry-1]= 1.0;
            }
	    
          }else {
            for (int ix=0;ix<inbndry;ix++){
              avec[ix]=dcomplex(0.,0.);
              bvec[ix]=dcomplex(1.,0.);
              cvec[ix]=dcomplex(0.,0.);
            }
          }
        }
	
        // Outer boundary
        if(flags & INVERT_DC_OUT_GRAD) {
          // Zero gradient at outer boundary

          if((flags & INVERT_OUT_SYM) && (outbndry > 1) && mesh->BoundaryOnCell) {
            // Use symmetric boundary to set zero-gradient
	    
            for (int ix=0;ix<outbndry-1;ix++) {
              avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
            }
            // Symmetric on last point
            int ix = outbndry-1;
            avec[ncx-ix] = 1.0; bvec[ncx-ix] = 0.0; cvec[ncx-ix] = -1.0;
	    
          }else {
            for (int ix=0;ix<outbndry;ix++){
              cvec[ncx-ix]=dcomplex(0.,0.);
              bvec[ncx-ix]=dcomplex(1.,0.);
              avec[ncx-ix]=dcomplex(-1.,0.);
            }
          }
        }else if(flags & INVERT_OUT_SET) {
          // Setting the values in the outer boundary
          for(int ix=0;ix<outbndry;ix++) {
            avec[ncx-ix] = 0.0;
            bvec[ncx-ix] = 1.0;
            cvec[ncx-ix] = 0.0;
          }
        }else {
          // Zero value at outer boundary
          if(flags & INVERT_OUT_SYM) {
            // Use anti-symmetric boundary to set zero-value
	    
            // Zero-gradient for first point(s)
            for(int ix=0;ix<outbndry-1;ix++) {
              avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
            }
            int ix = outbndry-1;
            if(mesh->BoundaryOnCell) {
              // Antisymmetric about boundary on cell
              avec[ncx-ix]=1.0; bvec[ncx-ix]=0.0; cvec[ncx-ix]= 1.0;
            }else { 
              // Antisymmetric across boundary between cells
              avec[ncx-ix]=1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
            }
          }else {
            for (int ix=0;ix<outbndry;ix++){
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
	  
          if((flags & INVERT_IN_SYM) && (inbndry > 1) && mesh->BoundaryOnCell) {
            // Use symmetric boundary to set zero-gradient
	    
            for (int ix=0;ix<inbndry-1;ix++) {
              avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
            }
            // Symmetric on last point
            avec[inbndry-1] = 1.0; bvec[inbndry-1] = 0.0; cvec[inbndry-1] = -1.0;
          }else {
            for (int ix=0;ix<inbndry;ix++){
              avec[ix]=dcomplex(0.,0.);
              bvec[ix]=dcomplex(1.,0.);
              cvec[ix]=dcomplex(-1.,0.);
            }
          }
        }else if(flags & INVERT_IN_SET) {
          // Setting the values in the boundary
          for(int ix=0;ix<inbndry;ix++) {
            avec[ix] = 0.0;
            bvec[ix] = 1.0;
            cvec[ix] = 0.0;
          }
        }else if(flags & INVERT_AC_IN_LAP) {
          // Use decaying zero-Laplacian solution in the boundary
          BoutReal kwave=iz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
          for (int ix=0;ix<inbndry;ix++) {
            avec[ix] = 0.0;
            bvec[ix] = -1.0;
            cvec[ix] = exp(-1.0*sqrt(mesh->g33[ix][jy]/mesh->g11[ix][jy])*kwave*mesh->dx[ix][jy]);
          }
        }else {
          // Zero value at inner boundary

          if(flags & INVERT_IN_SYM) {
            // Use anti-symmetric boundary to set zero-value
	    
            // Zero-gradient for first point(s)
            for(int ix=0;ix<inbndry-1;ix++) {
              avec[ix]=0.0; bvec[ix]=1.0; cvec[ix]= -1.0;
            }

            if(mesh->BoundaryOnCell) {
              // Antisymmetric about boundary on cell
              avec[inbndry-1]=1.0; bvec[inbndry-1]=0.0; cvec[inbndry-1]= 1.0;
            }else { 
              // Antisymmetric across boundary between cells
              avec[inbndry-1]=0.0; bvec[inbndry-1]=1.0; cvec[inbndry-1]= 1.0;
            }
	    
          }else {
            for (int ix=0;ix<inbndry;ix++){
              avec[ix]=dcomplex(0.,0.);
              bvec[ix]=dcomplex(1.,0.);
              cvec[ix]=dcomplex(0.,0.);
            }
          }
        }
	
        // Outer boundary
        if(flags & INVERT_AC_OUT_GRAD) {
          // Zero gradient at outer boundary
	  
          if((flags & INVERT_OUT_SYM) && (outbndry > 1) && mesh->BoundaryOnCell) {
            // Use symmetric boundary to set zero-gradient
	    
            for (int ix=0;ix<outbndry-1;ix++) {
              avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
            }
            // Symmetric on last point
            int ix = outbndry-1;
            avec[ncx-ix] = 1.0; bvec[ncx-ix] = 0.0; cvec[ncx-ix] = -1.0;
	    
          }else {
            for (int ix=0;ix<outbndry;ix++){
              cvec[ncx-ix]=dcomplex(0.,0.);
              bvec[ncx-ix]=dcomplex(1.,0.);
              avec[ncx-ix]=dcomplex(-1.,0.);
            }
          }
        }else if(flags & INVERT_AC_OUT_LAP) {
          // Use decaying zero-Laplacian solution in the boundary
          BoutReal kwave=iz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
          for (int ix=0;ix<outbndry;ix++) {
            avec[ncx-ix] = exp(-1.0*sqrt(mesh->g33[ncx-ix][jy]/mesh->g11[ncx-ix][jy])*kwave*mesh->dx[ncx-ix][jy]);;
            bvec[ncx-ix] = -1.0;
            cvec[ncx-ix] = 0.0;
          }
        }else if(flags & INVERT_OUT_SET) {
          // Setting the values in the outer boundary
          for(int ix=0;ix<outbndry;ix++) {
            avec[ncx-ix] = 0.0;
            bvec[ncx-ix] = 1.0;
            cvec[ncx-ix] = 0.0;
          }
        }else {
          // Zero value at outer boundary

          if(flags & INVERT_OUT_SYM) {
            // Use anti-symmetric boundary to set zero-value
	    
            // Zero-gradient for first point(s)
            for(int ix=0;ix<outbndry-1;ix++) {
              avec[ncx-ix]=-1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
            }
            int ix = outbndry-1;
            if(mesh->BoundaryOnCell) {
              // Antisymmetric about boundary on cell
              avec[ncx-ix]=1.0; bvec[ncx-ix]=0.0; cvec[ncx-ix]= 1.0;
            }else {
              // Antisymmetric across boundary between cells
              avec[ncx-ix]=1.0; bvec[ncx-ix]=1.0; cvec[ncx-ix]= 0.0;
            }
          }else {
            for (int ix=0;ix<outbndry;ix++){
              cvec[ncx-ix]=dcomplex(0.,0.);
              bvec[ncx-ix]=dcomplex(1.,0.);
              avec[ncx-ix]=dcomplex(0.,0.);
            }
          }
        }
      }
      
      // Call tridiagonal solver
      tridag(avec, bvec, cvec, bk1d, xk1d, mesh->ngx);
      
      if((flags & INVERT_IN_SYM) && (inbndry > 1)) {
        // (Anti-)symmetry on inner boundary. Nothing to do if only one boundary cell
        int xloc = 2*inbndry;
        if(!mesh->BoundaryOnCell)
          xloc--;
	  
        if( ((iz == 0) && (flags & INVERT_DC_IN_GRAD)) || ((iz != 0) && (flags & INVERT_AC_IN_GRAD)) ) {
          // Inner gradient zero - symmetric
          for(int ix=0;ix<inbndry-1;ix++)
            xk1d[ix] = xk1d[xloc-ix];
        }else {
          // Inner value zero - antisymmetric
          for(int ix=0;ix<inbndry-1;ix++)
            xk1d[ix] = -xk1d[xloc-ix];
        }
      }
      if((flags & INVERT_OUT_SYM) && (outbndry > 1)) {
        // (Anti-)symmetry on outer boundary. Nothing to do if only one boundary cell
	  
        int xloc =  mesh->ngx - 2*outbndry;
        if(mesh->BoundaryOnCell)
          xloc--;
	
        if( ((iz == 0) && (flags & INVERT_DC_IN_GRAD)) || ((iz != 0) && (flags & INVERT_AC_IN_GRAD)) ) {
          // Outer gradient zero - symmetric
          for(int ix=0;ix<outbndry-1;ix++)
            xk1d[ncx-ix] = xk1d[xloc + ix];
        }else {
          // Outer value zero - antisymmetric
          for(int ix=0;ix<outbndry-1;ix++)
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

  // Done inversion, transform back
  
  for(int ix=0; ix<=ncx; ix++){
    
    if(flags & INVERT_ZERO_DC)
      xk[ix][0] = 0.0;
    
    ZFFT_rev(xk[ix], mesh->zShift[ix][jy], x[ix]);
    
    x[ix][mesh->ngz-1] = x[ix][0]; // enforce periodicity
    
    for(int kz=0;kz<mesh->ngz;kz++)
      if(!finite(x[ix][kz]))
        throw BoutException("Non-finite at %d, %d, %d", ix, jy, kz);
  }

  return x;
}
