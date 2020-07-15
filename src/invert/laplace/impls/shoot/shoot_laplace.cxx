/*!
 * \file shoot_laplace.cxx
 *
 * \brief Laplacian solver using shooting method
 *  
 * CHANGELOG
 * =========
 * 
 * Feb 2014: Ben Dudson <benjamin.dudson@york.ac.uk>
 *         * Initial version
 * 
 **************************************************************************
 * Copyright 2014 B.D.Dudson
 *
 * Contact: Ben Dudson, benjamin.dudson@york.ac.uk
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

#include "shoot_laplace.hxx"
#include <bout/mesh.hxx>
#include <globals.hxx>
#include <fft.hxx>
#include <bout/constants.hxx>

LaplaceShoot::LaplaceShoot(Options *opt, const CELL_LOC loc, Mesh *mesh_in)
    : Laplacian(opt, loc, mesh_in), Acoef(0.0), Ccoef(1.0), Dcoef(1.0) {
  throw BoutException("LaplaceShoot is a test implementation and does not currently work. Please select a different implementation.");

  Acoef.setLocation(location);
  Ccoef.setLocation(location);
  Dcoef.setLocation(location);

  if(localmesh->periodicX) {
        throw BoutException("LaplaceShoot does not work with periodicity in the x direction (localmesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
  }

	
  nmode = maxmode + 1; // Number of Z modes. maxmode set in invert_laplace.cxx from options
  
  // Allocate memory
  int size = (localmesh->LocalNz)/2 + 1;
  km.reallocate(size);
  kc.reallocate(size);
  kp.reallocate(size);

  for(int i=0;i<size;i++) {
    km[i] = 0.0;
    kc[i] = 0.0;
    kp[i] = 0.0;
  }

  rhsk.reallocate(size);

  buffer.reallocate(4 * maxmode);
}

FieldPerp LaplaceShoot::solve(const FieldPerp& rhs) {
  ASSERT1(localmesh == rhs.getMesh());
  ASSERT1(rhs.getLocation() == location);

  FieldPerp x{emptyFrom(rhs)}; // Result
  
  int jy = rhs.getIndex();  // Get the Y index

  // Get the width of the boundary
  
  int inbndry = localmesh->xstart, outbndry=localmesh->xstart;
  // If the flags to assign that only one guard cell should be used is set
  if((global_flags & INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2))  {
    inbndry = outbndry = 1;
  }
  if(inner_boundary_flags & INVERT_BNDRY_ONE)
    inbndry = 1;
  if(outer_boundary_flags & INVERT_BNDRY_ONE)
    outbndry = 1;
  
  int xs, xe;
  xs = localmesh->xstart; // Starting X index
  if(localmesh->firstX())
    xs = inbndry;
  xe = localmesh->xend;  // Last X index
  if(localmesh->lastX())
    xe = localmesh->LocalNx-outbndry-1;

  if(localmesh->lastX()) {
    // Set initial value and gradient to zero
    // by setting kc and kp
    
    for(int i=0;i<maxmode;i++) {
      kc[i] = 0.0;
      kp[i] = 0.0;
    }
    
    for(int ix=xe;ix<localmesh->LocalNx;ix++)
      for(int iz=0;iz<localmesh->LocalNz;iz++) {
        x(ix, iz) = 0.0;
      }
      
  }else {
    // Wait for processor outer X
    comm_handle handle = localmesh->irecvXOut(std::begin(buffer), 4 * maxmode, jy);
    localmesh->wait(handle);
    
    // Copy into kc, kp
    for(int i=0;i<maxmode;i++) {
      kc[i] = dcomplex(buffer[4*i], buffer[4*i+1]);
      kp[i] = dcomplex(buffer[4*i+2], buffer[4*i+3]);
    }
    
    // Calculate solution at xe using kc
    irfft(std::begin(kc), localmesh->LocalNz, x[xe]);
  }
  
  // kc and kp now set to result at x and x+1 respectively
  // Use b at x to get km at x-1
  // Loop inwards from edge
  const BoutReal zlength = getUniform(coords->zlength());
  for(int ix=xe; ix >= xs; ix--) {
    rfft(rhs[ix], localmesh->LocalNz, std::begin(rhsk));

    for(int kz=0; kz<maxmode; kz++) {
      BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]

      // Get the coefficients
      dcomplex a,b,c;
      tridagCoefs(ix, jy, kwave, a, b, c, &Ccoef, &Dcoef);
      b += Acoef(ix, jy);

      // a*km + b*kc + c*kp = rhsk
      
      km[kz] = (rhsk[kz] - b*kc[kz] - c*kp[kz]) / a;
    }
    
    // Inverse FFT to get x[ix-1]
    irfft(std::begin(km), localmesh->LocalNz, x[ix - 1]);

    // Cycle km->kc->kp
    std::swap(kp, kc);
    std::swap(kc, km);
  }
  
  // Finished on this processor. Send data to next inner processor
  if(!localmesh->firstX()) {
    // Should be able to send dcomplex buffers. For now copy into BoutReal buffer
    for(int i=0;i<maxmode;i++) {
      buffer[4*i]     = kc[i].real();
      buffer[4*i + 1] = kc[i].imag();
      buffer[4*i + 2] = kp[i].real();
      buffer[4*i + 3] = kp[i].imag();
    }
    localmesh->sendXIn(std::begin(buffer), 4 * maxmode, jy);
  }else {
    // Set inner boundary
    for(int ix=xs-2;ix>=0;ix--) {
      for(int iz=0;iz<localmesh->LocalNz;iz++) {
        x(ix, iz) = x(xs - 1, iz);
      }
    }
  }
  
  return x;
}
