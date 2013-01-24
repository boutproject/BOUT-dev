/************************************************************************
 * Inversion of parallel derivatives
 * 
 * Inverts a matrix of the form 
 *
 * A + B * Grad2_par2
 * 
 * SERIAL ALGORITHM, for testing only
 *
 * Author: Ben Dudson, University of York, Oct 2011
 * 
 * Known issues:
 * ------------
 *
 * - Only works for NPE = 1
 * - Assumes all flux surfaces closed
 * - Coefficients A and B must be 2D (X,Y)
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
 ************************************************************************/

#include <globals.hxx>
#include <utils.hxx>
#include "serial.hxx"
#include <fft.hxx>
#include <lapack_routines.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>

#include <bout/surfaceiter.hxx>

#include <cmath>

InvertParSerial::InvertParSerial() {
  rhs = cmatrix(mesh->ngy, (mesh->ngz-1)/2 + 1);
  rhsk = new dcomplex[mesh->ngy-4];
  xk = new dcomplex[mesh->ngy-4];
  a = new dcomplex[mesh->ngy-4];
  b = new dcomplex[mesh->ngy-4];
  c = new dcomplex[mesh->ngy-4];
}

InvertParSerial::~InvertParSerial() {
  free_cmatrix(rhs);
  delete[] rhsk;
  delete[] xk;
  delete[] a;
  delete[] b;
  delete[] c;
}

const Field3D InvertParSerial::solve(const Field3D &f) {
#ifdef CHECK
  msg_stack.push("InvertParSerial::solve(Field3D)");
#endif
  
  Field3D result;
  result.allocate();
  
  // Loop over flux-surfaces
  SurfaceIter surf(mesh);
  for(surf.first(); !surf.isDone(); surf.next()) {
    int x = surf.xpos;
    BoutReal ts; // Twist-shift angle
    if(!surf.closed(ts))
      throw BoutException("InvertParSerial doesn't handle open surfaces");
    
    // Take Fourier transform 
    for(int y=0;y<mesh->ngy-4;y++)
      rfft(f[x][y+2], mesh->ngz-1, rhs[y]);
    
    // Set up tridiagonal system. Same for all k
    // except for phase shift at twist-shift
    for(int y=0;y<mesh->ngy-4;y++) {
      BoutReal acoef = A[x][y+2];
      BoutReal bcoef = B[x][y+2] / (mesh->g_22[x][y+2] * SQ(mesh->dy[x][y+2]));
      
      a[y] =            bcoef;
      b[y] = acoef - 2.*bcoef;
      c[y] =            bcoef;
    }
    
    dcomplex a0 = a[0]; // Will be modified later
    dcomplex cn = c[mesh->ngy-5];
    
    // Solve cyclic tridiagonal system for each k
    int nyq = (mesh->ngz-1)/2;
    for(int k=0;k<=nyq;k++) {
      // Copy component of rhs into 1D array
      for(int y=0;y<mesh->ngy-4;y++)
        rhsk[y] = rhs[y][k];
      
      // Modify coefficients across twist-shift
      BoutReal kwave=k*2.0*PI/mesh->zlength; // wave number is 1/[rad]
      dcomplex phase(cos(kwave*ts) , -sin(kwave*ts));
      a[0] = a0 * phase;
      c[mesh->ngy-5] = cn / phase;
      
      // Solve cyclic tridiagonal system
      cyclic_tridag(a, b, c, rhsk, xk, mesh->ngy-4);
      
      // Put back into rhs array
      for(int y=0;y<mesh->ngy-4;y++)
        rhs[y][k] = xk[y];
    }
    
    // Inverse Fourier transform 
    for(int y=0;y<mesh->ngy-4;y++)
      irfft(rhs[y], mesh->ngz-1, result[x][y+2]);
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}

