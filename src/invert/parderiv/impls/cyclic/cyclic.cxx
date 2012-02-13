/************************************************************************
 * Inversion of parallel derivatives
 * 
 * Inverts a matrix of the form 
 *
 * A + B * Grad2_par2
 * 
 * Parallel algorithm, using Cyclic Reduction
 *
 * Author: Ben Dudson, University of York, Oct 2011
 * 
 * Known issues:
 * ------------
 *
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
#include "cyclic.hxx"
#include <fft.hxx>
#include <boutexception.hxx>
#include <cyclic_reduction.hxx>

#include <cmath>

InvertParCR::InvertParCR() {
  // Number of k equations to solve for each x location
  nsys = 1 + (mesh->ngz-1)/2; 

  rhs = cmatrix(mesh->ngy, nsys);
  
  rhsk = cmatrix(nsys, mesh->ngy-4);
  xk = cmatrix(nsys, mesh->ngy-4);
  a = cmatrix(nsys, mesh->ngy-4);
  b = cmatrix(nsys, mesh->ngy-4);
  c = cmatrix(nsys, mesh->ngy-4);
}

InvertParCR::~InvertParCR() {
  free_cmatrix(rhs);
  delete[] rhsk;
  delete[] xk;
  delete[] a;
  delete[] b;
  delete[] c;
}

const Field3D InvertParCR::solve(const Field3D &f) {
#ifdef CHECK
  msg_stack.push("InvertParCR::solve(Field3D)");
#endif
  
  Field3D result;
  result.allocate();
  
  
  // Create cyclic reduction object
  CyclicReduce<dcomplex> *cr = 
    new CyclicReduce<dcomplex>();
  
  // Loop over flux-surfaces
  SurfaceIter *surf = mesh->iterateSurfaces();
  for(surf->first(); !surf->isDone(); surf->next()) {
    int x = surf->xpos;
    
    // Setup CyclicReduce object
    cr->setup(surf->communicator(), mesh->ngy-4);

    // Check if surface is periodic
    BoutReal ts; // Twist-shift angle
    cr->setPeriodic(surf->closed(ts));
    
    // Take Fourier transform 
    for(int y=0;y<mesh->ngy-4;y++)
      rfft(f[x][y+2], mesh->ngz-1, rhs[y]);
    
    // Set up tridiagonal system. Same for all k
    for(int k=0; k<nsys; k++) {
      for(int y=0;y<mesh->ngy-4;y++) {
	BoutReal acoef = A[x][y+2];
	BoutReal bcoef = B[x][y+2] / (mesh->g_22[x][y+2] * SQ(mesh->dy[x][y+2]));
	
	a[k][y] =            bcoef;
	b[k][y] = acoef - 2.*bcoef;
	c[k][y] =            bcoef;
	
	rhsk[k][y] = rhs[y][k]; // Transpose
      }
    }
    // Twist-shift
    int rank, np;
    MPI_Comm_rank(surf->communicator(), &rank);
    MPI_Comm_size(surf->communicator(), &np);
    if(rank == 0) {
      for(int k=0; k<nsys; k++) {
	BoutReal kwave=k*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	dcomplex phase(cos(kwave*ts) , -sin(kwave*ts));
	a[k][0] *= phase;
      }
    }
    if(rank == np-1) {
      for(int k=0; k<nsys; k++) {
	BoutReal kwave=k*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	dcomplex phase(cos(kwave*ts) , sin(kwave*ts));
	c[k][mesh->ngy-5] *= phase;
      }
    }
    
    // Solve cyclic tridiagonal system for each k
    cr->setCoefs(nsys, a, b, c);
    cr->solve(nsys, rhsk, xk);
    
    // Put back into rhs array
    for(int k=0;k<nsys;k++) {
      for(int y=0;y<mesh->ngy-4;y++)
        rhs[y][k] = xk[k][y];
    }
    
    // Inverse Fourier transform 
    for(int y=0;y<mesh->ngy-4;y++)
      irfft(rhs[y], mesh->ngz-1, result[x][y+2]);
    
  }
  
  // Delete cyclic reduction object
  delete cr;

#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}

