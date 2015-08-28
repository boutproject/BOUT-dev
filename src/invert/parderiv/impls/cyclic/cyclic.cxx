/************************************************************************
 * Inversion of parallel derivatives
 *
 * Inverts a matrix of the form
 *
 * A + B * Grad2_par2 + C*D2DYDZ + + D*D2DZ2 + E*DDY
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
#include <msg_stack.hxx>
#include <bout/constants.hxx>

#include <bout/surfaceiter.hxx>

#include <cmath>

InvertParCR::InvertParCR(Options *opt) : InvertPar(opt), A(1.0), B(0.0), C(0.0), D(0.0), E(0.0) {
  // Number of k equations to solve for each x location
  nsys = 1 + (mesh->ngz-1)/2;

  rhs = cmatrix(mesh->ngy, nsys);

  // Find out if we are on a boundary
  int size = mesh->ngy-4;
  SurfaceIter surf(mesh);
  for(surf.first(); !surf.isDone(); surf.next()) {
    BoutReal ts;
    int n = mesh->ngy-4;
    if(!surf.closed(ts)) {
      // Open field line
      if(surf.firstY())
        n += 2;
      if(surf.lastY())
        n += 2;

      if(n > size)
        size = n; // Maximum size
    }
  }

  rhsk = cmatrix(nsys, size);
  xk = cmatrix(nsys, size);
  a = cmatrix(nsys, size);
  b = cmatrix(nsys, size);
  c = cmatrix(nsys, size);
}

InvertParCR::~InvertParCR() {
  free_cmatrix(rhs);

  free_cmatrix(rhsk);
  free_cmatrix(xk);
  free_cmatrix(a);
  free_cmatrix(b);
  free_cmatrix(c);
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
  SurfaceIter surf(mesh);
  for(surf.first(); !surf.isDone(); surf.next()) {
    int x = surf.xpos;

    // Test if open or closed field-lines
    BoutReal ts;
    bool closed = surf.closed(ts);

    // Number of rows
    int y0 = 0, size = mesh->ngy-4; // If no boundaries
    if(!closed) {
      if(surf.firstY()) {
        y0 += 2;
        size += 2;
      }
      if(surf.lastY())
        size += 2;
    }

    // Setup CyclicReduce object
    cr->setup(surf.communicator(), size);
    cr->setPeriodic(closed);

    // Take Fourier transform
    for(int y=0;y<mesh->ngy-4;y++)
      rfft(f[x][y+2], mesh->ngz-1, rhs[y+y0]);

    // Set up tridiagonal system
    for(int k=0; k<nsys; k++) {
      BoutReal kwave=k*2.0*PI/mesh->zlength; // wave number is 1/[rad]
      for(int y=0;y<mesh->ngy-4;y++) {

        BoutReal acoef = A(x, y+2);                     // Constant
        BoutReal bcoef = B(x, y+2) / mesh->g_22(x,y+2); // d2dy2
        BoutReal ccoef = C(x, y+2);                     // d2dydz
        BoutReal dcoef = D(x, y+2);                     // d2dz2
        BoutReal ecoef = E(x, y+2);                     // ddy

        bcoef /= SQ(mesh->dy(x, y+2));
        ccoef /= mesh->dy(x,y+2)*mesh->dz;
        dcoef /= SQ(mesh->dz);
        ecoef /= mesh->dy(x,y+2);

        //           const     d2dy2        d2dydz             d2dz2           ddy
        //           -----     -----        ------             -----           ---
        a[k][y+y0] =            bcoef - 0.5*Im*kwave*ccoef                  -0.5*ecoef;
        b[k][y+y0] = acoef - 2.*bcoef                     - SQ(kwave)*dcoef;
        c[k][y+y0] =            bcoef + 0.5*Im*kwave*ccoef                  +0.5*ecoef;

        rhsk[k][y+y0] = rhs[y+y0][k]; // Transpose
      }
    }

    if(closed) {
      // Twist-shift
      int rank, np;
      MPI_Comm_rank(surf.communicator(), &rank);
      MPI_Comm_size(surf.communicator(), &np);
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
    }else {
      // Open surface, so may have boundaries
      if(surf.firstY()) {
        for(int k=0; k<nsys; k++) {
          for(int y=0;y<2;y++) {
            a[k][y] =  0.;
            b[k][y] =  1.;
            c[k][y] = -1.;

            rhsk[k][y] = 0.;
          }
        }
      }
      if(surf.lastY()) {
        for(int k=0; k<nsys; k++) {
          for(int y=size-2;y<size;y++) {
            a[k][y] = -1.;
            b[k][y] =  1.;
            c[k][y] =  0.;

            rhsk[k][y] = 0.;
          }
        }
      }
    }

    // Solve cyclic tridiagonal system for each k
    cr->setCoefs(nsys, a, b, c);
    cr->solve(nsys, rhsk, xk);

    // Put back into rhs array
    for(int k=0;k<nsys;k++) {
      for(int y=0;y<size;y++)
        rhs[y][k] = xk[k][y];
    }

    // Inverse Fourier transform
    for(int y=0;y<size;y++)
      irfft(rhs[y], mesh->ngz-1, result[x][y+2-y0]);

  }

  // Delete cyclic reduction object
  delete cr;

#ifdef CHECK
  msg_stack.pop();
#endif

  return result;
}

