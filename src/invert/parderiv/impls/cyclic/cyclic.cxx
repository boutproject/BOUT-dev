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

InvertParCR::InvertParCR(Options *opt, Mesh *mesh_in)
  : InvertPar(opt, mesh_in), A(1.0), B(0.0), C(0.0), D(0.0), E(0.0) {
#ifdef COORDINATES_USE_3D
  throw BoutException("Parallel cyclic solver does not support 3D metric yet.");
#endif

  // Number of k equations to solve for each x location
  nsys = 1 + (localmesh->LocalNz)/2; 
}

const Field3D InvertParCR::solve(const Field3D &f) {
#ifndef COORDINATES_USE_3D
  TRACE("InvertParCR::solve(Field3D)");
  ASSERT1(localmesh == f.getMesh());

  Field3D result = emptyFrom(f).setDirectionY(YDirectionType::Aligned);
  
  Coordinates *coord = f.getCoordinates();

  Field3D alignedField = toFieldAligned(f, "RGN_NOX");

  // Create cyclic reduction object
  auto cr = bout::utils::make_unique<CyclicReduce<dcomplex>>();

  // Find out if we are on a boundary
  int size = localmesh->LocalNy - 2 * localmesh->ystart;
  SurfaceIter surf(localmesh);
  for(surf.first(); !surf.isDone(); surf.next()) {
    int n = localmesh->LocalNy - 2 * localmesh->ystart;
    if (!surf.closed()) {
      // Open field line
      if (surf.firstY())
        n += localmesh->ystart;
      if (surf.lastY())
        n += localmesh->ystart;

      if (n > size)
        size = n; // Maximum size
    }
  }

  auto rhs = Matrix<dcomplex>(localmesh->LocalNy, nsys);
  auto rhsk = Matrix<dcomplex>(nsys, size);
  auto xk = Matrix<dcomplex>(nsys, size);
  auto a = Matrix<dcomplex>(nsys, size);
  auto b = Matrix<dcomplex>(nsys, size);
  auto c = Matrix<dcomplex>(nsys, size);

  // Loop over flux-surfaces
  for (surf.first(); !surf.isDone(); surf.next()) {
    int x = surf.xpos;
    
    // Test if open or closed field-lines
    BoutReal ts;
    bool closed = surf.closed(ts);

    // Number of rows
    int y0 = 0;
    size = localmesh->LocalNy - 2 * localmesh->ystart; // If no boundaries
    if(!closed) {
      if(surf.firstY()) {
        y0 += localmesh->ystart;
        size += localmesh->ystart;
      }
      if(surf.lastY())
        size += localmesh->ystart;
    }
      
    // Setup CyclicReduce object
    cr->setup(surf.communicator(), size);
    cr->setPeriodic(closed);
    
    // Take Fourier transform
    for (int y = 0; y < localmesh->LocalNy - 2 * localmesh->ystart; y++)
      rfft(alignedField(x, y + localmesh->ystart), localmesh->LocalNz, &rhs(y + y0, 0));

    // Set up tridiagonal system
    for(int k=0; k<nsys; k++) {
      BoutReal kwave=k*2.0*PI/coord->zlength(); // wave number is 1/[rad]
      for (int y = 0; y < localmesh->LocalNy - 2 * localmesh->ystart; y++) {

        BoutReal acoef = A(x, y + localmesh->ystart); // Constant
        BoutReal bcoef =
            B(x, y + localmesh->ystart) / coord->g_22(x, y + localmesh->ystart); // d2dy2
        BoutReal ccoef = C(x, y + localmesh->ystart);                            // d2dydz
        BoutReal dcoef = D(x, y + localmesh->ystart);                            // d2dz2
        BoutReal ecoef = E(x, y + localmesh->ystart);                            // ddy

        bcoef /= SQ(coord->dy(x, y + localmesh->ystart));
        ccoef /= coord->dy(x, y + localmesh->ystart) * coord->dz;
        dcoef /= SQ(coord->dz);
        ecoef /= coord->dy(x, y + localmesh->ystart);

        //           const       d2dy2        d2dydz              d2dz2           ddy
        //           -----       -----        ------              -----           ---
        a(k, y + y0) =           bcoef - 0.5 * Im * kwave * ccoef          - 0.5 * ecoef;
        b(k, y + y0) = acoef - 2. * bcoef           - SQ(kwave) * dcoef;
        c(k, y + y0) =           bcoef + 0.5 * Im * kwave * ccoef          + 0.5 * ecoef;

        rhsk(k, y + y0) = rhs(y + y0, k); // Transpose
      }
    }

    if(closed) {
      // Twist-shift
      int rank, np;
      MPI_Comm_rank(surf.communicator(), &rank);
      MPI_Comm_size(surf.communicator(), &np);
      if(rank == 0) {
        for(int k=0; k<nsys; k++) {
          BoutReal kwave=k*2.0*PI/coord->zlength(); // wave number is 1/[rad]
          dcomplex phase(cos(kwave*ts) , -sin(kwave*ts));
          a(k, 0) *= phase;
        }
      }
      if(rank == np-1) {
        for(int k=0; k<nsys; k++) {
          BoutReal kwave=k*2.0*PI/coord->zlength(); // wave number is 1/[rad]
          dcomplex phase(cos(kwave*ts) , sin(kwave*ts));
          c(k, localmesh->LocalNy - 2 * localmesh->ystart - 1) *= phase;
        }
      }
    }else {
      // Open surface, so may have boundaries
      if(surf.firstY()) {
        for(int k=0; k<nsys; k++) {
          for (int y = 0; y < localmesh->ystart; y++) {
            a(k, y) = 0.;
            b(k, y) = 1.;
            c(k, y) = -1.;

            rhsk(k, y) = 0.;
          }
        }
      }
      if(surf.lastY()) {
        for(int k=0; k<nsys; k++) {
          for (int y = size - localmesh->ystart; y < size; y++) {
            a(k, y) = -1.;
            b(k, y) = 1.;
            c(k, y) = 0.;

            rhsk(k, y) = 0.;
          }
        }
      }
    }
    
    // Solve cyclic tridiagonal system for each k
    cr->setCoefs(a, b, c);
    cr->solve(rhsk, xk);

    // Put back into rhs array
    for(int k=0;k<nsys;k++) {
      for(int y=0;y<size;y++)
        rhs(y, k) = xk(k, y);
    }
    
    // Inverse Fourier transform 
    for(int y=0;y<size;y++)
      irfft(&rhs(y, 0), localmesh->LocalNz, result(x, y + localmesh->ystart - y0));
  }

  return fromFieldAligned(result, "RGN_NOBNDRY");
#else
  return Field3D{};
#endif
};


