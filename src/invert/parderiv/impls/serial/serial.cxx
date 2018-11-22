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

InvertParSerial::InvertParSerial(Options *opt, Mesh *mesh_in)
    : InvertPar(opt, mesh_in), A(1.0), B(0.0), C(0.0), D(0.0), E(0.0) {
}

const Field3D InvertParSerial::solve(const Field3D &f) {
  TRACE("InvertParSerial::solve(Field3D)");

  ASSERT1(localmesh == f.getMesh());

  Field3D result(localmesh);
  result.allocate();
  result.setLocation(f.getLocation());
  
  Coordinates *coord = f.getCoordinates();

  // The following is the number of nonguard points
  const int ny = localmesh->LocalNy - 2 * localmesh->ystart;

  auto rhs = Matrix<dcomplex>(localmesh->LocalNy, (localmesh->LocalNz) / 2 + 1);
  auto rhsk = Array<dcomplex>(ny);
  auto xk = Array<dcomplex>(ny);
  auto a = Array<dcomplex>(ny);
  auto b = Array<dcomplex>(ny);
  auto c = Array<dcomplex>(ny);

  Field3D alignedField = localmesh->toFieldAligned(f);

  // Loop over flux-surfaces
  SurfaceIter surf(localmesh);
  for(surf.first(); !surf.isDone(); surf.next()) {
    int x = surf.xpos;
    BoutReal ts; // Twist-shift angle
    if(!surf.closed(ts))
      throw BoutException("InvertParSerial doesn't handle open surfaces");
    
    // Take Fourier transform
    for (int y = localmesh->ystart; y <= localmesh->yend; y++)
      rfft(alignedField(x, y), localmesh->LocalNz, &rhs(y - localmesh->ystart, 0));

    // Solve cyclic tridiagonal system for each k
    const int nyq = (localmesh->LocalNz) / 2;
    const BoutReal waveFac = TWOPI / coord->zlength();
    for(int k=0;k<=nyq;k++) {
      // Copy component of rhs into 1D array
      for (int y = 0; y < ny; y++)
        rhsk[y] = rhs(y, k);

      BoutReal kwave = k * waveFac; // wave number is 1/[rad]

      // Set up tridiagonal system
      for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
        BoutReal acoef = A(x, y);                     // Constant
        BoutReal bcoef = B(x, y) / coord->g_22(x, y); // d2dy2
        BoutReal ccoef = C(x, y);                     // d2dydz
        BoutReal dcoef = D(x, y);                     // d2dz2
        BoutReal ecoef = E(x, y);                     // ddy

        bcoef /= SQ(coord->dy(x, y));
        ccoef /= coord->dy(x, y) * coord->dz;
        dcoef /= SQ(coord->dz);
        ecoef /= coord->dy(x, y);

        //                       const     d2dy2        d2dydz             d2dz2 ddy
        //                       -----     -----        ------             ----- ---
        a[y - localmesh->ystart] = bcoef - 0.5 * Im * kwave * ccoef - 0.5 * ecoef;
        b[y - localmesh->ystart] = acoef - 2. * bcoef - SQ(kwave) * dcoef;
        c[y - localmesh->ystart] = bcoef + 0.5 * Im * kwave * ccoef + 0.5 * ecoef;
      }

      // Modify coefficients across twist-shift -- accounts for pseudo-periodicity
      // across y-boundary to produce periodic value
      dcomplex phase(cos(kwave*ts) , -sin(kwave*ts));
      a[0] *= phase;
      c[ny - 1] /= phase;

      // Solve cyclic tridiagonal system
      cyclic_tridag(std::begin(a), std::begin(b), std::begin(c), std::begin(rhsk),
                    std::begin(xk), ny);

      // Put back into rhs array
      for (int y = 0; y < ny; y++)
        rhs(y, k) = xk[y];
    }
    
    // Inverse Fourier transform
    for (int y = localmesh->ystart; y <= localmesh->yend; y++)
      irfft(&rhs(y - localmesh->ystart, 0), localmesh->LocalNz, result(x, y));
  }

  return localmesh->fromFieldAligned(result);
}

