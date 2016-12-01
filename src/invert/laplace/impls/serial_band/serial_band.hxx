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

class LaplaceSerialBand;

#ifndef __SERIAL_BAND_H__
#define __SERIAL_BAND_H__

#include <invert_laplace.hxx>
#include <dcomplex.hxx>
#include <options.hxx>

class LaplaceSerialBand : public Laplacian {
public:
  LaplaceSerialBand(Options *opt = NULL);
  ~LaplaceSerialBand();
  
  void setCoefA(const Field2D &val) { Acoef = val; }
  void setCoefC(const Field2D &val) { Ccoef = val; }
  void setCoefD(const Field2D &val) { Dcoef = val; }
  void setCoefEx(const Field2D &UNUSED(val)) { throw BoutException("LaplaceSPT does not have Ex coefficient"); }
  void setCoefEz(const Field2D &UNUSED(val)) { throw BoutException("LaplaceSPT does not have Ez coefficient"); }
  
  const FieldPerp solve(const FieldPerp &b);
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0);
private:
  Field2D Acoef, Ccoef, Dcoef;
  
  dcomplex **bk, *bk1d;
  dcomplex **xk, *xk1d;

  dcomplex **A;
};

#endif // __SERIAL_BAND_H__
