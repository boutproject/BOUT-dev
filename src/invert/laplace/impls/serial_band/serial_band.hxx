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

#include <bout/invert_laplace.hxx>
#include <bout/dcomplex.hxx>
#include <bout/options.hxx>

class LaplaceSerialBand : public Laplacian {
public:
  LaplaceSerialBand(Options *opt = NULL);
  ~LaplaceSerialBand();
  
  using Laplacian::setCoefA;
  void setCoefA(const Field2D &val) override { Acoef = val; }
  using Laplacian::setCoefC;
  void setCoefC(const Field2D &val) override { Ccoef = val; }
  using Laplacian::setCoefD;
  void setCoefD(const Field2D &val) override { Dcoef = val; }
  using Laplacian::setCoefEx;
  void setCoefEx(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceSerialBand does not have Ex coefficient");
  }
  using Laplacian::setCoefEz;
  void setCoefEz(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceSerialBand does not have Ez coefficient");
  }

  using Laplacian::solve;
  const FieldPerp solve(const FieldPerp &b) override;
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0) override;
private:
  Field2D Acoef, Ccoef, Dcoef;
  
  dcomplex **bk, *bk1d;
  dcomplex **xk, *xk1d;

  dcomplex **A;
};

#endif // __SERIAL_BAND_H__
