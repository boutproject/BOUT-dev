/**************************************************************************
 * Perpendicular Laplacian inversion in serial or parallel
 * 
 **************************************************************************
 * Copyright 2018 B.D.Dudson, M. Loiten, J. Omotani
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

class LaplaceNaulin;

#ifndef __LAP_NAULIN_H__
#define __LAP_NAULIN_H__

#include <invert_laplace.hxx>
#include <options.hxx>

/// Solves the 2D Laplacian equation
/*!
 * 
 */
class LaplaceNaulin : public Laplacian {
public:
  LaplaceNaulin(Options *opt = NULL);
  ~LaplaceNaulin();
  
  // ACoef is not implemented because the delp2solver that we use can probably
  // only handle a Field2D for Acoef, but we would have to pass Acoef/Dcoef,
  // where we allow Dcoef to be a Field3D
  using Laplacian::setCoefA;
  void setCoefA(const Field2D &val) override { throw BoutException("Acoef is not implemented"); }
  using Laplacian::setCoefC;
  void setCoefC(const Field2D &val) override { C1coef = val; C2coef = val; }
  using Laplacian::setCoefC1;
  void setCoefC1(const Field3D &val) override { C1coef = val; }
  void setCoefC1(const Field2D &val) override { C1coef = val; }
  using Laplacian::setCoefC2;
  void setCoefC2(const Field3D &val) override { C2coef = val; }
  void setCoefC2(const Field2D &val) override { C2coef = val; }
  using Laplacian::setCoefD;
  void setCoefD(const Field3D &val) override { Dcoef = val; }
  void setCoefD(const Field2D &val) override { Dcoef = val; }
  using Laplacian::setCoefEx;
  void setCoefEx(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceNaulin does not have Ex coefficient");
  }
  using Laplacian::setCoefEz;
  void setCoefEz(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceNaulin does not have Ez coefficient");
  }

  using Laplacian::solve;
  const FieldPerp solve(const FieldPerp &b) override {return solve(b,b);}
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0) override { throw BoutException("LaplaceNaulin has no solve(FieldPerp), must call solve(Field3D)"); }
  const Field3D solve(const Field3D &b, const Field3D&x0) override;
  const Field3D solve(const Field3D &b) override { return solve(b, Field3D(0.)); }
private:
  LaplaceNaulin(const LaplaceNaulin&);
  LaplaceNaulin& operator=(const LaplaceNaulin&);
  Field3D C1coef, C2coef, Dcoef;
  Laplacian* delp2solver;
  BoutReal rtol, atol;
  int maxits;

  void copy_x_boundaries(Field3D &x, const Field3D &x0, Mesh *mesh);
};

#endif // __LAP_NAULIN_H__
