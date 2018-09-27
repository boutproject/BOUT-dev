/**************************************************************************
 * Iterative solver to handle non-constant-in-z coefficients
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
  LaplaceNaulin(Options *opt = NULL, const CELL_LOC loc = CELL_CENTRE);
  ~LaplaceNaulin();
  
  // ACoef is not implemented because the delp2solver that we use can probably
  // only handle a Field2D for Acoef, but we would have to pass Acoef/Dcoef,
  // where we allow Dcoef to be a Field3D
  void setCoefA(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    Acoef = val;
  }
  void setCoefA(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    Acoef = val;
  }
  void setCoefC(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    setCoefC1(val);
    setCoefC2(val);
  }
  void setCoefC(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    setCoefC1(val);
    setCoefC2(val);
  }
  void setCoefC1(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    C1coef = val;
  }
  void setCoefC1(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    C1coef = val;
  }
  void setCoefC2(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    C2coef = val;
  }
  void setCoefC2(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    C2coef = val;
  }
  void setCoefD(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    Dcoef = val;
  }
  void setCoefD(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    Dcoef = val;
  }
  void setCoefEx(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceNaulin does not have Ex coefficient");
  }
  void setCoefEz(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceNaulin does not have Ez coefficient");
  }

  const FieldPerp solve(const FieldPerp &b) override {return solve(b,b);}
  const FieldPerp solve(const FieldPerp &UNUSED(b),
                        const FieldPerp &UNUSED(x0)) override {
    throw BoutException(
        "LaplaceNaulin has no solve(FieldPerp), must call solve(Field3D)");
  }
  const Field3D solve(const Field3D &b, const Field3D &x0) override;
  const Field3D solve(const Field3D &b) override {
    Field3D x0(b.getMesh());
    x0 = 0.;
    x0.setLocation(b.getLocation());
    return solve(b, Field3D(0.));
  }

  // Override flag-setting methods to set delp2solver's flags as well
  void setGlobalFlags(int f) override { Laplacian::setGlobalFlags(f); delp2solver->setGlobalFlags(f); }
  void setInnerBoundaryFlags(int f) override { Laplacian::setInnerBoundaryFlags(f); delp2solver->setInnerBoundaryFlags(f); }
  void setOuterBoundaryFlags(int f) override { Laplacian::setOuterBoundaryFlags(f); delp2solver->setOuterBoundaryFlags(f); }

  BoutReal getMeanIterations() const { return naulinsolver_mean_its; }
private:
  LaplaceNaulin(const LaplaceNaulin&);
  LaplaceNaulin& operator=(const LaplaceNaulin&);
  Field3D Acoef, C1coef, C2coef, Dcoef;
  Laplacian* delp2solver;
  BoutReal rtol, atol;
  int maxits;
  BoutReal naulinsolver_mean_its;
  int ncalls;

  void copy_x_boundaries(Field3D &x, const Field3D &x0, Mesh *mesh);
};

#endif // __LAP_NAULIN_H__
