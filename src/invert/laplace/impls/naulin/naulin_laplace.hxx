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
  LaplaceNaulin(Options *opt = NULL, const CELL_LOC loc = CELL_CENTRE, Mesh *mesh_in = nullptr);
  ~LaplaceNaulin();
  
  using Laplacian::setCoefA;
  using Laplacian::setCoefC;
  using Laplacian::setCoefC1;
  using Laplacian::setCoefC2;
  using Laplacian::setCoefD;
  using Laplacian::setCoefEx;
  using Laplacian::setCoefEz;

  // ACoef is not implemented because the delp2solver that we use can probably
  // only handle a Field2D for Acoef, but we would have to pass Acoef/Dcoef,
  // where we allow Dcoef to be a Field3D
  void setCoefA(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Acoef = val;
  }
  void setCoefA(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Acoef = val;
  }
  void setCoefC(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    setCoefC1(val);
    setCoefC2(val);
  }
  void setCoefC(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    setCoefC1(val);
    setCoefC2(val);
  }
  void setCoefC1(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1coef = val;
  }
  void setCoefC1(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1coef = val;
  }
  void setCoefC2(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C2coef = val;
  }
  void setCoefC2(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C2coef = val;
  }
  void setCoefD(const Field3D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Dcoef = val;
  }
  void setCoefD(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Dcoef = val;
  }
  void setCoefEx(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceNaulin does not have Ex coefficient");
  }
  void setCoefEz(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceNaulin does not have Ez coefficient");
  }

  bool uses3DCoefs() const override { return true; }

  using Laplacian::solve;

  FieldPerp solve(const FieldPerp &b) override {return solve(b,b);}
  FieldPerp solve(const FieldPerp &UNUSED(b),
                        const FieldPerp &UNUSED(x0)) override {
    throw BoutException(
        "LaplaceNaulin has no solve(FieldPerp), must call solve(Field3D)");
  }
  Field3D solve(const Field3D &b, const Field3D &x0) override;
  Field3D solve(const Field3D &b) override {
    return solve(b, zeroFrom(b));
  }

  // Override flag-setting methods to set delp2solver's flags as well
  void setGlobalFlags(int f) override { Laplacian::setGlobalFlags(f); delp2solver->setGlobalFlags(f); }
  void setInnerBoundaryFlags(int f) override { Laplacian::setInnerBoundaryFlags(f); delp2solver->setInnerBoundaryFlags(f); }
  void setOuterBoundaryFlags(int f) override { Laplacian::setOuterBoundaryFlags(f); delp2solver->setOuterBoundaryFlags(f); }

  BoutReal getMeanIterations() const { return naulinsolver_mean_its; }
  void resetMeanIterations() { naulinsolver_mean_its = 0; }
private:
  LaplaceNaulin(const LaplaceNaulin&);
  LaplaceNaulin& operator=(const LaplaceNaulin&);
  Field3D Acoef, C1coef, C2coef, Dcoef;

  /// Laplacian solver used to solve the equation with constant-in-z coefficients
  Laplacian* delp2solver;

  /// Solver tolerances
  BoutReal rtol, atol;

  /// Maximum number of iterations
  int maxits;

  /// Initial choice for under-relaxation factor, should be greater than 0 and
  /// less than or equal to 1. Value of 1 means no underrelaxation
  BoutReal initial_underrelax_factor{1.};

  /// Mean number of iterations taken by the solver
  BoutReal naulinsolver_mean_its;

  /// Mean number of times the underrelaxation factor is reduced
  BoutReal naulinsolver_mean_underrelax_counts{0.};

  /// Counter for the number of times the solver has been called
  int ncalls;

  /// Copy the boundary guard cells from the input 'initial guess' x0 into x.
  /// These may be used to set non-zero-value boundary conditions
  void copy_x_boundaries(Field3D &x, const Field3D &x0, Mesh *mesh);
};

#endif // __LAP_NAULIN_H__
