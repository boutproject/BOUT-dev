/**************************************************************************
 * Perpendicular Laplacian inversion in serial or parallel
 * 
 * Simplified version of serial_tri and spt algorithms. Uses FFTs and
 * the cyclic reduction class to solve tridiagonal systems.
 *
 **************************************************************************
 * Copyright 2013 B.D.Dudson
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

class LaplaceCyclic;

#ifndef __LAP_CYCLIC_H__
#define __LAP_CYCLIC_H__

#include <invert_laplace.hxx>
#include <cyclic_reduction.hxx>
#include <dcomplex.hxx>
#include <options.hxx>

#include "utils.hxx"

#include "../../laplacefactory.hxx"

namespace {
RegisterLaplace<LaplaceCyclic> registerlaplacecycle(LAPLACE_CYCLIC);
}

/// Solves the 2D Laplacian equation using the CyclicReduce class
/*!
 * 
 */
class LaplaceCyclic : public Laplacian {
public:
  LaplaceCyclic(Options *opt = nullptr, const CELL_LOC loc = CELL_CENTRE, Mesh *mesh_in = nullptr);
  ~LaplaceCyclic();
  
  using Laplacian::setCoefA;
  void setCoefA(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Acoef = val;
  }
  using Laplacian::setCoefC;
  void setCoefC(const Field2D &val) override {
    setCoefC1(val);
    setCoefC2(val);
  }
  using Laplacian::setCoefC1;
  void setCoefC1(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C1coef = val;
  }
  using Laplacian::setCoefC2;
  void setCoefC2(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    C2coef = val;
  }
  using Laplacian::setCoefD;
  void setCoefD(const Field2D &val) override {
    ASSERT1(val.getLocation() == location);
    ASSERT1(localmesh == val.getMesh());
    Dcoef = val;
  }
  using Laplacian::setCoefEx;
  void setCoefEx(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceCyclic does not have Ex coefficient");
  }
  using Laplacian::setCoefEz;
  void setCoefEz(const Field2D &UNUSED(val)) override {
    throw BoutException("LaplaceCyclic does not have Ez coefficient");
  }

  using Laplacian::solve;
  FieldPerp solve(const FieldPerp &b) override {return solve(b,b);}
  FieldPerp solve(const FieldPerp &b, const FieldPerp &x0) override;

  Field3D solve(const Field3D &b) override {return solve(b,b);}
  Field3D solve(const Field3D &b, const Field3D &x0) override;
private:
  Field2D Acoef, C1coef, C2coef, Dcoef;
  
  int nmode;  // Number of modes being solved
  int xs, xe; // Start and end X indices
  Matrix<dcomplex> a, b, c, bcmplx, xcmplx;
  
  bool dst;
  
  CyclicReduce<dcomplex> *cr; ///< Tridiagonal solver
};

#endif // __SPT_H__
