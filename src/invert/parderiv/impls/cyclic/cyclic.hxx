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
 * - For CELL_YLOW implementation, boundary conditions are only 1st order accurate.
 *   Should be OK for preconditioners, which are allowed to be less accurate.
 *   Only 1st-order accurate one-sided derivative is possible in a tri-diagonal matrix and
 *   staggered mesh requires one-sided derivative as boundary condition.
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

#ifndef __INV_PAR_CR_H__
#define __INV_PAR_CR_H__

#include "invert_parderiv.hxx"
#include "dcomplex.hxx"
#include <globals.hxx>
#include "utils.hxx"

class InvertParCR : public InvertPar {
public:
  InvertParCR(Options *opt, CELL_LOC location = CELL_CENTRE,
              Mesh *mesh_in = bout::globals::mesh);

  using InvertPar::solve;
  const Field3D solve(const Field3D &f) override;

  using InvertPar::setCoefA;
  void setCoefA(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    A = f;
  }
  using InvertPar::setCoefB;
  void setCoefB(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    B = f;
  }
  using InvertPar::setCoefC;
  void setCoefC(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    C = f;
  }
  using InvertPar::setCoefD;
  void setCoefD(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    D = f;
  }
  using InvertPar::setCoefE;
  void setCoefE(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    E = f;
  }

private:
  Field2D A{0.0}, B{0.0}, C{0.0}, D{0.0}, E{0.0};
  Field2D sg; // Coefficient of DDY contribution to Grad2_par2
  
  int nsys;
};

namespace {
RegisterInvertPar<InvertParCR> registerinvertparcyclic{PARDERIVCYCLIC};
}

#endif // __INV_PAR_CR_H__
