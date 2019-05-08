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

#ifndef __INV_PAR_CR_H__
#define __INV_PAR_CR_H__

#include "invert_parderiv.hxx"
#include "dcomplex.hxx"
#include <globals.hxx>
#include "utils.hxx"

class InvertParCR : public InvertPar {
public:
  InvertParCR(Options *opt, Mesh *mesh_in = bout::globals::mesh);
  ~InvertParCR();

  using InvertPar::solve;
  const Field3D solve(const Field3D &f) override;

  using InvertPar::setCoefA;
  void setCoefA(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    A = f;
  }
  using InvertPar::setCoefB;
  void setCoefB(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    B = f;
  }
  using InvertPar::setCoefC;
  void setCoefC(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    C = f;
  }
  using InvertPar::setCoefD;
  void setCoefD(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    D = f;
  }
  using InvertPar::setCoefE;
  void setCoefE(const Field2D &f) override {
    ASSERT1(localmesh == f.getMesh());
    E = f;
  }

private:
  Field2D A, B, C, D, E;
  
  int nsys;
};


#endif // __INV_PAR_CR_H__
