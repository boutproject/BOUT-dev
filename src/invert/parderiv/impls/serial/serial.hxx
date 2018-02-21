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

#ifndef __INV_PAR_SERIAL_H__
#define __INV_PAR_SERIAL_H__

#include "bout/invert_parderiv.hxx"
#include "bout/dcomplex.hxx"

class InvertParSerial : public InvertPar {
public:
  InvertParSerial(Options* opt);
  ~InvertParSerial();

  using InvertPar::solve;
  const Field3D solve(const Field3D &f) override;

  using InvertPar::setCoefA;
  void setCoefA(const Field2D &f) override { A = f; }
  using InvertPar::setCoefB;
  void setCoefB(const Field2D &f) override { B = f; }
  using InvertPar::setCoefC;
  void setCoefC(const Field2D &f) override { C = f; }
  using InvertPar::setCoefD;
  void setCoefD(const Field2D &f) override { D = f; }
  using InvertPar::setCoefE;
  void setCoefE(const Field2D &f) override { E = f; }

private:
  Field2D A, B, C, D, E;
  
  dcomplex **rhs;
  dcomplex *rhsk;
  dcomplex *xk;
  dcomplex *a, *b, *c; // Matrix coefficients
};


#endif // __INV_PAR_SERIAL_H__
