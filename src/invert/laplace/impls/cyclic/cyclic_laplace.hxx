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

/// 
/*!
 * 
 */
class LaplaceCyclic : public Laplacian {
public:
  LaplaceCyclic(Options *opt = NULL);
  ~LaplaceCyclic();
  
  void setCoefA(const Field2D &val) { A = val; }
  void setCoefC(const Field2D &val) { C = val; }
  void setCoefD(const Field2D &val) { D = val; }
  void setCoefEx(const Field2D &val) { bout_error("LaplaceCyclic does not have Ex coefficient"); }
  void setCoefEz(const Field2D &val) { bout_error("LaplaceCyclic does not have Ez coefficient"); }
  
  const FieldPerp solve(const FieldPerp &b) {return solve(b,b);}
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0);
private:
  Field2D A, C, D;
  
  int nmode;  // Number of modes being solved
  int xs, xe; // Start and end X indices
  dcomplex **a, **b, **c, **bcmplx, **xcmplx;
  dcomplex *k1d;
  
  bool dst;
  
  CyclicReduce<dcomplex> *cr; ///< Tridiagonal solver
};

#endif // __SPT_H__
