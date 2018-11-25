/**************************************************************************
 * Calculate the coefficients for a 1d cubic spline interpolation
 * by using the 4th order derivatives at the points
 *
 **************************************************************************
 * Copyright 2012 J.T.Omotani
 *
 * Contact: John Omotani, john.omotani@york.ac.uk
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

#ifndef __CUBICSPLINE_H__
#define __CUBICSPLINE_H__

#include <bout/constants.hxx>
#include <bout.hxx>
#include <bout_types.hxx>
#include <derivs.hxx>
#include <cmath>

class CubicSpline {
public:
  CubicSpline();
  ~CubicSpline();
  void initialise(const char &direction, const bool pass_include_boundary_guard_cells=true, const bool pass_staggered=false);
  void initialise(const CubicSpline &clone_from);
  void calculate(const Field3D &pass_input);
  BoutReal * coefficients(bindex* position);
  BoutReal * coefficients(const int &jx, const int &jy, const int &jz);
  
private:
  char direction;
  bool include_boundary_guard_cells;
  bool upper_boundary_offset;
  Field3D input;
  Field3D grad_input;
  BoutReal * coeffs;
};

#endif // __CUBICSPLINE_H__
