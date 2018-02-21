/**************************************************************************
 * Calculate the coefficients for a 1d cubic spline interpolation (in index space)
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

#include <bout/globals.hxx>
#include <bout/bout_types.hxx>

#include "cubic_spline_local.hxx"

CubicSpline::CubicSpline() {
}

CubicSpline::~CubicSpline() {
  delete [] coeffs;
}

void CubicSpline::initialise(const char &inputdirection, const bool pass_include_boundary_guard_cells, const bool pass_upper_boundary_offset) {
  direction = inputdirection;
  if (direction == 'y'){
    
  }
  else throw BoutException("CubicSpline: direction is not valid");
  include_boundary_guard_cells = pass_include_boundary_guard_cells;
  upper_boundary_offset=pass_upper_boundary_offset;
  coeffs = new BoutReal[4];
  for (int i=0; i<4; i++) coeffs[i] = 1./0.;
}

void CubicSpline::initialise(const CubicSpline &clone_from) {
  direction = clone_from.direction;
  include_boundary_guard_cells = clone_from.include_boundary_guard_cells;
  upper_boundary_offset = clone_from.upper_boundary_offset;
  coeffs = new BoutReal[4];
  for (int i=0; i<4; i++) coeffs[i] = NAN;
}

void CubicSpline::calculate(const Field3D &pass_input) {
  if (direction=='y') {
    input = pass_input;
    if (!include_boundary_guard_cells) {
      for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
	for (int jz=0; jz<mesh->LocalNz; jz++)
	  for (int jy=mesh->ystart-1; jy>=0; jy--)
	    input(rlow.ind,jy,jz) = 0.;
      for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
	for (int jz=0; jz<mesh->LocalNz; jz++)
// 	  for (int jy=mesh->yend+1-staggered; jy<mesh->LocalNy; jy++)
	  for (int jy=mesh->yend+1-upper_boundary_offset; jy<mesh->LocalNy; jy++)
	    input(rup.ind,jy,jz) = 0.;
	    
      grad_input = mesh->indexDDY(input,input.getLocation(), DIFF_DEFAULT); // derivative in index space
      
      // values of grad_input that depended on the guard cells are wrong, replace them with forward/backward derivatives
      for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
	for (int jz=0; jz<mesh->LocalNz; jz++)
	  for (int jy=mesh->ystart+1; jy>=0; jy--) {
	    grad_input(rlow.ind,jy,jz) = (-11.*input(rlow.ind,jy,jz) + 18.*input(rlow.ind,jy+1,jz) - 9.*input(rlow.ind,jy+2,jz) + 2.*input(rlow.ind,jy+3,jz)) / 6.;
	  }
      for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
	for (int jz=0; jz<mesh->LocalNz; jz++)
// 	  for (int jy=mesh->yend-1-staggered; jy<mesh->LocalNy; jy++) {
	  for (int jy=mesh->yend-1-upper_boundary_offset; jy<mesh->LocalNy; jy++) {
	    grad_input(rup.ind,jy,jz) = (11.*input(rup.ind,jy,jz) - 18.*input(rup.ind,jy-1,jz) + 9.*input(rup.ind,jy-2,jz) - 2.*input(rup.ind,jy-3,jz)) / 6.;
	  }
    }
    else {
      grad_input = mesh->indexDDY(input,input.getLocation(), DIFF_DEFAULT); // derivative in index space
      
      // May need values in the guard cells as well, calculate with forward/backward derivatives
      for (RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
	for (int jz=0; jz<mesh->LocalNz; jz++)
	  for (int jy=mesh->ystart-1; jy>=0; jy--) {
	    grad_input(rlow.ind,jy,jz) = (-11.*input(rlow.ind,jy,jz) + 18.*input(rlow.ind,jy+1,jz) - 9.*input(rlow.ind,jy+2,jz) + 2.*input(rlow.ind,jy+3,jz)) / 6.;
	  }
      for (RangeIterator rup = mesh->iterateBndryUpperY(); !rup.isDone(); rup++)
	for (int jz=0; jz<mesh->LocalNz; jz++)
// 	  for (int jy=mesh->yend+1-staggered; jy<mesh->LocalNy; jy++) {
	  for (int jy=mesh->yend+1-upper_boundary_offset; jy<mesh->LocalNy; jy++) {
	    grad_input(rup.ind,jy,jz) = (11.*input(rup.ind,jy,jz) - 18.*input(rup.ind,jy-1,jz) + 9.*input(rup.ind,jy-2,jz) - 2.*input(rup.ind,jy-3,jz)) / 6.;
	  }
    }
    mesh->communicate(grad_input);
  }
}

BoutReal* CubicSpline::coefficients (bindex* position) {
  coeffs[0] = input(position->jx,position->jy,position->jz);
  coeffs[1] = grad_input(position->jx,position->jy,position->jz);
  coeffs[2] = 3.*(input(position->jx,position->jyp,position->jz)-input(position->jx,position->jy,position->jz))
    - 2.*grad_input(position->jx,position->jy,position->jz) - grad_input(position->jx,position->jyp,position->jz);
  coeffs[3] = 2.*(input(position->jx,position->jy,position->jz)-input(position->jx,position->jyp,position->jz))
    + grad_input(position->jx,position->jy,position->jz) + grad_input(position->jx,position->jyp,position->jz);
  return coeffs;
}

BoutReal* CubicSpline::coefficients (const int &jx, const int &jy, const int &jz) {
  coeffs[0] = input(jx,jy,jz);
  coeffs[1] = grad_input(jx,jy,jz);
  coeffs[2] = 3.*(input(jx,jy+1,jz)-input(jx,jy,jz))
    - 2.*grad_input(jx,jy,jz) - grad_input(jx,jy+1,jz);
  coeffs[3] = 2.*(input(jx,jy,jz)-input(jx,jy+1,jz))
    + grad_input(jx,jy,jz) + grad_input(jx,jy+1,jz);
  return coeffs;
}
