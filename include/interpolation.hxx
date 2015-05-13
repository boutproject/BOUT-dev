/**************************************************************************
 * Functions to interpolate between cell locations (e.g. lower Y and centred)
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
 **************************************************************************/

#ifndef __INTERP_H__
#define __INTERP_H__

#include "field3d.hxx"
#include "bout_types.hxx"

/// Interpolate to a give cell location
const Field3D interp_to(const Field3D &var, CELL_LOC loc);
const Field2D interp_to(const Field2D &var, CELL_LOC loc);

/// Print out the cell location (for debugging)
void printLocation(const Field3D &var);

const char* strLocation(CELL_LOC loc);


/// Interpolate a field onto a perturbed set of points
const Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z);

const Field3D interpolate(const Field2D &f, const Field3D &delta_x, const Field3D &delta_z);
const Field3D interpolate(const Field2D &f, const Field3D &delta_x);

////////////////////////////////////////

class Interpolation {
public:
  virtual void calcWeights(const Field3D &delta_x, const Field3D &delta_z) = 0;

  virtual const Field3D interpolate(const Field3D& f) const = 0;
  virtual const Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z) = 0;
};

class HermiteSpline : public Interpolation {
public:
  HermiteSpline(int y_offset=0);

  void calcWeights(const Field3D &delta_x, const Field3D &delta_z);

  // Use precalculated weights
  const Field3D interpolate(const Field3D& f) const;
  // Calculate weights and interpolate
  const Field3D interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z);

  int*** i_corner;      // x-index of bottom-left grid point
  int*** k_corner;      // z-index of bottom-left grid point

private:
  // Interpolate using the field at (x,y+y_offset,z), rather than (x,y,z)
  int y_offset;

  // 3D vector of points to skip (true -> skip this point)
  BoutMask skip_mask;

  // Basis functions for cubic Hermite spline interpolation
  //    see http://en.wikipedia.org/wiki/Cubic_Hermite_spline
  // The h00 and h01 basis functions are applied to the function itself
  // and the h10 and h11 basis functions are applied to its derivative
  // along the interpolation direction.

  Field3D h00_x;
  Field3D h01_x;
  Field3D h10_x;
  Field3D h11_x;
  Field3D h00_z;
  Field3D h01_z;
  Field3D h10_z;
  Field3D h11_z;
};

#endif // __INTERP_H__
