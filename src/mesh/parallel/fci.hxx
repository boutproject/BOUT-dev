
/**************************************************************************
 * Flux-coordinate Independent parallel derivatives
 *
 **************************************************************************
 * Copyright 2014 B.D.Dudson, P. Hill
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

#include <bout/paralleltransform.hxx>

// Field line map - contains the coefficients for interpolation
class FCIMap {
  // Private constructor - must be initialised with mesh
  FCIMap();
public:
  // dir MUST be either +1 or -1
  FCIMap(Mesh& mesh, int dir);

  int*** i_corner;				// x-index of bottom-left grid point
  int*** k_corner;				// z-index of bottom-left grid point

  // Basis functions for cubic Hermite spline interpolation
  //	see http://en.wikipedia.org/wiki/Cubic_Hermite_spline
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

class FCITransform : public ParallelTransform {
public:
  FCITransform(Mesh *m);

  void calcYUpDown(Field3D &f);
  
  const Field3D toFieldAligned(const Field3D &f) {
    throw BoutException("FCI method cannot transform into field aligned grid");
  }

  const Field3D fromFieldAligned(const Field3D &f) {
    throw BoutException("FCI method cannot transform into field aligned grid");
  }
private:
  FCITransform();

  // The FCI object is tied to the mesh it was created on, so the mesh must
  // not change
  Mesh &mesh;

  // The maps hold the grid indices of the field line end-points and
  // associated interpolation coefficients
  const FCIMap forward_map;
  const FCIMap backward_map;
  
  // Interpolate field in direction DIR
  void interpolate(Field3D &f, Field3D &f_next, const FCIMap &fcimap, int dir);
};

