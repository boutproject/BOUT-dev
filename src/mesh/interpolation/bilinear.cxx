/**************************************************************************
 * Copyright 2015 B.D.Dudson, P. Hill
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

#include "bout/mesh.hxx"
#include "globals.hxx"
#include "interpolation.hxx"

#include <string>
#include <vector>

Bilinear::Bilinear(int y_offset, Mesh *mesh)
  : Interpolation(y_offset, mesh),
    w0(localmesh), w1(localmesh), w2(localmesh), w3(localmesh) {

  // Index arrays contain guard cells in order to get subscripts right
  i_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);
  k_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);

  // Allocate Field3D members
  w0.allocate();
  w1.allocate();
  w2.allocate();
  w3.allocate();
}

void Bilinear::calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                           const std::string& region) {

  BOUT_FOR(i, delta_x.getRegion(region)) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    if (skip_mask(x, y, z)) continue;

    // The integer part of xt_prime, zt_prime are the indices of the cell
    // containing the field line end-point
    i_corner(x, y, z) = static_cast<int>(floor(delta_x(x, y, z)));
    k_corner(x, y, z) = static_cast<int>(floor(delta_z(x, y, z)));

    // t_x, t_z are the normalised coordinates \in [0,1) within the cell
    // calculated by taking the remainder of the floating point index
    BoutReal t_x = delta_x(x, y, z) - static_cast<BoutReal>(i_corner(x, y, z));
    BoutReal t_z = delta_z(x, y, z) - static_cast<BoutReal>(k_corner(x, y, z));
    BoutReal t_x1 = 1.0 - t_x;
    BoutReal t_z1 = 1.0 - t_z;

    // Check that t_x and t_z are in range
    if ((t_x < 0.0) || (t_x > 1.0)) {
      throw BoutException(
          "t_x={:e} out of range at ({:d},{:d},{:d}) (delta_x={:e}, i_corner={:d})", t_x, x, y,
          z, delta_x(x, y, z), i_corner(x, y, z));
    }

    if ((t_z < 0.0) || (t_z > 1.0)) {
      throw BoutException(
          "t_z={:e} out of range at ({:d},{:d},{:d}) (delta_z={:e}, k_corner={:d})", t_z, x, y,
          z, delta_z(x, y, z), k_corner(x, y, z));
    }

    w0(x,y,z) = t_x1 * t_z1;
    w1(x,y,z) = t_x  * t_z1;
    w2(x,y,z) = t_x1 * t_z;
    w3(x,y,z) = t_x  * t_z;

  }
}

void Bilinear::calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                           const BoutMask &mask, const std::string& region) {
  skip_mask = mask;
  calcWeights(delta_x, delta_z, region);
}

Field3D Bilinear::interpolate(const Field3D& f, const std::string& region) const {
  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  BOUT_FOR(i, f.getRegion(region)) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    if (skip_mask(x, y, z)) continue;

    int y_next = y + y_offset;
    // Due to lack of guard cells in z-direction, we need to ensure z-index
    // wraps around
    int ncz = localmesh->LocalNz;
    int z_mod = ((k_corner(x, y, z) % ncz) + ncz) % ncz;
    int z_mod_p1 = (z_mod + 1) % ncz;

    f_interp(x, y_next, z) = f(i_corner(x, y, z), y_next, z_mod) * w0(x, y, z) +
      f(i_corner(x, y, z) + 1, y_next, z_mod) * w1(x, y, z) +
      f(i_corner(x, y, z), y_next, z_mod_p1) * w2(x, y, z) +
      f(i_corner(x, y, z) + 1, y_next, z_mod_p1) * w3(x, y, z);
  }
  return f_interp;
}

Field3D Bilinear::interpolate(const Field3D& f, const Field3D &delta_x,
                              const Field3D &delta_z, const std::string& region) {
  calcWeights(delta_x, delta_z, region);
  return interpolate(f, region);
}

Field3D Bilinear::interpolate(const Field3D& f, const Field3D &delta_x,
                              const Field3D &delta_z, const BoutMask &mask,
                              const std::string& region) {
  calcWeights(delta_x, delta_z, mask, region);
  return interpolate(f, region);
}
