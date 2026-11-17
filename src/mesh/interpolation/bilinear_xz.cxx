/**************************************************************************
 * Copyright 2015 - 2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/mesh.hxx"
#include "bout/region.hxx"

#include <cmath>
#include <string>

XZBilinear::XZBilinear(int y_offset, Mesh* mesh)
    : XZInterpolation(y_offset, mesh), w0(localmesh), w1(localmesh), w2(localmesh),
      w3(localmesh) {

  if (localmesh->getNXPE() > 1) {
    throw BoutException("XZBilinear interpolation does not support MPI splitting in X");
  }

  // Index arrays contain guard cells in order to get subscripts right
  i_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);
  k_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);

  // Allocate Field3D members
  w0.allocate();
  w1.allocate();
  w2.allocate();
  w3.allocate();
}

void XZBilinear::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                             const std::string& region) {

  const auto curregion{getRegion(region)};
  BOUT_FOR(i, curregion) {
    // The integer part of xt_prime, zt_prime are the indices of the cell
    // containing the field line end-point
    i_corner[i] = static_cast<int>(floor(delta_x[i]));
    k_corner[i] = static_cast<int>(floor(delta_z[i]));

    // t_x, t_z are the normalised coordinates \in [0,1) within the cell
    // calculated by taking the remainder of the floating point index
    const BoutReal t_x = delta_x[i] - static_cast<BoutReal>(i_corner[i]);
    const BoutReal t_z = delta_z[i] - static_cast<BoutReal>(k_corner[i]);
    const BoutReal t_x1 = 1.0 - t_x;
    const BoutReal t_z1 = 1.0 - t_z;

    // Check that t_x and t_z are in range
    if ((t_x < 0.0) || (t_x > 1.0)) {
      throw BoutException(
          "t_x={:e} out of range at ({:d},{:d},{:d}) (delta_x={:e}, i_corner={:d})", t_x,
          i.x(), i.y(), i.z(), delta_x[i], i_corner[i]);
    }

    if ((t_z < 0.0) || (t_z > 1.0)) {
      throw BoutException(
          "t_z={:e} out of range at ({:d},{:d},{:d}) (delta_z={:e}, k_corner={:d})", t_z,
          i.x(), i.y(), i.z(), delta_z[i], k_corner[i]);
    }

    w0[i] = t_x1 * t_z1;
    w1[i] = t_x * t_z1;
    w2[i] = t_x1 * t_z;
    w3[i] = t_x * t_z;
  }
}

void XZBilinear::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                             const BoutMask& mask, const std::string& region) {
  setMask(mask);
  calcWeights(delta_x, delta_z, region);
}

Field3D XZBilinear::interpolate(const Field3D& f, const std::string& region) const {
  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  const auto curregion{getRegion(region)};
  const int ncz = localmesh->LocalNz;
  BOUT_FOR(i, curregion) {
    const auto iyp = i.yp(y_offset);
    // The stored corner coordinates define the structured source mapping.
    // Wrap Z over the full storage line, including any Z guards.
    const int z_mod = ((k_corner[i] % ncz) + ncz) % ncz;
    const int x_corner = i_corner[i];
    ASSERT3(x_corner >= 0 && x_corner < localmesh->LocalNx - 1);
    const auto source = f.indexAt(x_corner, iyp.y(), z_mod);
    const auto source_xp = source.xp();
    const auto source_zp = source.zp();
    const auto source_xpzp = source_xp.zp();

    f_interp[iyp] = (f[source] * w0[i]) + (f[source_xp] * w1[i]) + (f[source_zp] * w2[i])
                    + (f[source_xpzp] * w3[i]);
  }
  return f_interp;
}

Field3D XZBilinear::interpolate(const Field3D& f, const Field3D& delta_x,
                                const Field3D& delta_z, const std::string& region) {
  calcWeights(delta_x, delta_z, region);
  return interpolate(f, region);
}

Field3D XZBilinear::interpolate(const Field3D& f, const Field3D& delta_x,
                                const Field3D& delta_z, const BoutMask& mask,
                                const std::string& region) {
  calcWeights(delta_x, delta_z, mask, region);
  return interpolate(f, region);
}
