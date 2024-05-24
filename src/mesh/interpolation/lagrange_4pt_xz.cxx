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

#include "bout/bout_types.hxx"
#include "bout/globals.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/mesh.hxx"

#include <array>

namespace {
// 4-point Lagrangian interpolation
// offset must be between 0 and 1
BoutReal lagrange_4pt(const BoutReal v2m, const BoutReal v1m, const BoutReal v1p,
                      const BoutReal v2p, const BoutReal offset) {
  return -offset * (offset - 1.0) * (offset - 2.0) * v2m / 6.0
         + 0.5 * (offset * offset - 1.0) * (offset - 2.0) * v1m
         - 0.5 * offset * (offset + 1.0) * (offset - 2.0) * v1p
         + offset * (offset * offset - 1.0) * v2p / 6.0;
}
// Convenience helper
BoutReal lagrange_4pt(const std::array<BoutReal, 4>& x, const BoutReal offset) {
  return lagrange_4pt(x[0], x[1], x[2], x[3], offset);
}

} // namespace

XZLagrange4pt::XZLagrange4pt(int y_offset, Mesh* mesh)
    : XZInterpolation(y_offset, mesh), t_x(localmesh), t_z(localmesh) {

  // Index arrays contain guard cells in order to get subscripts right
  i_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);
  k_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);

  t_x.allocate();
  t_z.allocate();
}

void XZLagrange4pt::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                                const std::string& region) {
  const auto curregion{getRegion(region)};
  BOUT_FOR(i, curregion) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    // The integer part of xt_prime, zt_prime are the indices of the cell
    // containing the field line end-point
    i_corner(x, y, z) = static_cast<int>(floor(delta_x(x, y, z)));
    k_corner(x, y, z) = static_cast<int>(floor(delta_z(x, y, z)));

    // t_x, t_z are the normalised coordinates \in [0,1) within the cell
    // calculated by taking the remainder of the floating point index
    t_x(x, y, z) = delta_x(x, y, z) - static_cast<BoutReal>(i_corner(x, y, z));
    t_z(x, y, z) = delta_z(x, y, z) - static_cast<BoutReal>(k_corner(x, y, z));

    // NOTE: A (small) hack to avoid one-sided differences
    if (i_corner(x, y, z) == localmesh->xend) {
      i_corner(x, y, z) -= 1;
      t_x(x, y, z) = 1.0;
    }

    // Check that t_x and t_z are in range
    if ((t_x(x, y, z) < 0.0) || (t_x(x, y, z) > 1.0)) {
      throw BoutException(
          "t_x={:e} out of range at ({:d},{:d},{:d}) (delta_x={:e}, i_corner={:d})",
          t_x(x, y, z), x, y, z, delta_x(x, y, z), i_corner(x, y, z));
    }
    if ((t_z(x, y, z) < 0.0) || (t_z(x, y, z) > 1.0)) {
      throw BoutException(
          "t_z={:e} out of range at ({:d},{:d},{:d}) (delta_z={:e}, k_corner={:d})",
          t_z(x, y, z), x, y, z, delta_z(x, y, z), k_corner(x, y, z));
    }
  }
}

Field3D XZLagrange4pt::interpolate(const Field3D& f, const std::string& region) const {

  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  const auto curregion{getRegion(region)};
  BOUT_FOR(i, curregion) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    const int jx = i_corner(x, y, z);
    const int jx2mnew = (jx == 0) ? 0 : (jx - 1);
    const int jxpnew = jx + 1;
    const int jx2pnew = (jx == (localmesh->LocalNx - 2)) ? jxpnew : (jxpnew + 1);

    const int ncz = localmesh->LocalNz;

    // Get the 4 Z points
    const int jz = ((k_corner(x, y, z) % ncz) + ncz) % ncz;

    const int jzpnew = (jz + 1) % ncz;
    const int jz2pnew = (jz + 2) % ncz;
    const int jz2mnew = (jz - 1 + ncz) % ncz;

    // Interpolate in Z first
    const int y_next = y + y_offset;

    const std::array<BoutReal, 4> xvals{
        lagrange_4pt(f(jx2mnew, y_next, jz2mnew), f(jx2mnew, y_next, jz),
                     f(jx2mnew, y_next, jzpnew), f(jx2mnew, y_next, jz2pnew),
                     t_z(x, y, z)),

        lagrange_4pt(f(jx, y_next, jz2mnew), f(jx, y_next, jz), f(jx, y_next, jzpnew),
                     f(jx, y_next, jz2pnew), t_z(x, y, z)),

        lagrange_4pt(f(jxpnew, y_next, jz2mnew), f(jxpnew, y_next, jz),
                     f(jxpnew, y_next, jzpnew), f(jxpnew, y_next, jz2pnew), t_z(x, y, z)),

        lagrange_4pt(f(jx2pnew, y_next, jz2mnew), f(jx2pnew, y_next, jz),
                     f(jx2pnew, y_next, jzpnew), f(jx2pnew, y_next, jz2pnew),
                     t_z(x, y, z)),
    };
    // Then in X
    f_interp(x, y_next, z) = lagrange_4pt(xvals, t_x(x, y, z));
  }
  return f_interp;
}
