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
#include "interpolation_xz.hxx"

#include <vector>

XZLagrange4pt::XZLagrange4pt(int y_offset, Mesh *mesh)
    : XZInterpolation(y_offset, mesh), t_x(localmesh), t_z(localmesh) {

  // Index arrays contain guard cells in order to get subscripts right
  i_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);
  k_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);

  t_x.allocate();
  t_z.allocate();
}

void XZLagrange4pt::calcWeights(const Field3D &delta_x, const Field3D &delta_z) {

  for (int x = localmesh->xstart; x <= localmesh->xend; x++) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      for (int z = 0; z < localmesh->LocalNz; z++) {

        if (skip_mask(x, y, z))
          continue;

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
            "t_x={:e} out of range at ({:d},{:d},{:d}) (delta_x={:e}, i_corner={:d})", t_x(x, y, z), x, y,
              z, delta_x(x, y, z), i_corner(x, y, z));
        }
        if ((t_z(x, y, z) < 0.0) || (t_z(x, y, z) > 1.0)) {
          throw BoutException(
            "t_z={:e} out of range at ({:d},{:d},{:d}) (delta_z={:e}, k_corner={:d})", t_z(x, y, z), x, y,
              z, delta_z(x, y, z), k_corner(x, y, z));
        }
      }
    }
  }
}

void XZLagrange4pt::calcWeights(const Field3D &delta_x, const Field3D &delta_z,
                                const BoutMask &mask) {
  skip_mask = mask;
  calcWeights(delta_x, delta_z);
}

Field3D XZLagrange4pt::interpolate(const Field3D &f) const {

  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  for (int x = localmesh->xstart; x <= localmesh->xend; x++) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      for (int z = 0; z < localmesh->LocalNz; z++) {

        if (skip_mask(x, y, z))
          continue;

        int jx = i_corner(x, y, z);
        int jx2mnew = (jx == 0) ? 0 : (jx - 1);
        int jxpnew = jx + 1;
        int jx2pnew = (jx == (localmesh->LocalNx - 2)) ? jxpnew : (jxpnew + 1);

        int ncz = localmesh->LocalNz;

        // Get the 4 Z points
        int jz = ((k_corner(x, y, z) % ncz) + ncz) % ncz;

        int jzpnew = (jz + 1) % ncz;
        int jz2pnew = (jz + 2) % ncz;
        int jz2mnew = (jz - 1 + ncz) % ncz;

        // Interpolate in Z first
        BoutReal xvals[4];

        int y_next = y + y_offset;

        xvals[0] = lagrange_4pt(f(jx2mnew, y_next, jz2mnew), f(jx2mnew, y_next, jz),
                                f(jx2mnew, y_next, jzpnew), f(jx2mnew, y_next, jz2pnew),
                                t_z(x, y, z));

        xvals[1] =
            lagrange_4pt(f(jx, y_next, jz2mnew), f(jx, y_next, jz), f(jx, y_next, jzpnew),
                         f(jx, y_next, jz2pnew), t_z(x, y, z));

        xvals[2] = lagrange_4pt(f(jxpnew, y_next, jz2mnew), f(jxpnew, y_next, jz),
                                f(jxpnew, y_next, jzpnew), f(jxpnew, y_next, jz2pnew),
                                t_z(x, y, z));

        xvals[3] = lagrange_4pt(f(jx2pnew, y_next, jz2mnew), f(jx2pnew, y_next, jz),
                                f(jx2pnew, y_next, jzpnew), f(jx2pnew, y_next, jz2pnew),
                                t_z(x, y, z));

        // Then in X
        f_interp(x, y_next, z) = lagrange_4pt(xvals, t_x(x, y, z));
      }
    }
  }
  return f_interp;
}

Field3D XZLagrange4pt::interpolate(const Field3D &f, const Field3D &delta_x,
                                   const Field3D &delta_z) {
  calcWeights(delta_x, delta_z);
  return interpolate(f);
}

Field3D XZLagrange4pt::interpolate(const Field3D &f, const Field3D &delta_x,
                                   const Field3D &delta_z, const BoutMask &mask) {
  calcWeights(delta_x, delta_z, mask);
  return interpolate(f);
}

// 4-point Lagrangian interpolation
// offset must be between 0 and 1
BoutReal XZLagrange4pt::lagrange_4pt(const BoutReal v2m, const BoutReal vm,
                                     const BoutReal vp, const BoutReal v2p,
                                     const BoutReal offset) const {
  return -offset * (offset - 1.0) * (offset - 2.0) * v2m / 6.0 +
         0.5 * (offset * offset - 1.0) * (offset - 2.0) * vm -
         0.5 * offset * (offset + 1.0) * (offset - 2.0) * vp +
         offset * (offset * offset - 1.0) * v2p / 6.0;
}

BoutReal XZLagrange4pt::lagrange_4pt(const BoutReal v[], const BoutReal offset) const {
  return lagrange_4pt(v[0], v[1], v[2], v[3], offset);
}
