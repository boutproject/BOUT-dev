/**************************************************************************
 * Copyright 2018 B.D.Dudson, P. Hill
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

#include "globals.hxx"
#include "interpolation.hxx"
#include "output.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/mesh.hxx"

#include <vector>

Field3D MonotonicHermiteSpline::interpolate(const Field3D &f,
                                            const std::string& region) const {
  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp(f.getMesh());
  f_interp.allocate();

  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  Field3D fx = bout::derivatives::index::DDX(f, CELL_DEFAULT, "DEFAULT");
  localmesh->communicateXZ(fx);
  Field3D fz = bout::derivatives::index::DDZ(f, CELL_DEFAULT, "DEFAULT", "RGN_WITH_XBNDRIES");
  localmesh->communicateXZ(fz);
  Field3D fxz = bout::derivatives::index::DDX(fz, CELL_DEFAULT, "DEFAULT");
  localmesh->communicateXZ(fxz);

  BOUT_FOR(i, f.getRegion(region)) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    if (skip_mask(x, y, z))
      continue;

    // Due to lack of guard cells in z-direction, we need to ensure z-index
    // wraps around
    int ncz = localmesh->LocalNz;
    int z_mod = ((k_corner(x, y, z) % ncz) + ncz) % ncz;
    int z_mod_p1 = (z_mod + 1) % ncz;

    int y_next = y + y_offset;

    // Interpolate f in X at Z
    BoutReal f_z = f(i_corner(x, y, z), y_next, z_mod) * h00_x(x, y, z) +
      f(i_corner(x, y, z) + 1, y_next, z_mod) * h01_x(x, y, z) +
      fx(i_corner(x, y, z), y_next, z_mod) * h10_x(x, y, z) +
      fx(i_corner(x, y, z) + 1, y_next, z_mod) * h11_x(x, y, z);

    // Interpolate f in X at Z+1
    BoutReal f_zp1 = f(i_corner(x, y, z), y_next, z_mod_p1) * h00_x(x, y, z) +
      f(i_corner(x, y, z) + 1, y_next, z_mod_p1) * h01_x(x, y, z) +
      fx(i_corner(x, y, z), y_next, z_mod_p1) * h10_x(x, y, z) +
      fx(i_corner(x, y, z) + 1, y_next, z_mod_p1) * h11_x(x, y, z);

    // Interpolate fz in X at Z
    BoutReal fz_z = fz(i_corner(x, y, z), y_next, z_mod) * h00_x(x, y, z) +
      fz(i_corner(x, y, z) + 1, y_next, z_mod) * h01_x(x, y, z) +
      fxz(i_corner(x, y, z), y_next, z_mod) * h10_x(x, y, z) +
      fxz(i_corner(x, y, z) + 1, y_next, z_mod) * h11_x(x, y, z);

    // Interpolate fz in X at Z+1
    BoutReal fz_zp1 = fz(i_corner(x, y, z), y_next, z_mod_p1) * h00_x(x, y, z) +
      fz(i_corner(x, y, z) + 1, y_next, z_mod_p1) * h01_x(x, y, z) +
      fxz(i_corner(x, y, z), y_next, z_mod_p1) * h10_x(x, y, z) +
      fxz(i_corner(x, y, z) + 1, y_next, z_mod_p1) * h11_x(x, y, z);

    // Interpolate in Z
    BoutReal result = +f_z * h00_z(x, y, z) + f_zp1 * h01_z(x, y, z) +
      fz_z * h10_z(x, y, z) + fz_zp1 * h11_z(x, y, z);

    ASSERT2(finite(result) || x < localmesh->xstart || x > localmesh->xend);

    // Monotonicity
    // Force the interpolated result to be in the range of the
    // neighbouring cell values. This prevents unphysical overshoots,
    // but also degrades accuracy near maxima and minima.
    // Perhaps should only impose near boundaries, since that is where
    // problems most obviously occur.
    BoutReal localmax = BOUTMAX(f(i_corner(x, y, z), y_next, z_mod),
        f(i_corner(x, y, z)+1, y_next, z_mod),
        f(i_corner(x, y, z), y_next, z_mod_p1),
        f(i_corner(x, y, z)+1, y_next, z_mod_p1));

    BoutReal localmin = BOUTMIN(f(i_corner(x, y, z), y_next, z_mod),
        f(i_corner(x, y, z)+1, y_next, z_mod),
        f(i_corner(x, y, z), y_next, z_mod_p1),
        f(i_corner(x, y, z)+1, y_next, z_mod_p1));

    ASSERT2(finite(localmax) || x < localmesh->xstart || x > localmesh->xend);
    ASSERT2(finite(localmin) || x < localmesh->xstart || x > localmesh->xend);

    if (result > localmax) {
      result = localmax;
    }
    if (result < localmin) {
      result = localmin;
    }

    f_interp(x, y_next, z) = result;

  }
  return f_interp;
}
