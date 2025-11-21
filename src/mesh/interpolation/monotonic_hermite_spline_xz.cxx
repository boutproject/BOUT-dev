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

#include "bout/globals.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/mesh.hxx"

#include <vector>

Field3D XZMonotonicHermiteSpline::interpolate(const Field3D& f,
                                              const std::string& region) const {
  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp(f.getMesh());
  f_interp.allocate();

  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  Field3D fx = bout::derivatives::index::DDX(f, CELL_DEFAULT, "DEFAULT");
  Field3D fz = bout::derivatives::index::DDZ(f, CELL_DEFAULT, "DEFAULT", "RGN_ALL");
  localmesh->communicate_no_slices(fx, fz);
  Field3D fxz = bout::derivatives::index::DDX(fz, CELL_DEFAULT, "DEFAULT");
  localmesh->communicate_no_slices(fxz);

  const auto curregion{getRegion(region)};
  BOUT_FOR(i, curregion) {
    const auto iyp = i.yp(y_offset);

    const auto ic = i_corner[i];
    const auto iczp = ic.zp();
    const auto icxp = ic.xp();
    const auto icxpzp = iczp.xp();

    // Interpolate f in X at Z
    const BoutReal f_z =
        f[ic] * h00_x[i] + f[icxp] * h01_x[i] + fx[ic] * h10_x[i] + fx[icxp] * h11_x[i];

    // Interpolate f in X at Z+1
    const BoutReal f_zp1 = f[iczp] * h00_x[i] + f[icxpzp] * h01_x[i] + fx[iczp] * h10_x[i]
                           + fx[icxpzp] * h11_x[i];

    // Interpolate fz in X at Z
    const BoutReal fz_z = fz[ic] * h00_x[i] + fz[icxp] * h01_x[i] + fxz[ic] * h10_x[i]
                          + fxz[icxp] * h11_x[i];

    // Interpolate fz in X at Z+1
    const BoutReal fz_zp1 = fz[iczp] * h00_x[i] + fz[icxpzp] * h01_x[i]
                            + fxz[iczp] * h10_x[i] + fxz[icxpzp] * h11_x[i];

    // Interpolate in Z
    BoutReal result =
        +f_z * h00_z[i] + f_zp1 * h01_z[i] + fz_z * h10_z[i] + fz_zp1 * h11_z[i];

    ASSERT2(std::isfinite(result) || i.x() < localmesh->xstart
            || i.x() > localmesh->xend);

    // Monotonicity
    // Force the interpolated result to be in the range of the
    // neighbouring cell values. This prevents unphysical overshoots,
    // but also degrades accuracy near maxima and minima.
    // Perhaps should only impose near boundaries, since that is where
    // problems most obviously occur.
    const BoutReal localmax = BOUTMAX(f[ic], f[icxp], f[iczp], f[icxpzp]);

    const BoutReal localmin = BOUTMIN(f[ic], f[icxp], f[iczp], f[icxpzp]);

    ASSERT2(std::isfinite(localmax) || i.x() < localmesh->xstart
            || i.x() > localmesh->xend);
    ASSERT2(std::isfinite(localmin) || i.x() < localmesh->xstart
            || i.x() > localmesh->xend);

    if (result > localmax) {
      result = localmax;
    }
    if (result < localmin) {
      result = localmin;
    }

    f_interp[iyp] = result;
  }
  return f_interp;
}
