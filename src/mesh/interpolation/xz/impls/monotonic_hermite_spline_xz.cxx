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

#include "monotonic_hermite_spline_xz.hxx"

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/globals.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"

#include <algorithm>
#include <cmath>
#include <string>

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

    auto [ic, iczp, icxp, icxpzp, result] = interpolate_point(f, fx, fz, fxz, i);
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

    const auto diff = ((localmax - localmin) * rtol) + atol;

    result = std::min(result, localmax + diff);
    result = std::max(result, localmin - diff);

    f_interp[iyp] = result;
  }
  return f_interp;
}
