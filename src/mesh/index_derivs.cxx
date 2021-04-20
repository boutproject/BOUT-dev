/**************************************************************************
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

#include "bout/build_config.hxx"

#include "bout/traits.hxx"
#include <bout/index_derivs.hxx>
#include <bout/mesh.hxx>
#include <msg_stack.hxx>
#include <unused.hxx>

/*******************************************************************************
 * Helper routines
 *******************************************************************************/

/// Initialise the derivative methods. Must be called before any derivatives are used
void Mesh::derivs_init(Options* options) {
  TRACE("Initialising derivatives");
  // For each direction need to set what the default method is for each type
  // of derivative.
  DerivativeStore<Field3D>::getInstance().initialise(options);
  DerivativeStore<Field2D>::getInstance().initialise(options);
  // Get the fraction of modes filtered out in FFT derivatives
  options->getSection("ddz")->get("fft_filter", fft_derivs_filter, 0.0);
}

STAGGER Mesh::getStagger(const CELL_LOC inloc, const CELL_LOC outloc,
                         const CELL_LOC allowedStaggerLoc) const {
  TRACE("Mesh::getStagger -- three arguments");
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == allowedStaggerLoc)
          || (outloc == allowedStaggerLoc && inloc == CELL_CENTRE));

  if ((!StaggerGrids) || outloc == inloc)
    return STAGGER::None;
  if (outloc == allowedStaggerLoc) {
    return STAGGER::C2L;
  } else {
    return STAGGER::L2C;
  }
}

STAGGER Mesh::getStagger(const CELL_LOC vloc, MAYBE_UNUSED(const CELL_LOC inloc),
                         const CELL_LOC outloc, const CELL_LOC allowedStaggerLoc) const {
  TRACE("Mesh::getStagger -- four arguments");
  ASSERT1(inloc == outloc);
  ASSERT1(vloc == inloc || (vloc == CELL_CENTRE && inloc == allowedStaggerLoc)
          || (vloc == allowedStaggerLoc && inloc == CELL_CENTRE));
  return getStagger(vloc, outloc, allowedStaggerLoc);
}

std::tuple<BoutReal, BoutReal> vUpDown(BoutReal v) {
  return std::tuple<BoutReal, BoutReal>{0.5 * (v + fabs(v)), 0.5 * (v - fabs(v))};
}
