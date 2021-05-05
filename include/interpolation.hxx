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

#include "bout/mesh.hxx"

/// Perform interpolation between centre -> shifted or vice-versa
/*!
  Interpolate using 4th-order staggered formula

  @param[in] s  Input stencil. mm -> -3/2, m -> -1/2, p -> +1/2, pp -> +3/2
*/
inline BoutReal interp(const stencil& s) {
  return (9. * (s.m + s.p) - s.mm - s.pp) / 16.;
}

/// Interpolate to a give cell location
/*!
  Interpolate between different cell locations

  NOTE: This requires communication if the result is required in guard cells
  NOTE: Since corner guard cells cannot be communicated, it never makes sense
  to calculate interpolation in guard cells. If guard cell values are required,
  we must communicate (unless interpolating in z). Since mesh->communicate()
  communicates both x- and y-guard cells by default, there is no difference
  between RGN_ALL, RGN_NOX and RGN_NOY.

  @param[in]   var  Input variable
  @param[in]   loc  Location of output values
  @param[in]   region  Region where output will be calculated
*/
template <typename T>
const T interp_to(const T& var, CELL_LOC loc, const std::string region = "RGN_ALL") {
  AUTO_TRACE();
  static_assert(bout::utils::is_Field2D<T>::value || bout::utils::is_Field3D<T>::value,
                "interp_to must be templated with one of Field2D or Field3D.");
  ASSERT1(loc != CELL_DEFAULT); // doesn't make sense to interplote to CELL_DEFAULT

  Mesh* fieldmesh = var.getMesh();

  if ((loc != CELL_CENTRE) && (fieldmesh->StaggerGrids == false)) {
    throw BoutException("Asked to interpolate, but StaggerGrids is disabled!");
  }

  if (var.getLocation() == loc) {
    // Nothing to do - just return unchanged
    return var;
  }

  // NOTE: invalidateGuards() is called in Field3D::alloctate() if the data
  // block is not already allocated, so will be called here if
  // region==RGN_NOBNDRY
  T result{emptyFrom(var).setLocation(loc)};

  // Staggered grids enabled, and need to perform interpolation
  TRACE("Interpolating {} -> {}", toString(var.getLocation()), toString(loc));

  if (region != "RGN_NOBNDRY") {
    // result is requested in some boundary region(s)
    result = var; // NOTE: This is just for boundaries. FIX!
    result.setLocation(loc); // location gets reset when assigning from var
    result.allocate();
  }

  // Cell location of the input field
  const CELL_LOC location = var.getLocation();

  if ((location == CELL_CENTRE) || (loc == CELL_CENTRE)) {
    // Going between centred and shifted

    // Get the non-centre location for interpolation direction
    const CELL_LOC dir = (loc == CELL_CENTRE) ? location : loc;

    switch (dir) {
    case CELL_XLOW: {
      // At least 2 boundary cells needed for interpolation in x-direction
      ASSERT0(fieldmesh->xstart >= 2);

      if ((location == CELL_CENTRE) && (loc == CELL_XLOW)) { // C2L
        BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
          // Producing a stencil centred around a lower X value
          result[i] = interp(populateStencil<DIRECTION::X, STAGGER::C2L, 2>(var, i));
        }
      } else if (location == CELL_XLOW) { // L2C
        BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
          // Stencil centred around a cell centre
          result[i] = interp(populateStencil<DIRECTION::X, STAGGER::L2C, 2>(var, i));
        }
      }

      break;
    }
    case CELL_YLOW: {
      // At least 2 boundary cells needed for interpolation in y-direction
      ASSERT0(fieldmesh->ystart >= 2);

      // We can't interpolate in y unless we're field-aligned
      const bool is_unaligned = (var.getDirectionY() == YDirectionType::Standard);
      const T var_fa = is_unaligned ? toFieldAligned(var, "RGN_NOX") : var;

      if (not std::is_base_of<Field2D, T>::value) {
        // Field2D is axisymmetric, so YDirectionType::Standard and
        // YDirectionType::Aligned are equivalent, but trying to set
        // YDirectionType::Aligned explicitly is an error
        result.setDirectionY(YDirectionType::Aligned);
      }

      if (region != "RGN_NOBNDRY") {
        // repeat the hack above for boundary points
        // this avoids a duplicate toFieldAligned call if we had called
        // result = toFieldAligned(result)
        // to get the boundary cells
        //
        // result is requested in some boundary region(s)
        result = var_fa; // NOTE: This is just for boundaries. FIX!
        result.setLocation(loc); // location gets reset when assigning from var
        result.allocate();
      }

      if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) { // C2L
        BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
          // Producing a stencil centred around a lower X value
          result[i] =
              interp(populateStencil<DIRECTION::YAligned, STAGGER::C2L, 2>(var_fa, i));
        }
      } else if (location == CELL_YLOW) { // L2C
        BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
          // Stencil centred around a cell centre
          result[i] =
              interp(populateStencil<DIRECTION::YAligned, STAGGER::L2C, 2>(var_fa, i));
        }
      }

      if (is_unaligned) {
        result = fromFieldAligned(result, "RGN_NOBNDRY");
      }

      break;
    }
    case CELL_ZLOW: {

      if ((location == CELL_CENTRE) && (loc == CELL_ZLOW)) { // C2L
        BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
          // Producing a stencil centred around a lower X value
          result[i] = interp(populateStencil<DIRECTION::Z, STAGGER::C2L, 2>(var, i));
        }
      } else if (location == CELL_ZLOW) { // L2C
        BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
          // Stencil centred around a cell centre
          result[i] = interp(populateStencil<DIRECTION::Z, STAGGER::L2C, 2>(var, i));
        }
      }
      break;
    }
    default: {
      // This should never happen
      throw BoutException("Unsupported direction of interpolation\n"
                          " - don't know how to interpolate to {:s}",
                          toString(loc));
    }
    };

    if ((dir != CELL_ZLOW) && (region != "RGN_NOBNDRY")) {
      fieldmesh->communicate(result);
    }

  } else {
    // Shifted -> shifted
    // For now, shift to centre then to final location loc
    // We probably should not rely on this, but it might work if one of the
    // shifts is in the z-direction where guard cells aren't needed.
    result = interp_to(interp_to(var, CELL_CENTRE), loc, region);
  }
  return result;
}

#endif // __INTERP_H__
