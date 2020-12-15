/**************************************************************************
 * Base class for fields
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

//#include <globals.hxx>

#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <boutexception.hxx>
#include <field.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <utils.hxx>

namespace bout {
/// Make sure \p location is a sensible value for \p mesh
///
/// Throws if checks are enabled and trying to use a staggered
/// location on a non-staggered mesh
CELL_LOC normaliseLocation(CELL_LOC location, Mesh* mesh) {
  AUTO_TRACE();

  // CELL_DEFAULT always means CELL_CENTRE
  if (location == CELL_DEFAULT) {
    return CELL_CENTRE;
  }

  // No mesh means we can't check if we're using staggered grids, so
  // we'll have to trust the user in this case. This can happen if
  // we're making a field before the global mesh has been initialised
  // -- probably not good, but possible.
  if (mesh == nullptr) {
    return location;
  }

  if (mesh->StaggerGrids) {
    if (location == CELL_VSHIFT) {
      throw BoutException(
          "Field: CELL_VSHIFT cell location only makes sense for vectors");
    }
    return location;
  } else {
#if CHECK > 0
    if (location != CELL_CENTRE) {
      throw BoutException("Field: Trying to set off-centre location on "
                          "non-staggered grid\n"
                          "         Did you mean to enable staggered grids?");
    }
#endif
    return CELL_CENTRE;
  }
}
} // namespace bout

Field::Field(Mesh* localmesh, CELL_LOC location_in, DirectionTypes directions_in)
    : fieldmesh(localmesh == nullptr ? bout::globals::mesh : localmesh),
      location(bout::normaliseLocation(location_in, fieldmesh)),
      directions(directions_in) {
}

void Field::setLocation(CELL_LOC new_location) {
  AUTO_TRACE();

  location = bout::normaliseLocation(new_location, getMesh());

  fieldCoordinates = nullptr;
}

CELL_LOC Field::getLocation() const {
  AUTO_TRACE();
  return location;
}

Coordinates *Field::getCoordinates() const {
  if (fieldCoordinates) {
    return fieldCoordinates.get();
  } else {
    fieldCoordinates = getMesh()->getCoordinatesSmart(getLocation());
    return fieldCoordinates.get();
  }
}

Coordinates *Field::getCoordinates(CELL_LOC loc) const {
  if (loc == CELL_DEFAULT) return getCoordinates();  
  return getMesh()->getCoordinates(loc);
}

int Field::getNx() const{
  return getMesh()->LocalNx;
}

int Field::getNy() const{
  return getMesh()->LocalNy;
}

int Field::getNz() const{
  return getMesh()->LocalNz;
}
