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

#include <field.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <bout/mesh.hxx>

Field::Field(Mesh *localmesh)
  : fieldmesh(localmesh==nullptr ? bout::globals::mesh : localmesh) {

// Note we would like to do `fieldCoordinates = getCoordinates();` here but can't
// currently as this would lead to circular/recursive behaviour (getCoordinates would
// call fieldmesh->coordinates, which would create fields, which would then call
// getCoordinates again etc.). This also requires care in the derived class
// constructors.
}

void Field::setLocation(CELL_LOC new_location) {
  AUTO_TRACE();
  if (getMesh()->StaggerGrids) {
    if (new_location == CELL_VSHIFT) {
      throw BoutException(
          "Field: CELL_VSHIFT cell location only makes sense for vectors");
    }
    if (new_location == CELL_DEFAULT) {
      new_location = CELL_CENTRE;
    }

    location = new_location;
  } else {
#if CHECK > 0
    if (new_location != CELL_CENTRE && new_location != CELL_DEFAULT) {
      throw BoutException("Field: Trying to set off-centre location on "
                          "non-staggered grid\n"
                          "         Did you mean to enable staggered grids?");
    }
#endif
    location = CELL_CENTRE;
  }

  // Ensures Coordinates object is initialized for this Field's location
  getCoordinates();
}

CELL_LOC Field::getLocation() const {
  AUTO_TRACE();
  return location;
}

Coordinates *Field::getCoordinates() const {
  if (fieldCoordinates) {
    return fieldCoordinates;    
  } else {
    fieldCoordinates = getMesh()->getCoordinates(getLocation());
    return fieldCoordinates;
  }
}

Coordinates *Field::getCoordinates(CELL_LOC loc) const {
  if (loc == CELL_DEFAULT) return getCoordinates();  
  return getMesh()->getCoordinates(loc);
}

int Field::getNx() const{
  return getMesh()->LocalNx;
};

int Field::getNy() const{
  return getMesh()->LocalNy;
};

int Field::getNz() const{
  return getMesh()->LocalNz;
};

