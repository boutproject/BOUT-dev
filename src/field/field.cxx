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

Field::Field(Mesh* localmesh, CELL_LOC location_in, DirectionTypes directions_in)
    : FieldData(localmesh, location_in), directions(directions_in) {}

int Field::getNx() const{
  return getMesh()->LocalNx;
}

int Field::getNy() const{
  return getMesh()->LocalNy;
}

int Field::getNz() const{
  return getMesh()->LocalNz;
}
