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

#include <stdarg.h>

#include <field.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <bout/mesh.hxx>

Field::Field() : fieldmesh(nullptr), fieldCoordinates(nullptr) {
#if CHECK > 0
  bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true;
#endif
}

Field::Field(Mesh *localmesh) : fieldmesh(localmesh), fieldCoordinates(nullptr) {
  if (fieldmesh == nullptr) {
    fieldmesh = mesh;
  }

/// Note we would like to do `fieldCoordinates = getCoordinates();` here but can't
/// currently as this would lead to circular/recursive behaviour (getCoordinates would
/// call fieldmesh->coordinates, which would create fields, which would then call
/// getCoordinates again etc.).
#if CHECK > 0
  bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true;
#endif
}

Coordinates *Field::getCoordinates() const {
  if (fieldCoordinates) {
    return fieldCoordinates;    
  } else {
    /// Note we don't set fieldCoordinates here as the routine would become
    /// non-const limiting our ability to call it in places like derivatives
    /// where fields are passed as const.
    fieldCoordinates = getMesh()->coordinates(getLocation());
    return fieldCoordinates;
  }
}

Coordinates *Field::getCoordinates(CELL_LOC loc) const {
  return getMesh()->coordinates(loc);
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

/////////////////// PROTECTED ////////////////////


// Report an error occurring
void Field::error(const char *s, ...) const {
  int buf_len=512;
  char * err_buffer=new char[buf_len];

  if (s == nullptr) {
    output_error.write("Unspecified error in field\n");
  } else {

    bout_vsnprintf(err_buffer,buf_len, s);

#ifdef TRACK
      output_error.write("Error in '%s': %s", name.c_str(), err_buffer);
#else
      output_error.write("Error in field: %s", err_buffer);
#endif
  }
  std::string msg="Error in field: ";
  msg+=err_buffer;
  delete[] err_buffer;
  throw BoutException(msg);
}

