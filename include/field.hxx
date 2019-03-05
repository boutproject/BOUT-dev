/*!
 * \file field.hxx
 * \brief field base class definition for differencing methods
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
 */

class Field;

#ifndef __FIELD_H__
#define __FIELD_H__

#include <cstdio>
#include <memory>

#include "bout_types.hxx"
#include "boutexception.hxx"
#include <globals.hxx>
#include "msg_stack.hxx"
#include "stencils.hxx"
#include <bout/rvec.hxx>

#include "unused.hxx"

class Mesh;
class Coordinates;

#ifdef TRACK
#include <string>
#endif

/*!
 * \brief Base class for fields
 *
 * Defines the virtual function SetStencil, used by differencing methods
 */
class Field {
 public:
  Field() = default;

  Field(Mesh* localmesh, CELL_LOC location_in, DIRECTION xDirectionType_in,
      DIRECTION yDirectionType_in, DIRECTION zDirectionType_in);

  // Copy constructor
  Field(const Field& f)
    : name(f.name), fieldmesh(f.fieldmesh),
      fieldCoordinates(f.fieldCoordinates), location(f.location),
      xDirectionType(f.xDirectionType), yDirectionType(f.yDirectionType),
      zDirectionType(f.zDirectionType) {}

  virtual ~Field() { }

  /// Set variable location for staggered grids to @param new_location
  ///
  /// Throws BoutException if new_location is not `CELL_CENTRE` and
  /// staggered grids are turned off and checks are on. If checks are
  /// off, silently sets location to ``CELL_CENTRE`` instead.
  void setLocation(CELL_LOC new_location);
  /// Get variable location
  CELL_LOC getLocation() const;

  /// Getters for DIRECTION types
  DIRECTION getDirectionX() const {
    return xDirectionType;
  }
  DIRECTION getDirectionY() const {
    return yDirectionType;
  }
  DIRECTION getDirectionZ() const {
    return zDirectionType;
  }

  std::string name;

#if CHECK > 0
  // Routines to test guard/boundary cells set
  
  virtual bool bndryValid() {
    if(!bndry_xin)
      throw BoutException("Inner X guard cells not set\n");
    if(!bndry_xout)
      throw BoutException("Outer X guard cells not set\n");
    if(!bndry_yup)
      throw BoutException("Upper y guard cells not set\n");
    if(!bndry_ydown)
      throw BoutException("Lower y guard cells not set\n");
    return true;
  }
  
  /// Status of the 4 boundaries
  bool bndry_xin{true}, bndry_xout{true}, bndry_yup{true}, bndry_ydown{true};
#endif

  virtual Mesh * getMesh() const{
    if (fieldmesh){
      return fieldmesh;
    } else {
      // Don't set fieldmesh=mesh here, so that fieldmesh==nullptr until
      // allocate() is called in one of the derived classes. fieldmesh==nullptr
      // indicates that some initialization that would be done in the
      // constructor if fieldmesh was a valid Mesh object still needs to be
      // done.
      return bout::globals::mesh;
    }
  }

  /// Returns a pointer to the coordinates object at this field's
  /// location from the mesh this field is on.
  virtual Coordinates *getCoordinates() const;
  
  /// Returns a pointer to the coordinates object at the requested
  /// location from the mesh this field is on. If location is CELL_DEFAULT
  /// then return coordinates at field location
  virtual Coordinates *getCoordinates(CELL_LOC loc) const;
  
  /*!
   * Return the number of nx points
   */
  virtual int getNx() const;
  /*!
   * Return the number of ny points
   */
  virtual int getNy() const;
  /*!
   * Return the number of nz points
   */
  virtual int getNz() const;

  friend void swap(Field& first, Field& second) noexcept {
    using std::swap;
    swap(first.name, second.name);
    swap(first.fieldmesh, second.fieldmesh);
    swap(first.fieldCoordinates, second.fieldCoordinates);
    swap(first.location, second.location);
    swap(first.xDirectionType, second.xDirectionType);
    swap(first.yDirectionType, second.yDirectionType);
    swap(first.zDirectionType, second.zDirectionType);
  }

  friend bool fieldsCompatible(const Field& field1, const Field& field2) {
    return
        // The following is a possible alternative to
        // checking the two meshes and location are the same
        // It is slightly stricter if we decide that coordinates
        // could differ in more than just location.
        field1.getCoordinates() == field2.getCoordinates() &&
        // In the unit tests fieldCoordinates get set to nullptr, so we still
        // need to check fieldmesh and location, at least for now
        field1.getMesh() == field2.getMesh() &&
        field1.getLocation() == field2.getLocation() &&
        // Compatible directions
        compatibleDirections(field1.xDirectionType, field2.xDirectionType)
        && compatibleDirections(field1.yDirectionType, field2.yDirectionType)
        && compatibleDirections(field1.zDirectionType, field2.zDirectionType);
  }
protected:
  Mesh* fieldmesh{nullptr};
  mutable std::shared_ptr<Coordinates> fieldCoordinates{nullptr};

  /// Location of the variable in the cell
  CELL_LOC location{CELL_CENTRE};

  /// Set any direction types which are DIRECTION::Null to default values from
  /// fieldmesh.
  void setNullDirectionTypesToDefault();

  /// Copy the members from another Field
  void copyFieldMembers(const Field& f) {
    name = f.name;
    fieldmesh = f.fieldmesh;
    fieldCoordinates = f.fieldCoordinates;
    location = f.location;
    xDirectionType = f.xDirectionType;
    yDirectionType = f.yDirectionType;
    zDirectionType = f.zDirectionType;
  }

  /// Setters for *DirectionType
  void setDirectionX(DIRECTION d) {
    xDirectionType = d;
  }
  void setDirectionY(DIRECTION d) {
    yDirectionType = d;
  }
  void setDirectionZ(DIRECTION d) {
    zDirectionType = d;
  }

private:
  DIRECTION xDirectionType{DIRECTION::Null};
  DIRECTION yDirectionType{DIRECTION::Null};
  DIRECTION zDirectionType{DIRECTION::Null};
};

/// Return an empty shell field of some type derived from Field, with metadata
/// copied but empty data array
template<typename T>
inline T emptyFrom(const T& f) {
  static_assert(std::is_base_of<Field, T>::value, "emptyFrom only works on Fields");
  return T(f.getMesh(), f.getLocation(), f.getDirectionX(), f.getDirectionY(), f.getDirectionZ()).allocate();
}

/// Unary + operator. This doesn't do anything
template<typename T>
T operator+(const T& f) {return f;}

#endif /* __FIELD_H__ */
