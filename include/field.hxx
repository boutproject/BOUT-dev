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

#include <stdio.h>

#include "bout_types.hxx"
#include "stencils.hxx"
#include <bout/rvec.hxx>
#include "boutexception.hxx"

#include "bout/deprecated.hxx"

#include "bout/dataiterator.hxx"

#include "unused.hxx"

class Mesh;
extern Mesh * mesh;

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
  Field();
  Field(Mesh * localmesh);
  virtual ~Field() { }

  // Data access
  virtual const BoutReal& operator[](const Indices &i) const = 0;

  virtual void setLocation(CELL_LOC loc) {
    if (loc != CELL_CENTRE)
      throw BoutException("not implemented!");
  }
  virtual CELL_LOC getLocation() const {
    return CELL_CENTRE;
  }

#ifdef TRACK
  std::string getName() const { return name; }
  void setName(std::string s) { name = s; }
#else
  std::string getName() const { return ""; }
  void setName(std::string UNUSED(s)) {}
#endif
  std::string name;

#if CHECK > 0
  // Routines to test guard/boundary cells set
  
  virtual bool bndryValid() {
    if(!bndry_xin)
      error("Inner X guard cells not set\n");
    if(!bndry_xout)
      error("Outer X guard cells not set\n");
    if(!bndry_yup)
      error("Upper y guard cells not set\n");
    if(!bndry_ydown)
      error("Lower y guard cells not set\n");
    return true;
  }
  
  // Status of the 4 boundaries
  bool bndry_xin, bndry_xout, bndry_yup, bndry_ydown;
#endif
  virtual Mesh * getMesh() const{
    if (fieldmesh){
      return fieldmesh;
    } else {
      return mesh;
    }
  }
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

  /// Make region mendatory for all fields
  virtual const IndexRange region(REGION rgn) const = 0;
 protected:
  Mesh * fieldmesh;
  /// Supplies an error method. Currently just prints and exits, but
  /// should do something more cunning...
  void error(const char *s, ...) const;
};

#endif /* __FIELD_H__ */
