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

#include "bout/deprecated.hxx"

/*!
 * \brief Base class for fields
 *
 * Defines the virtual function SetStencil, used by differencing methods
 */
class Field {
 public:
  Field();
  virtual ~Field() { }
  
  /// Clone function
  /*!
    C++ doesn't support statements like:
      Field3D a;
      a = ...
      Field b = a;  // Error
    The clone function implements a second-best
      Field3D a;
      a = ...
      Field *b = a.clone();
    This is useful where a const field must be modified within a function
    and is used in the deriv.cpp functions.
  */
  DEPRECATED(virtual Field* clone() const) = 0;

  virtual void shiftToReal(bool toBoutReal) {
    // Does nothing by default
  }

  virtual void setStencil(bstencil *val, bindex *bx) const = 0;

  // These routines only set a stencil in one dimension
  // Should be faster, and replaces the above SetStencil function.
  virtual void setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const = 0;
  virtual void setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const = 0;
  virtual void setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const = 0;

  virtual void setLocation(CELL_LOC loc) { }
  virtual CELL_LOC getLocation() const {
    return CELL_CENTRE;
  }

  virtual void getXArray(int y, int z, rvec &xv) const {
    error("Field: Base class does not implement getXarray");
  }
  virtual void getYArray(int x, int z, rvec &yv) const {
    error("Field: Base class does not implement getYarray");
  }
  virtual void getZArray(int x, int y, rvec &zv) const {
    error("Field: Base class does not implement getZarray");
  }

  virtual void setXArray(int y, int z, const rvec &xv) {
    error("Field: Base class does not implement setXarray");
  }
  virtual void setYArray(int x, int z, const rvec &yv) {
    error("Field: Base class does not implement setYarray");
  }
  virtual void setZArray(int x, int y, const rvec &zv) {
    error("Field: Base class does not implement setZarray");
  }
    
#ifdef TRACK
  string getName() const { return name; }
  void setName(string s) { name = s; }

  string name;
#endif

#ifdef CHECK
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
 protected:
  /// Supplies an error method. Currently just prints and exits, but
  /// should do something more cunning...
  void error(const char *s, ...) const;
};

#ifndef ddt

/// Concise way to write time-derivatives
#define ddt(f) (*((f).timeDeriv()))

#endif // ddt

#endif /* __FIELD_H__ */
