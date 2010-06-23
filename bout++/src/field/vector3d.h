/*!
 * \file vector2d.h
 *
 * \brief Class for 3D vectors. Built on the Field3D class.
 *
 * \author B. Dudson, October 2007
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

class Vector3D;

#ifndef __VECTOR3D_H__
#define __VECTOR3D_H__

#include "field2d.h"
#include "field3d.h"

#include "vector2d.h"

class Vector3D : public FieldData {
 public:
  Vector3D();
  Vector3D(const Vector3D &f);
  ~Vector3D();

  Field3D x, y, z; // components

  bool covariant; // true if the components are covariant (default)

  void toCovariant();
  void toContravariant();

  /// Return a pointer to the time-derivative field
  Vector3D* timeDeriv();

  // Assignment
  Vector3D & operator=(const Vector3D &rhs);
  Vector3D & operator=(const Vector2D &rhs);
  real operator=(const real val);
  
  // Operators
  Vector3D & operator+=(const Vector3D &rhs);
  Vector3D & operator+=(const Vector2D &rhs);

  const Vector3D operator-() const;
  Vector3D & operator-=(const Vector3D &rhs);
  Vector3D & operator-=(const Vector2D &rhs);
  
  Vector3D & operator*=(const real rhs);
  Vector3D & operator*=(const Field2D &rhs);
  Vector3D & operator*=(const Field3D &rhs);
  
  Vector3D & operator/=(const real rhs);
  Vector3D & operator/=(const Field2D &rhs);
  Vector3D & operator/=(const Field3D &rhs);

  Vector3D & operator^=(const Vector3D &rhs); // Cross product
  Vector3D & operator^=(const Vector2D &rhs);
  
  // Binary operators

  const Vector3D operator+(const Vector3D &rhs) const;
  const Vector3D operator+(const Vector2D &rhs) const;

  const Vector3D operator-(const Vector3D &rhs) const;
  const Vector3D operator-(const Vector2D &rhs) const;

  const Vector3D operator*(const real rhs) const;
  const Vector3D operator*(const Field2D &rhs) const;
  const Vector3D operator*(const Field3D &rhs) const;

  const Vector3D operator/(const real rhs) const;
  const Vector3D operator/(const Field2D &rhs) const;
  const Vector3D operator/(const Field3D &rhs) const;

  const Field3D operator*(const Vector3D &rhs) const; // Dot product
  const Field3D operator*(const Vector2D &rhs) const;
  
  const Vector3D operator^(const Vector3D &rhs) const; // Cross product
  const Vector3D operator^(const Vector2D &rhs) const;

  // Z shifting
  void shiftZ(int jx, int jy, double zangle) {
      x.shiftZ(jx, jy, zangle);
      y.shiftZ(jx, jy, zangle);
      z.shiftZ(jx, jy, zangle);
  }
  const Vector3D shiftZ(const real zangle) const;

  // Non-member functions
  friend const Field3D abs(const Vector3D &v);

  // FieldData virtual functions
  
  bool isReal() const   { return true; }
  bool is3D() const     { return true; }
  int  byteSize() const { return 3*sizeof(real); }
  int  realSize() const { return 3; }
  int  getData(int jx, int jy, int jz, void *vptr) const;
  int  getData(int jx, int jy, int jz, real *rptr) const;
  int  setData(int jx, int jy, int jz, void *vptr);
  int  setData(int jx, int jy, int jz, real *rptr);

  bool ioSupport() { return true; }
  const string getSuffix(int component) const {
    if(covariant) {
      switch(component) {
      case 0:
	return string("_x");
      case 1:
	return string("_y");
      case 2:
	return string("_z");
      }
    }else {
      switch(component) {
      case 0:
	return string("x");
      case 1:
	return string("y");
      case 2:
	return string("z");
      }
    }
    return string("");
  }
  void *getMark() const {
    bool *c = new bool;
    *c = covariant;
    return (void*) c;
  }
  void setMark(void *setting) {
    bool *c = (bool*) setting;
    if(*c) {
      toCovariant();
    }else
      toContravariant();
  }
  real *getData(int component) {
    switch(component) {
    case 0:
      return **(x.getData());
    case 1:
      return **(y.getData());
    case 2:
      return **(z.getData());
    }
    return NULL;
  }
  void zeroComponent(int component) {
    switch(component) {
    case 0:
      x = 0.0;
    case 1:
      y = 0.0;
    case 2:
      z = 0.0;
    }
  }
  
 private:
  Vector3D *ddt; ///< Time-derivative, can be NULL
};

// Non-member overloaded operators

const Vector3D operator*(const real lhs, const Vector3D &rhs);
const Vector3D operator*(const Field2D &lhs, const Vector3D &rhs);
const Vector3D operator*(const Field3D &lhs, const Vector3D &rhs);

#endif // __VECTOR3D_H__

