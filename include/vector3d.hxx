/*!
 * \file vector2d.hxx
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

#pragma once
#ifndef __VECTOR3D_H__
#define __VECTOR3D_H__

class Field2D; //#include "field2d.hxx"
#include "field3d.hxx"

class Vector2D; //#include "vector2d.hxx"

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
  BoutReal operator=(const BoutReal val);
  
  // Operators
  Vector3D & operator+=(const Vector3D &rhs);
  Vector3D & operator+=(const Vector2D &rhs);

  const Vector3D operator-() const;
  Vector3D & operator-=(const Vector3D &rhs);
  Vector3D & operator-=(const Vector2D &rhs);
  
  Vector3D & operator*=(const BoutReal rhs);
  Vector3D & operator*=(const Field2D &rhs);
  Vector3D & operator*=(const Field3D &rhs);
  
  Vector3D & operator/=(const BoutReal rhs);
  Vector3D & operator/=(const Field2D &rhs);
  Vector3D & operator/=(const Field3D &rhs);

  Vector3D & operator^=(const Vector3D &rhs); // Cross product
  Vector3D & operator^=(const Vector2D &rhs);
  
  // Binary operators

  const Vector3D operator+(const Vector3D &rhs) const;
  const Vector3D operator+(const Vector2D &rhs) const;

  const Vector3D operator-(const Vector3D &rhs) const;
  const Vector3D operator-(const Vector2D &rhs) const;

  const Vector3D operator*(const BoutReal rhs) const;
  const Vector3D operator*(const Field2D &rhs) const;
  const Vector3D operator*(const Field3D &rhs) const;

  const Vector3D operator/(const BoutReal rhs) const;
  const Vector3D operator/(const Field2D &rhs) const;
  const Vector3D operator/(const Field3D &rhs) const;

  const Field3D operator*(const Vector3D &rhs) const; // Dot product
  const Field3D operator*(const Vector2D &rhs) const;
  
  const Vector3D operator^(const Vector3D &rhs) const; // Cross product
  const Vector3D operator^(const Vector2D &rhs) const;
  
  void setLocation(CELL_LOC loc); // Set variable location
  
  // Z shifting
  void shiftZ(int jx, int jy, double zangle) {
      x.shiftZ(jx, jy, zangle);
      y.shiftZ(jx, jy, zangle);
      z.shiftZ(jx, jy, zangle);
  }
  const Vector3D shiftZ(const BoutReal zangle) const;

  // Non-member functions
  friend const Field3D abs(const Vector3D &v);

  // FieldData virtual functions
  
  bool isReal() const   { return true; }
  bool is3D() const     { return true; }
  int  byteSize() const { return 3*sizeof(BoutReal); }
  int  BoutRealSize() const { return 3; }
  int  getData(int jx, int jy, int jz, void *vptr) const;
  int  getData(int jx, int jy, int jz, BoutReal *rptr) const;
  int  setData(int jx, int jy, int jz, void *vptr);
  int  setData(int jx, int jy, int jz, BoutReal *rptr);
  
  void applyBoundary(bool init=false);
  void applyTDerivBoundary();
 private:
  Vector3D *deriv; ///< Time-derivative, can be NULL
};

// Non-member overloaded operators

const Vector3D operator*(const BoutReal lhs, const Vector3D &rhs);
const Vector3D operator*(const Field2D &lhs, const Vector3D &rhs);
const Vector3D operator*(const Field3D &lhs, const Vector3D &rhs);

/*!
 * @brief Time derivative of 3D vector field
 */
inline Vector3D& ddt(Vector3D &f) {
  return *(f.timeDeriv());
}

#endif // __VECTOR3D_H__

