/*!
 * \file vector3d.hxx
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

/*!
 * Represents a 3D vector, with x,y,z components
 * stored as separate Field3D objects
 *
 * Example
 * -------
 *
 * Vector3D f;
 * 
 * a.x; // Error! a.x not allocated
 *
 * 
 */ 
class Vector3D : public FieldData {
 public:
  /*!
   * Constructor. Just sets covariant = true and deriv = NULL
   *
   * Does not initialise any of the fields
   */
  Vector3D();
  
  /*!
   * Copy constructor. After this the components (x,y,z)
   * will refer to the same data as f.(x,y,z)
   */
  Vector3D(const Vector3D &f);

  /*!
   * Destructor. If the time derivative has been
   * used, then some book-keeping is needed to ensure
   * that fields are only destroyed once.
   */
  ~Vector3D();

  /*!
   * The components of the vector. These can be 
   * either co- or contra-variant, depending on
   * the boolean flag "covariant"
   */ 
  Field3D x, y, z;

  /*!
   * Flag to specify whether the components (x,y,z)
   * are co- or contra-variant.
   * 
   * true if the components are covariant (default)
   * false if the components are contravariant
   *
   * Conversion between forms should be done by calling
   * the toContravariant and toCovariant methods.
   *
   * Only modify this variable directly if you know what you are doing!
   * 
   */ 
  bool covariant;

  /*!
   * In-place conversion to covariant form. 
   * 
   * If already covariant (covariant = true) then does nothing
   * If contravariant, multiplies by metric tensor g_{ij}
   */
  void toCovariant();
  
  /*!
   * In-place conversion to contravariant form
   *
   * If already contravariant (covariant = false) then does nothing
   * If covariant, multiplies by metric tensor g^{ij}
   * 
   */
  void toContravariant();
  
  /*!
   * Return a pointer to the time-derivative field
   *
   * The first time this is called, a new Vector3D object is created.
   * Subsequent calls return a pointer to this same object
   *
   * For convenience, a standalone function "ddt" exists, so that 
   *
   * ddt(v) is equivalent to *(v.timeDeriv())
   * 
   * This does some book-keeping to ensure that the time derivative
   * of the components is the same as the components of the time derivative
   *
   * ddt(v).x == ddt(v.x) 
   */
  Vector3D* timeDeriv();

  // Assignment
  Vector3D & operator=(const Vector3D &rhs);
  Vector3D & operator=(const Vector2D &rhs);
  BoutReal operator=(BoutReal val);
  
  // Operators
  Vector3D & operator+=(const Vector3D &rhs);
  Vector3D & operator+=(const Vector2D &rhs);

  const Vector3D operator-() const;
  Vector3D & operator-=(const Vector3D &rhs);
  Vector3D & operator-=(const Vector2D &rhs);
  
  Vector3D & operator*=(BoutReal rhs);
  Vector3D & operator*=(const Field2D &rhs);
  Vector3D & operator*=(const Field3D &rhs);
  
  Vector3D & operator/=(BoutReal rhs);
  Vector3D & operator/=(const Field2D &rhs);
  Vector3D & operator/=(const Field3D &rhs);

  Vector3D & operator^=(const Vector3D &rhs); // Cross product
  Vector3D & operator^=(const Vector2D &rhs);
  
  // Binary operators

  const Vector3D operator+(const Vector3D &rhs) const;
  const Vector3D operator+(const Vector2D &rhs) const;

  const Vector3D operator-(const Vector3D &rhs) const;
  const Vector3D operator-(const Vector2D &rhs) const;

  const Vector3D operator*(BoutReal rhs) const;
  const Vector3D operator*(const Field2D &rhs) const;
  const Vector3D operator*(const Field3D &rhs) const;

  const Vector3D operator/(BoutReal rhs) const;
  const Vector3D operator/(const Field2D &rhs) const;
  const Vector3D operator/(const Field3D &rhs) const;

  const Field3D operator*(const Vector3D &rhs) const; // Dot product
  const Field3D operator*(const Vector2D &rhs) const;
  
  /*!
   * Cross product
   *
   * Note: This operator has low precedence in C++,
   *       lower than multiplication for example
   */ 
  const Vector3D operator^(const Vector3D &rhs) const; 
  
  /*!
   * Cross product
   *
   * Note: This operator has low precedence in C++,
   *       lower than multiplication for example
   */
  const Vector3D operator^(const Vector2D &rhs) const;
  
  /*!
   * Set variable cell location
   */ 
  void setLocation(CELL_LOC loc); 

  /// Visitor pattern support
  void accept(FieldVisitor &v) override;
  
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

const Vector3D operator*(BoutReal lhs, const Vector3D &rhs);
const Vector3D operator*(const Field2D &lhs, const Vector3D &rhs);
const Vector3D operator*(const Field3D &lhs, const Vector3D &rhs);

/*!
 * Absolute magnitude (modulus) of a vector  |v|
 * 
 * sqrt( v.x^2 + v.y^2 + v.z^2 )
 */ 
const Field3D abs(const Vector3D &v);

/*!
 * @brief Time derivative of 3D vector field
 */
inline Vector3D& ddt(Vector3D &f) {
  return *(f.timeDeriv());
}

#endif // __VECTOR3D_H__

