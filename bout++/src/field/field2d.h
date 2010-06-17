/*!
 * \file field2d.h
 *
 * \brief Definition of 2D scalar field class
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
class Field2D;

#ifndef __FIELD2D_H__
#define __FIELD2D_H__

#include "field.h"
#include "field_data.h"
#include "field3d.h"
#include "fieldperp.h"
#include "stencils.h"

/*!
 * \brief 2D X-Y scalar fields
 *
 * Handles data for axisymmetric quantities. Essentially the same
 * as the Field3D class.
 */
class Field2D : public Field, public FieldData {
 public:
  Field2D();
  Field2D(const Field2D& f);
  Field2D(real val);
  ~Field2D();

  Field2D* clone() const;

  real **getData() const; // Remove this!
  
  /// Ensure data is allocated
  void Allocate();
  bool isAllocated() { return data !=  NULL; } ///< Test if data is allocated

  // Operators

  Field2D & operator=(const Field2D &rhs);
  Field2D & operator=(const real rhs);

  real* operator[](int jx) const;

  Field2D & operator+=(const Field2D &rhs);
  Field2D & operator+=(const real rhs);
  Field2D & operator-=(const Field2D &rhs);
  Field2D & operator-=(const real rhs);
  Field2D & operator*=(const Field2D &rhs);
  Field2D & operator*=(const real rhs);
  Field2D & operator/=(const Field2D &rhs);
  Field2D & operator/=(const real rhs);
  Field2D & operator^=(const Field2D &rhs);
  Field2D & operator^=(const real rhs);
  
  // Binary operators

  const Field2D operator+(const Field2D &other) const;
  const Field2D operator+(const real rhs) const;
  const Field2D operator-() const;
  const Field2D operator-(const Field2D &other) const;
  const Field2D operator-(const real rhs) const;
  const Field2D operator*(const Field2D &other) const;
  const Field2D operator*(const real rhs) const;
  const Field2D operator/(const Field2D &other) const;
  const Field2D operator/(const real rhs) const;
  const Field2D operator^(const Field2D &other) const;
  const Field2D operator^(const real rhs) const;

  // Left binary operators

  const Field3D operator+(const Field3D &other) const;
  const Field3D operator-(const Field3D &other) const;
  const Field3D operator*(const Field3D &other) const;
  const Field3D operator/(const Field3D &other) const;
  const Field3D operator^(const Field3D &other) const;

  const FieldPerp operator+(const FieldPerp &other) const;
  const FieldPerp operator-(const FieldPerp &other) const;
  const FieldPerp operator*(const FieldPerp &other) const;
  const FieldPerp operator/(const FieldPerp &other) const;
  friend const Field2D operator/(const real lhs, const Field2D &rhs);
  const FieldPerp operator^(const FieldPerp &other) const;
  friend const Field2D operator^(const real lhs, const Field2D &rhs);
    
  // Stencils

  void getXarray(int y, int z, rvec &xv) const;
  void getYarray(int x, int z, rvec &yv) const;
  void getZarray(int x, int y, rvec &zv) const;

  void setXarray(int y, int z, const rvec &xv);
  void setYarray(int x, int z, const rvec &yv);

  void SetStencil(bstencil *fval, bindex *bx) const;
  void SetXStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void SetYStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void SetZStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;

  // Functions
  
  const Field2D Sqrt() const;
  const Field2D Abs() const;
  real Min(bool allpe=false) const;
  real Max(bool allpe=false) const;
  bool Finite() const;
  
  friend const Field2D sin(const Field2D &f);
  friend const Field2D cos(const Field2D &f);
  friend const Field2D tan(const Field2D &f);

  friend const Field2D sinh(const Field2D &f);
  friend const Field2D cosh(const Field2D &f);
  friend const Field2D tanh(const Field2D &f);

  bool is_const;
  real value;

  // FieldData virtual functions
  
  bool isReal() const   { return true; }         // Consists of real values
  bool is3D() const     { return false; }        // Field is 2D
  int  byteSize() const { return sizeof(real); } // Just one real
  int  realSize() const { return 1; }
  int  getData(int x, int y, int z, void *vptr) const;
  int  getData(int x, int y, int z, real *rptr) const;
  int  setData(int x, int y, int z, void *vptr);
  int  setData(int x, int y, int z, real *rptr);

  bool ioSupport() { return true; } ///< This class supports I/O operations
  real *getData(int component) { 
    real **d = getData();
    return *d;
  }
  void zeroComponent(int component){
    *this = 0.0;
  }

#ifdef CHECK
  bool check_data(bool vital = true) const; ///< Checks if the data is all valid.
  void doneComms() {bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true; }
#endif
  
 private:
  real **data;

  // Data stack: Blocks of memory for this class
  static int nblocks, max_blocks;
  static real ***block; // Pointer to blocks of memory

  void alloc_data();
  void free_data();
};

// Non-member overloaded operators

const Field2D operator+(const real lhs, const Field2D &rhs);
const Field2D operator-(const real lhs, const Field2D &rhs);
const Field2D operator*(const real lhs, const Field2D &rhs);
const Field2D operator/(const real lhs, const Field2D &rhs);
const Field2D operator^(const real lhs, const Field2D &rhs);

// Non-member functions
const Field2D sqrt(const Field2D &f);
const Field2D abs(const Field2D &f);
real min(const Field2D &f, bool allpe=false);
real max(const Field2D &f, bool allpe=false);
bool finite(const Field2D &f);

#endif /* __FIELD2D_H__ */
