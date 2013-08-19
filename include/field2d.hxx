/*!
 * \file field2d.hxx
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

#pragma once
#ifndef __FIELD2D_H__
#define __FIELD2D_H__

#include "field.hxx"
#include "field_data.hxx"
class Field3D; //#include "field3d.hxx"
#include "fieldperp.hxx"
#include "stencils.hxx"

#include "bout/deprecated.hxx"

#include <stack>
using std::stack;

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
  Field2D(BoutReal val);
  ~Field2D();

  DEPRECATED(Field2D* clone() const);

  BoutReal **getData() const; // Remove this!
  
  static void cleanup(); // Frees all memory

  /// Ensure data is allocated
  void allocate();
  bool isAllocated() { return data !=  NULL; } ///< Test if data is allocated

  /// Return a pointer to the time-derivative field
  Field2D* timeDeriv();

  // Operators

  Field2D & operator=(const Field2D &rhs);
  Field2D & operator=(const BoutReal rhs);

  // Data indexing
  BoutReal* operator[](int jx) const;
  BoutReal& operator()(int jx, int jy);
  const BoutReal& operator()(int jx, int jy) const;

  Field2D & operator+=(const Field2D &rhs);
  Field2D & operator+=(const BoutReal rhs);
  Field2D & operator-=(const Field2D &rhs);
  Field2D & operator-=(const BoutReal rhs);
  Field2D & operator*=(const Field2D &rhs);
  Field2D & operator*=(const BoutReal rhs);
  Field2D & operator/=(const Field2D &rhs);
  Field2D & operator/=(const BoutReal rhs);
  Field2D & operator^=(const Field2D &rhs);
  Field2D & operator^=(const BoutReal rhs);
  
  // Binary operators

  const Field2D operator+(const Field2D &other) const;
  const Field2D operator+(const BoutReal rhs) const;
  const Field2D operator-() const;
  const Field2D operator-(const Field2D &other) const;
  const Field2D operator-(const BoutReal rhs) const;
  const Field2D operator*(const Field2D &other) const;
  const Field2D operator*(const BoutReal rhs) const;
  const Field2D operator/(const Field2D &other) const;
  const Field2D operator/(const BoutReal rhs) const;
  const Field2D operator^(const Field2D &other) const;
  const Field2D operator^(const BoutReal rhs) const;

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
  friend const Field2D operator/(const BoutReal lhs, const Field2D &rhs);
  const FieldPerp operator^(const FieldPerp &other) const;
  friend const Field2D operator^(const BoutReal lhs, const Field2D &rhs);
    
  // Stencils

  void getXArray(int y, int z, rvec &xv) const;
  void getYArray(int x, int z, rvec &yv) const;
  void getZArray(int x, int y, rvec &zv) const;

  void setXArray(int y, int z, const rvec &xv);
  void setYArray(int x, int z, const rvec &yv);

  void setStencil(bstencil *fval, bindex *bx) const;
  void setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;

  // Functions
  
  const Field2D sqrt() const;
  const Field2D abs() const;
  BoutReal min(bool allpe=false) const;
  BoutReal max(bool allpe=false) const;
  bool finite() const;
  
  friend const Field2D exp(const Field2D &f);
  friend const Field2D log(const Field2D &f);
  
  friend const Field2D sin(const Field2D &f);
  friend const Field2D cos(const Field2D &f);
  friend const Field2D tan(const Field2D &f);

  friend const Field2D sinh(const Field2D &f);
  friend const Field2D cosh(const Field2D &f);
  friend const Field2D tanh(const Field2D &f);

  bool is_const;
  BoutReal value;

  // FieldData virtual functions
  
  bool isReal() const   { return true; }         // Consists of BoutReal values
  bool is3D() const     { return false; }        // Field is 2D
  int  byteSize() const { return sizeof(BoutReal); } // Just one BoutReal
  int  BoutRealSize() const { return 1; }
  int  getData(int x, int y, int z, void *vptr) const;
  int  getData(int x, int y, int z, BoutReal *rptr) const;
  int  setData(int x, int y, int z, void *vptr);
  int  setData(int x, int y, int z, BoutReal *rptr);

  bool ioSupport() { return true; } ///< This class supports I/O operations
  BoutReal *getData(int component) { 
    BoutReal **d = getData();
    return *d;
  }
  void zeroComponent(int component){
    *this = 0.0;
  }

#ifdef CHECK
  bool checkData(bool vital = true) const; ///< Checks if the data is all valid.
  void doneComms() {bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true; }
#endif
  
  friend class Vector2D;
  
  void applyBoundary();
  void applyBoundary(const string &condition);
  void applyBoundary(const string &region, const string &condition);
  void applyTDerivBoundary();
  void setBoundaryTo(const Field2D &f2d); ///< Copy the boundary region
  
 private:
  BoutReal **data;

  // Data stack: Blocks of memory for this class
  static stack<BoutReal**> block;
  static bool recycle; ///< Re-use blocks rather than freeing/allocating

  void allocData();
  void freeData();
  
  Field2D *deriv; ///< Time-derivative, can be NULL
};

// Non-member overloaded operators

const Field2D operator+(const BoutReal lhs, const Field2D &rhs);
const Field2D operator-(const BoutReal lhs, const Field2D &rhs);
const Field2D operator*(const BoutReal lhs, const Field2D &rhs);
const Field2D operator/(const BoutReal lhs, const Field2D &rhs);
const Field2D operator^(const BoutReal lhs, const Field2D &rhs);

// Non-member functions
const Field2D SQ(const Field2D &f);
const Field2D sqrt(const Field2D &f);
const Field2D abs(const Field2D &f);
BoutReal min(const Field2D &f, bool allpe=false);
BoutReal max(const Field2D &f, bool allpe=false);
bool finite(const Field2D &f);

const Field2D copy(const Field2D &f);

const Field2D floor(const Field2D &var, BoutReal f);

#endif /* __FIELD2D_H__ */
