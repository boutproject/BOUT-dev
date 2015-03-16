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

#include "bout/dataiterator.hxx"

#include "bout/deprecated.hxx"

#include "bout/field_visitor.hxx"

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

  /// Data type
  using value_type = BoutReal;

  
  DEPRECATED(Field2D* clone() const);
  
  static void cleanup(); // Frees all memory

  /// Ensure data is allocated
  void allocate();
  bool isAllocated() const { return data !=  NULL; } ///< Test if data is allocated

  /// Return a pointer to the time-derivative field
  Field2D* timeDeriv();

  // Operators

  Field2D & operator=(const Field2D &rhs);
  Field2D & operator=(const BoutReal rhs);

  /////////////////////////////////////////////////////////
  // Data access

  const DataIterator iterator() const;

  const DataIterator begin() const;
  const DataIterator end() const;

  inline BoutReal& operator[](DataIterator &d) {
    return operator()(d.x, d.y);
  }
  inline const BoutReal& operator[](DataIterator &d) const {
    return operator()(d.x, d.y);
  }
  inline BoutReal& operator[](Indices &i) {
    return operator()(i.x, i.y);
  }
  inline const BoutReal& operator[](Indices &i) const {
    return operator()(i.x, i.y);
  }
  
  inline BoutReal& operator()(int jx, int jy) {
#if CHECK > 2
    if(!isAllocated())
      throw BoutException("Field2D: () operator on empty data");
    if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) )
      throw BoutException("Field2D: (%d, %d) index out of bounds (%d , %d)\n", 
                          jx, jy, mesh->ngx, mesh->ngy);
#endif
  
    return data[jx][jy];
  }
  inline const BoutReal& operator()(int jx, int jy) const {
    #if CHECK > 2
    if(!isAllocated())
      throw BoutException("Field2D: () operator on empty data");
    if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) )
      throw BoutException("Field2D: (%d, %d) index out of bounds (%d , %d)\n", 
                          jx, jy, mesh->ngx, mesh->ngy);
#endif
  
    return data[jx][jy];
  }

  Field2D & operator+=(const Field2D &rhs);
  Field2D & operator+=(const BoutReal rhs);
  Field2D & operator-=(const Field2D &rhs);
  Field2D & operator-=(const BoutReal rhs);
  Field2D & operator*=(const Field2D &rhs);
  Field2D & operator*=(const BoutReal rhs);
  Field2D & operator/=(const Field2D &rhs);
  Field2D & operator/=(const BoutReal rhs);
  
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

  // Left binary operators

  const Field3D operator+(const Field3D &other) const;
  const Field3D operator-(const Field3D &other) const;
  const Field3D operator*(const Field3D &other) const;
  const Field3D operator/(const Field3D &other) const;

  const FieldPerp operator+(const FieldPerp &other) const;
  const FieldPerp operator-(const FieldPerp &other) const;
  const FieldPerp operator*(const FieldPerp &other) const;
  const FieldPerp operator/(const FieldPerp &other) const;
  friend const Field2D operator/(const BoutReal lhs, const Field2D &rhs);
    
  // Stencils

  void getXArray(int y, int z, rvec &xv) const;
  void getYArray(int x, int z, rvec &yv) const;
  void getZArray(int x, int y, rvec &zv) const;

  void setXArray(int y, int z, const rvec &xv);
  void setYArray(int x, int z, const rvec &yv);

  //void setStencil(bstencil *fval, bindex *bx) const;
  void setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setXStencil(forward_stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setXStencil(backward_stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(forward_stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(backward_stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  
  // FieldData virtual functions

  /// Visitor pattern support
  void accept(FieldVisitor &v) override {v.accept(*this);}
  
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
    return data[0];
  }
  void zeroComponent(int component){
    *this = 0.0;
  }

#ifdef CHECK
  bool checkData(bool vital = true) const; ///< Checks if the data is all valid.
  void doneComms() { bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true; }
#endif
  
  friend class Vector2D;
  
  void applyBoundary(bool init=false);
  void applyBoundary(const string &condition);
  void applyBoundary(const char* condition) { applyBoundary(string(condition)); }
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


/// Square
const Field2D SQ(const Field2D &f);

/// Square root
const Field2D sqrt(const Field2D &f);

/// Absolute value
const Field2D abs(const Field2D &f);

/// Minimum over field. By default only on this processor
BoutReal min(const Field2D &f, bool allpe=false);

/// Maximum over field. By default only on this processor
BoutReal max(const Field2D &f, bool allpe=false);

/// Test if all values of this field are finite
bool finite(const Field2D &f);

/// Exponential
const Field2D exp(const Field2D &f);

/// Natural logarithm
const Field2D log(const Field2D &f);
  
const Field2D sin(const Field2D &f);
const Field2D cos(const Field2D &f);
const Field2D tan(const Field2D &f);

const Field2D sinh(const Field2D &f);
const Field2D cosh(const Field2D &f);
const Field2D tanh(const Field2D &f);

/// Make an independent copy of field f
const Field2D copy(const Field2D &f);

/// Sets a floor on var, so minimum of the return value is >= f
const Field2D floor(const Field2D &var, BoutReal f);

/// Power, lhs ** rhs
Field2D pow(const Field2D &lhs, const Field2D &rhs);
Field2D pow(const Field2D &lhs, BoutReal rhs);
Field2D pow(BoutReal lhs, const Field2D &rhs);

/*!
 * @brief Returns a reference to the time-derivative of a field
 * 
 * Wrapper around member function f.timeDeriv()
 *
 */
inline Field2D& ddt(Field2D &f) {
  return *(f.timeDeriv());
}

#endif /* __FIELD2D_H__ */
