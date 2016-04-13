/**************************************************************************
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

class Field3D;

#pragma once
#ifndef __FIELD3D_H__
#define __FIELD3D_H__

class Mesh;  // #include "bout/mesh.hxx"
#include "field.hxx"
#include "field2d.hxx"
#include "fieldperp.hxx"
#include "stencils.hxx"
#include "bout_types.hxx"

#include "bout/dataiterator.hxx"

#include "bout/array.hxx"

#include "bout/deprecated.hxx"

#include "bout/field_visitor.hxx"

/// Class for 3D X-Y-Z scalar fields
/*!
  Recycles memory using a global stack. Overloaded operators
  provide operations on data.

  July 2008: Added FieldData virtual functions
  May 2008: Added reference counting to reduce memory copying
 */
class Field3D : public Field, public FieldData {
 public:
  /*!
   * Constructor
   *
   * Note: the global "mesh" can't be passed here because
   * fields may be created before the mesh is.
   */
  Field3D(Mesh *msh = nullptr);
  
  /*! 
   * Copy constructor
   */
  Field3D(const Field3D& f);
  
  /// Constructor from 2D field
  Field3D(const Field2D& f);
  /// Constructor from value
  Field3D(const BoutReal val);
  /// Destructor
  ~Field3D();

  /// Data type
  using value_type = BoutReal;

  /*!
   * Ensures that memory is allocated and unique
   */
  void allocate();
  
  /*!
   * Test if data is allocated
   */
  bool isAllocated() const { return !data.empty(); } 
  
  /*!
   * Return a pointer to the time-derivative field
   */
  Field3D* timeDeriv();

  /*!
   * Ensure that this field has separate fields
   * for yup and ydown.
   */
  void splitYupYdown();

  /*!
   * Ensure that yup and ydown refer to this field
   */
  void mergeYupYdown();
  
  /// Flux Coordinate Independent (FCI) method
  Field3D& yup() { return *yup_field; }
  const Field3D& yup() const { return *yup_field; }
  
  Field3D& ydown() { return *ydown_field; }
  const Field3D& ydown() const { return *ydown_field; }

  /// Return yup if dir=+1, and ydown if dir=-1
  Field3D& ynext(int dir);
  const Field3D& ynext(int dir) const;

  // Staggered grids
  void setLocation(CELL_LOC loc); // Set variable location
  CELL_LOC getLocation() const; // Variable location
  
  /////////////////////////////////////////////////////////
  // Data access
  
  const DataIterator iterator() const;

  const DataIterator begin() const;
  const DataIterator end() const;
  
  /*
   * Returns a range of indices which can be iterated over
   * Uses the REGION flags in bout_types.hxx
   */
  const IndexRange region(REGION rgn) const;

  BoutReal& operator[](DataIterator &d) {
    return operator()(d.x, d.y, d.z);
  }
  const BoutReal& operator[](DataIterator &d) const {
    return operator()(d.x, d.y, d.z);
  }
  BoutReal& operator[](const Indices &i) {
    return operator()(i.x, i.y, i.z);
  }
  const BoutReal& operator[](const Indices &i) const {
    return operator()(i.x, i.y, i.z);
  }
  
  BoutReal& operator[](bindex &bx) {
    return operator()(bx.jx, bx.jy, bx.jz);
  }
  const BoutReal& operator[](bindex &bx) const {
    return operator()(bx.jx, bx.jy, bx.jz);
  }
  
  inline BoutReal& operator()(int jx, int jy, int jz) {
#if CHECK > 2
    // Perform bounds checking
    if(data.empty())
      throw BoutException("Field3D: () operator on empty data");
    
    if((jx < 0) || (jx >= nx) || 
       (jy < 0) || (jy >= ny) || 
       (jz < 0) || (jz >= nz))
      throw BoutException("Field3D: (%d, %d, %d) operator out of bounds (%d, %d, %d)", 
			  jx, jy, jz, nx, ny, nz);
#endif
    return data[(jx*ny +jy)*nz + jz];
  }
  
  inline const BoutReal& operator()(int jx, int jy, int jz) const {
#if CHECK > 2
    if(data.empty())
      throw BoutException("Field3D: () operator on empty data");
    
    if((jx < 0) || (jx >= nx) || 
       (jy < 0) || (jy >= ny) || 
       (jz < 0) || (jz >= nz))
      throw BoutException("Field3D: (%d, %d, %d) operator out of bounds (%d, %d, %d)", 
			  jx, jy, jz, nx, ny, nz);
#endif
    return data[(jx*ny +jy)*nz + jz];
  }

  inline const BoutReal* operator()(int jx, int jy) const {
#if CHECK > 2
    if(data.empty())
      throw BoutException("Field3D: () operator on empty data");

    if((jx < 0) || (jx >= nx) ||
       (jy < 0) || (jy >= ny))
      throw BoutException("Field3D: (%d, %d) operator out of bounds (%d, %d)",
                          jx, jy, nx, ny);
#endif
    return &data[(jx*ny +jy)*nz];
  }

  inline BoutReal* operator()(int jx, int jy) {
#if CHECK > 2
    if(data.empty())
      throw BoutException("Field3D: () operator on empty data");

    if((jx < 0) || (jx >= nx) ||
       (jy < 0) || (jy >= ny))
      throw BoutException("Field3D: (%d, %d) operator out of bounds (%d, %d)",
                          jx, jy, nx, ny);
#endif
    return &data[(jx*ny +jy)*nz];
  }
  
  /////////////////////////////////////////////////////////
  // Operators
  
  const Field3D operator+() {return *this;}
  
  /// Assignment operators
  Field3D & operator=(const Field3D &rhs);
  Field3D & operator=(const Field2D &rhs);
  Field3D & operator=(const FieldPerp &rhs);
  const bvalue & operator=(const bvalue &val);
  BoutReal operator=(const BoutReal val);

  /// Addition operators
  Field3D & operator+=(const Field3D &rhs);
  Field3D & operator+=(const Field2D &rhs);
  Field3D & operator+=(const BoutReal &rhs);
  
  /// Subtraction
  Field3D & operator-=(const Field3D &rhs);
  Field3D & operator-=(const Field2D &rhs);
  Field3D & operator-=(const BoutReal &rhs);

  /// Multiplication
  Field3D & operator*=(const Field3D &rhs);
  Field3D & operator*=(const Field2D &rhs);
  Field3D & operator*=(const BoutReal &rhs);
  
  /// Division
  Field3D & operator/=(const Field3D &rhs);
  Field3D & operator/=(const Field2D &rhs);
  Field3D & operator/=(const BoutReal &rhs);

  // Stencils for differencing
  
  void setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setXStencil(forward_stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setXStencil(backward_stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(forward_stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(backward_stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  
  // FieldData virtual functions
  
  bool isReal() const   { return true; }         // Consists of BoutReal values
  bool is3D() const     { return true; }         // Field is 3D
  int  byteSize() const { return sizeof(BoutReal); } // Just one BoutReal
  int  BoutRealSize() const { return 1; }
  int  getData(int x, int y, int z, void *vptr) const;
  int  getData(int x, int y, int z, BoutReal *rptr) const;
  int  setData(int x, int y, int z, void *vptr);
  int  setData(int x, int y, int z, BoutReal *rptr);

  /// Visitor pattern support
  void accept(FieldVisitor &v) override { v.accept(*this); }
  
#ifdef CHECK
  void doneComms() { bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true; }
#else
  void doneComms() {}
#endif

  friend class Vector3D;

  void setBackground(const Field2D &f2d); // Boundary is applied to the total of this and f2d
  void applyBoundary(bool init=false);
  void applyBoundary(BoutReal t);
  void applyBoundary(const string &condition);
  void applyBoundary(const char* condition) { applyBoundary(string(condition)); }
  void applyBoundary(const string &region, const string &condition);
  void applyTDerivBoundary();
  void setBoundaryTo(const Field3D &f3d); ///< Copy the boundary region

  void applyParallelBoundary();
  void applyParallelBoundary(BoutReal t);
  void applyParallelBoundary(const string &condition);
  void applyParallelBoundary(const char* condition) { applyParallelBoundary(string(condition)); }
  void applyParallelBoundary(const string &region, const string &condition);
  void applyParallelBoundary(const string &region, const string &condition, Field3D *f);
  
private:
  /// Boundary - add a 2D field
  const Field2D *background;

  Mesh *fieldmesh; ///< The mesh over which the field is defined
  int nx, ny, nz;  ///< Array sizes (from fieldmesh). These are valid only if fieldmesh is not null
  
  /// Internal data array. Handles allocation/freeing of memory
  Array<BoutReal> data;
  
  CELL_LOC location; // Location of the variable in the cell
  
  Field3D *deriv; ///< Time derivative (may be NULL)

  /// Pointers to fields containing values along Y
  Field3D *yup_field, *ydown_field;
};

// Non-member overloaded operators

// Binary operators
const FieldPerp operator+(const Field3D &lhs, const FieldPerp &rhs);
const FieldPerp operator-(const Field3D &lhs, const FieldPerp &rhs);
const FieldPerp operator*(const Field3D &lhs, const FieldPerp &rhs);
const FieldPerp operator/(const Field3D &lhs, const FieldPerp &rhs);

const Field3D operator+(const Field3D &lhs, const Field3D &rhs);
const Field3D operator-(const Field3D &lhs, const Field3D &rhs);
const Field3D operator*(const Field3D &lhs, const Field3D &rhs);
const Field3D operator/(const Field3D &lhs, const Field3D &rhs);

const Field3D operator+(const Field3D &lhs, const Field2D &rhs);
const Field3D operator-(const Field3D &lhs, const Field2D &rhs);
const Field3D operator*(const Field3D &lhs, const Field2D &rhs);
const Field3D operator/(const Field3D &lhs, const Field2D &rhs);

const Field3D operator+(const Field3D &lhs, BoutReal rhs);
const Field3D operator-(const Field3D &lhs, BoutReal rhs);
const Field3D operator*(const Field3D &lhs, BoutReal rhs);
const Field3D operator/(const Field3D &lhs, BoutReal rhs);

const Field3D operator+(BoutReal lhs, const Field3D &rhs);
const Field3D operator-(BoutReal lhs, const Field3D &rhs);
const Field3D operator*(BoutReal lhs, const Field3D &rhs);
const Field3D operator/(BoutReal lhs, const Field3D &rhs);

// Unary operators
const Field3D operator-(const Field3D &f);

// Non-member functions
BoutReal min(const Field3D &f, bool allpe=false);
BoutReal max(const Field3D &f, bool allpe=false);

Field3D pow(const Field3D &lhs, const Field3D &rhs);
Field3D pow(const Field3D &lhs, const Field2D &rhs);
Field3D pow(const Field3D &lhs, const FieldPerp &rhs);
Field3D pow(const Field3D &f, BoutReal rhs);
Field3D pow(BoutReal lhs, const Field3D &rhs);

const Field3D SQ(const Field3D &f);
const Field3D sqrt(const Field3D &f);
const Field3D abs(const Field3D &f);

const Field3D exp(const Field3D &f);
const Field3D log(const Field3D &f);

const Field3D sin(const Field3D &f);
const Field3D cos(const Field3D &f);
const Field3D tan(const Field3D &f);

const Field3D sinh(const Field3D &f);
const Field3D cosh(const Field3D &f);
const Field3D tanh(const Field3D &f);

bool finite(const Field3D &var);

void checkData(const Field3D &f); ///< Checks if the data is valid. 
const Field3D copy(const Field3D &f);

const Field3D floor(const Field3D &var, BoutReal f);

const Field3D filter(const Field3D &var, int N0);
const Field3D lowPass(const Field3D &var, int zmax);
const Field3D lowPass(const Field3D &var, int zmax, int zmin);

Field2D DC(const Field3D &f);

/*!
 * @brief Returns a reference to the time-derivative of a field
 * 
 * Wrapper around member function f.timeDeriv()
 *
 */
inline Field3D& ddt(Field3D &f) {
  return *(f.timeDeriv());
}

#endif /* __FIELD3D_H__ */
