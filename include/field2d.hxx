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

class Mesh;
#include "field.hxx"
#include "field_data.hxx"
class Field3D; //#include "field3d.hxx"
#include "fieldperp.hxx"
#include "stencils.hxx"

#include "bout/dataiterator.hxx"

#include "bout/deprecated.hxx"

#include "bout/field_visitor.hxx"

#include "bout/array.hxx"

#include "unused.hxx"

/*!
 * \brief 2D X-Y scalar fields
 *
 * Handles data for axisymmetric quantities. Essentially the same
 * as the Field3D class.
 */
class Field2D : public Field, public FieldData {
 public:
  /*!
   * Constructor, taking an optional mesh pointer
   * This mesh pointer is not used until the data is allocated,
   * since Field2D objects can be globals, created before a mesh
   * has been created.
   *
   * @param[in] msh  The mesh which defines the field size. 
   * 
   * By default the global Mesh pointer (mesh) is used.
   */ 
  Field2D(Mesh *msh = nullptr);

  /*!
   * Copy constructor. After this both fields
   * will share the same underlying data.
   */
  Field2D(const Field2D& f);

  /*!
   * Constructor. This creates a Field2D using the global Mesh pointer (mesh)
   * allocates data, and assigns the value \p val to all points including
   * boundary cells.
   */ 
  Field2D(BoutReal val);

  /*!
   * Destructor
   */ 
  ~Field2D();

  /// Data type
  using value_type = BoutReal;

  /// Ensure data is allocated
  void allocate();
  bool isAllocated() const { return !data.empty(); } ///< Test if data is allocated

  /// Return a pointer to the time-derivative field
  Field2D* timeDeriv();

  // Operators

  /*!
   * Assignment from Field2D. After this both fields will
   * share the same underlying data. To make a true copy,
   * call .allocate() after assignment, or use the copy()
   * function.
   */
  Field2D & operator=(const Field2D &rhs);

  /*!
   * Allocates data if not already allocated, then
   * sets all cells to \p rhs
   */ 
  Field2D & operator=(BoutReal rhs);

  /////////////////////////////////////////////////////////
  // Data access

  /// Iterator over the Field2D indices
  const DataIterator iterator() const;

  const DataIterator begin() const;
  const DataIterator end() const;
  
  /*!
   * Returns a range of indices which can be iterated over
   * Uses the REGION flags in bout_types.hxx
   */
  const IndexRange region(REGION rgn) const;

  /*!
   * Direct access to the data array. Since operator() is used
   * to implement this, no checks are performed if CHECK <= 2
   */
  inline BoutReal& operator[](const DataIterator &d) {
    return operator()(d.x, d.y);
  }

  /// Const access to data array
  inline const BoutReal& operator[](const DataIterator &d) const {
    return operator()(d.x, d.y);
  }

  /// Indices are also used as a lightweight way to specify indexing
  /// for example DataIterator offsets (xp, xm, yp etc.) return Indices
  inline BoutReal& operator[](const Indices &i) {
    return operator()(i.x, i.y);
  }
  /// const Indices data access
  inline const BoutReal& operator[](const Indices &i) const {
    return operator()(i.x, i.y);
  }

  /*!
   * Access to the underlying data array. 
   * 
   * If CHECK <= 2 then no checks are performed
   *
   * If CHECK > 2 then both \p jx and \p jy are bounds checked. This will
   * significantly reduce performance.
   */
  inline BoutReal& operator()(int jx, int jy) {
#if CHECK > 2
    if(!isAllocated())
      throw BoutException("Field2D: () operator on empty data");
    
    if((jx < 0) || (jx >= nx) || 
       (jy < 0) || (jy >= ny) )
      throw BoutException("Field2D: (%d, %d) index out of bounds (%d , %d)\n", 
                          jx, jy, nx, ny);
#endif
  
    return data[jx*ny + jy];
  }
  inline const BoutReal& operator()(int jx, int jy) const {
#if CHECK > 2
    if(!isAllocated())
      throw BoutException("Field2D: () operator on empty data");
    
    if((jx < 0) || (jx >= nx) || 
       (jy < 0) || (jy >= ny) )
      throw BoutException("Field2D: (%d, %d) index out of bounds (%d , %d)\n", 
                          jx, jy, nx, ny);
#endif
  
    return data[jx*ny + jy];
  }

  /*!
   * DIrect access to underlying array. This version is for compatibility
   * with Field3D objects
   */
  BoutReal& operator()(int jx, int jy, int UNUSED(jz)) {
    return operator()(jx, jy);
  }
  const BoutReal& operator()(int jx, int jy, int UNUSED(jz)) const {
    return operator()(jx, jy);
  }
  
  Field2D & operator+=(const Field2D &rhs); ///< In-place addition. Copy-on-write used if data is shared
  Field2D & operator+=(BoutReal rhs);       ///< In-place addition. Copy-on-write used if data is shared
  Field2D & operator-=(const Field2D &rhs); ///< In-place subtraction. Copy-on-write used if data is shared
  Field2D & operator-=(BoutReal rhs);       ///< In-place subtraction. Copy-on-write used if data is shared
  Field2D & operator*=(const Field2D &rhs); ///< In-place multiplication. Copy-on-write used if data is shared
  Field2D & operator*=(BoutReal rhs);       ///< In-place multiplication. Copy-on-write used if data is shared
  Field2D & operator/=(const Field2D &rhs); ///< In-place division. Copy-on-write used if data is shared
  Field2D & operator/=(BoutReal rhs);       ///< In-place division. Copy-on-write used if data is shared

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

  DEPRECATED(int getData(int x, int y, int z, void *vptr) const);
  DEPRECATED(int getData(int x, int y, int z, BoutReal *rptr) const);
  DEPRECATED(int setData(int x, int y, int z, void *vptr));
  DEPRECATED(int setData(int x, int y, int z, BoutReal *rptr));
  
#ifdef CHECK
  void doneComms() { bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true; }
#else
  void doneComms() {}
#endif

  friend class Vector2D;
  
  void applyBoundary(bool init=false);
  void applyBoundary(const string &condition);
  void applyBoundary(const char* condition) { applyBoundary(string(condition)); }
  void applyBoundary(const string &region, const string &condition);
  void applyTDerivBoundary();
  void setBoundaryTo(const Field2D &f2d); ///< Copy the boundary region
  
 private:
  Mesh *fieldmesh; ///< The mesh over which the field is defined
  int nx, ny;      ///< Array sizes (from fieldmesh). These are valid only if fieldmesh is not null
  
  /// Internal data array. Handles allocation/freeing of memory
  Array<BoutReal> data;
  
  Field2D *deriv; ///< Time-derivative, can be NULL
};

// Non-member overloaded operators

const Field2D operator+(const Field2D &lhs, const Field2D &rhs);
const Field2D operator-(const Field2D &lhs, const Field2D &rhs);
const Field2D operator*(const Field2D &lhs, const Field2D &rhs);
const Field2D operator/(const Field2D &lhs, const Field2D &rhs);

const Field3D operator+(const Field2D &lhs, const Field3D &rhs);
const Field3D operator-(const Field2D &lhs, const Field3D &rhs);
const Field3D operator*(const Field2D &lhs, const Field3D &rhs);
const Field3D operator/(const Field2D &lhs, const Field3D &rhs);

const Field2D operator+(const Field2D &lhs, BoutReal rhs);
const Field2D operator-(const Field2D &lhs, BoutReal rhs);
const Field2D operator*(const Field2D &lhs, BoutReal rhs);
const Field2D operator/(const Field2D &lhs, BoutReal rhs);

const Field2D operator+(BoutReal lhs, const Field2D &rhs);
const Field2D operator-(BoutReal lhs, const Field2D &rhs);
const Field2D operator*(BoutReal lhs, const Field2D &rhs);
const Field2D operator/(BoutReal lhs, const Field2D &rhs);

/*!
 * Unary minus. Returns the negative of given field,
 * iterates over whole domain including guard/boundary cells.
 */
const Field2D operator-(const Field2D &f);

// Non-member functions

/// Square root
const Field2D sqrt(const Field2D &f);

/// Absolute value
const Field2D abs(const Field2D &f);

/*!
 * Calculates the minimum of a field, excluding
 * the boundary/guard cells. 
 * By default this is only on the local processor,
 * but setting allpe=true does a collective Allreduce
 * over all processors.
 *
 * @param[in] f  The field to loop over
 * @param[in] allpe  Minimum over all processors?
 * 
 */
BoutReal min(const Field2D &f, bool allpe=false);

/*!
 * Calculates the maximum of a field, excluding
 * the boundary/guard cells. 
 * By default this is only on the local processor,
 * but setting allpe=true does a collective Allreduce
 * over all processors.
 *
 * @param[in] f  The field to loop over
 * @param[in] allpe  Minimum over all processors?
 * 
 */
BoutReal max(const Field2D &f, bool allpe=false);

/*!
 * Test if all values of this field are finite
 * Loops over the entire domain including boundaries
 */
bool finite(const Field2D &f);

/// Exponential
const Field2D exp(const Field2D &f);

/// Natural logarithm
const Field2D log(const Field2D &f);

/*!
 * Sine trigonometric function. 
 *
 * @param[in] f  Angle in radians
 *
 * This loops over the entire domain, including guard/boundary cells
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D sin(const Field2D &f);

/*!
 * Cosine trigonometric function. 
 *
 * @param[in] f  Angle in radians
 *
 * This loops over the entire domain, including guard/boundary cells
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D cos(const Field2D &f);

/*!
 * Tangent trigonometric function. 
 *
 * @param[in] f  Angle in radians
 *
 * This loops over the entire domain, including guard/boundary cells
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D tan(const Field2D &f);

/*!
 * Hyperbolic sine function. 
 *
 * This loops over the entire domain, including guard/boundary cells
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D sinh(const Field2D &f);

/*!
 * Hyperbolic cosine function. 
 *
 * This loops over the entire domain, including guard/boundary cells
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D cosh(const Field2D &f);

/*!
 * Hyperbolic tangent function. 
 *
 * This loops over the entire domain, including guard/boundary cells
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D tanh(const Field2D &f);

/// Make an independent copy of field f
const Field2D copy(const Field2D &f);

/// Sets a floor on var, so minimum of the return value is >= f
const Field2D floor(const Field2D &var, BoutReal f);

/// Power, lhs ** rhs
Field2D pow(const Field2D &lhs, const Field2D &rhs);
Field2D pow(const Field2D &lhs, BoutReal rhs);
Field2D pow(BoutReal lhs, const Field2D &rhs);

#ifdef CHECK
void checkData(const Field2D &f);
#else
inline void checkData(const Field2D &f) {}
#endif


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
