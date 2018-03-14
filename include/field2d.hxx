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

#include "bout/field_visitor.hxx"

#include "bout/array.hxx"
#include "bout/region.hxx"

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
   * @param[in] localmesh  The mesh which defines the field size. 
   * 
   * By default the global Mesh pointer (mesh) is used.
   */ 
  Field2D(Mesh *localmesh = nullptr);

  /*!
   * Copy constructor. After this both fields
   * will share the same underlying data.
   */
  Field2D(const Field2D& f);

  /*!
   * Move constructor
   */
  Field2D(Field2D&& f) = default;

  /*!
   * Constructor. This creates a Field2D using the global Mesh pointer (mesh)
   * allocates data, and assigns the value \p val to all points including
   * boundary cells.
   */ 
  Field2D(BoutReal val, Mesh *localmesh = nullptr);

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

  /*!
   * Return the number of nx points
   */
  int getNx() const override {return nx;};
  /*!
   * Return the number of ny points
   */
  int getNy() const override {return ny;};
  /*!
   * Return the number of nz points
   */
  int getNz() const override {return 1;};

  // Operators

  /*!
   * Assignment from Field2D. After this both fields will
   * share the same underlying data. To make a true copy,
   * call .allocate() after assignment, or use the copy()
   * function.
   */
  Field2D & operator=(const Field2D &rhs);
  Field2D & operator=(Field2D &&rhs) = default;

  /*!
   * Allocates data if not already allocated, then
   * sets all cells to \p rhs
   */ 
  Field2D & operator=(BoutReal rhs);

  /// Set variable location for staggered grids to @param new_location
  ///
  /// Throws BoutException if new_location is not `CELL_CENTRE` and
  /// staggered grids are turned off and checks are on. If checks are
  /// off, silently sets location to ``CELL_CENTRE`` instead.
  void setLocation(CELL_LOC new_location) override;
  /// Get variable location
  CELL_LOC getLocation() const override;

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

  BoutReal& operator[](const Ind2D &d) {
    return data[d.ind];
  }
  const BoutReal& operator[](const Ind2D &d) const {
    return data[d.ind];
  }
  BoutReal& operator[](const Ind3D &d); 
  const BoutReal& operator[](const Ind3D &d) const;

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
  inline BoutReal& operator[](const Indices &i) override {
    return operator()(i.x, i.y);
  }
  /// const Indices data access
  inline const BoutReal& operator[](const Indices &i) const override {
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

  // FieldData virtual functions

  /// Visitor pattern support
  void accept(FieldVisitor &v) override {v.accept(*this);}
  
  bool isReal() const override  { return true; }         // Consists of BoutReal values
  bool is3D() const override    { return false; }        // Field is 2D
  int  byteSize() const override { return sizeof(BoutReal); } // Just one BoutReal
  int  BoutRealSize() const override { return 1; }

#if CHECK > 0
  void doneComms() override { bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true; }
#else
  void doneComms() override {}
#endif

  friend class Vector2D;
  
  void applyBoundary(bool init=false) override;
  void applyBoundary(const string &condition);
  void applyBoundary(const char* condition) { applyBoundary(string(condition)); }
  void applyBoundary(const string &region, const string &condition);
  void applyTDerivBoundary() override;
  void setBoundaryTo(const Field2D &f2d); ///< Copy the boundary region
  
 private:
  int nx, ny;      ///< Array sizes (from fieldmesh). These are valid only if fieldmesh is not null
  
  /// Internal data array. Handles allocation/freeing of memory
  Array<BoutReal> data;
  
  CELL_LOC location; ///< Location of the variable in the cell

  Field2D *deriv; ///< Time-derivative, can be NULL
};

// Non-member overloaded operators

Field2D operator+(const Field2D &lhs, const Field2D &rhs);
Field2D operator-(const Field2D &lhs, const Field2D &rhs);
Field2D operator*(const Field2D &lhs, const Field2D &rhs);
Field2D operator/(const Field2D &lhs, const Field2D &rhs);

Field3D operator+(const Field2D &lhs, const Field3D &rhs);
Field3D operator-(const Field2D &lhs, const Field3D &rhs);
Field3D operator*(const Field2D &lhs, const Field3D &rhs);
Field3D operator/(const Field2D &lhs, const Field3D &rhs);

Field2D operator+(const Field2D &lhs, BoutReal rhs);
Field2D operator-(const Field2D &lhs, BoutReal rhs);
Field2D operator*(const Field2D &lhs, BoutReal rhs);
Field2D operator/(const Field2D &lhs, BoutReal rhs);

Field2D operator+(BoutReal lhs, const Field2D &rhs);
Field2D operator-(BoutReal lhs, const Field2D &rhs);
Field2D operator*(BoutReal lhs, const Field2D &rhs);
Field2D operator/(BoutReal lhs, const Field2D &rhs);

/*!
 * Unary minus. Returns the negative of given field,
 * iterates over whole domain including guard/boundary cells.
 */
Field2D operator-(const Field2D &f);

// Non-member functions

/// Square root
const Field2D sqrt(const Field2D &f, REGION rgn=RGN_ALL);

/// Absolute value
const Field2D abs(const Field2D &f, REGION rgn=RGN_ALL);

/*!
 * Calculates the minimum of a field, excluding
 * the boundary/guard cells by default (this can be
 * changed with the rgn argument).
 * By default this is only on the local processor,
 * but setting allpe=true does a collective Allreduce
 * over all processors.
 *
 * @param[in] f  The field to loop over
 * @param[in] allpe  Minimum over all processors?
 * @param[in] rgn  the boundaries that should be ignored
 * 
 */
BoutReal min(const Field2D &f, bool allpe=false, REGION rgn=RGN_NOBNDRY);

/*!
 * Calculates the maximum of a field, excluding
 * the boundary/guard cells by default (this can be
 * changed with the rgn argument).
 * By default this is only on the local processor,
 * but setting allpe=true does a collective Allreduce
 * over all processors.
 *
 * @param[in] f  The field to loop over
 * @param[in] allpe  Minimum over all processors?
 * @param[in] rgn  the boundaries that should be ignored
 * 
 */
BoutReal max(const Field2D &f, bool allpe=false, REGION rgn=RGN_NOBNDRY);

/*!
 * Test if all values of this field are finite
 * Loops over the entire domain including boundaries by
 * default (can be changed using the rgn argument)
 */
bool finite(const Field2D &f, REGION rgn=RGN_ALL);

/// Exponential
const Field2D exp(const Field2D &f, REGION rgn=RGN_ALL);

/// Natural logarithm
const Field2D log(const Field2D &f, REGION rgn=RGN_ALL);

/*!
 * Sine trigonometric function. 
 *
 * @param[in] f  Angle in radians
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D sin(const Field2D &f, REGION rgn=RGN_ALL);

/*!
 * Cosine trigonometric function. 
 *
 * @param[in] f  Angle in radians
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D cos(const Field2D &f, REGION rgn=RGN_ALL);

/*!
 * Tangent trigonometric function. 
 *
 * @param[in] f  Angle in radians
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D tan(const Field2D &f, REGION rgn=RGN_ALL);

/*!
 * Hyperbolic sine function. 
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D sinh(const Field2D &f, REGION rgn=RGN_ALL);

/*!
 * Hyperbolic cosine function. 
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D cosh(const Field2D &f, REGION rgn=RGN_ALL);

/*!
 * Hyperbolic tangent function. 
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field2D tanh(const Field2D &f, REGION rgn=RGN_ALL);

/// Make an independent copy of field f
const Field2D copy(const Field2D &f);

/// Sets a floor on var, so minimum of the return value is >= f
const Field2D floor(const Field2D &var, BoutReal f, REGION rgn=RGN_ALL);

/// Power, lhs ** rhs
Field2D pow(const Field2D &lhs, const Field2D &rhs, REGION rgn=RGN_ALL);
Field2D pow(const Field2D &lhs, BoutReal rhs, REGION rgn=RGN_ALL);
Field2D pow(BoutReal lhs, const Field2D &rhs, REGION rgn=RGN_ALL);

#if CHECK > 0
void checkData(const Field2D &f, REGION region = RGN_NOBNDRY);
#else
inline void checkData(const Field2D &UNUSED(f), REGION UNUSED(region) = RGN_NOBNDRY) {}
#endif

/*!
 * Force guard cells of passed field to nan
 */ 
void invalidateGuards(Field2D &var);

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
