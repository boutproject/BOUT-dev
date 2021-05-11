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

#include "bout/field_visitor.hxx"

#include "bout/array.hxx"
#include "bout/region.hxx"
#include "utils.hxx"

#include "unused.hxx"

/*!
 * \brief 2D X-Y scalar fields
 *
 * Handles data for axisymmetric quantities. Essentially the same
 * as the Field3D class.
 */
class Field2D : public Field {
 public:
  using ind_type = Ind2D;    
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
  Field2D(Mesh *localmesh = nullptr, CELL_LOC location_in=CELL_CENTRE,
          DirectionTypes directions_in =
          {YDirectionType::Standard, ZDirectionType::Average});

  /*!
   * Copy constructor. After this both fields
   * will share the same underlying data.
   */
  Field2D(const Field2D& f);

  /*!
   * Move constructor
   */
  Field2D(Field2D&& f) noexcept { swap(*this, f); };

  /*!
   * Constructor. This creates a Field2D using the global Mesh pointer (mesh)
   * allocates data, and assigns the value \p val to all points including
   * boundary cells.
   */ 
  Field2D(BoutReal val, Mesh *localmesh = nullptr);
  
  /// Constructor from Array and Mesh
  Field2D(Array<BoutReal> data, Mesh* localmesh, CELL_LOC location = CELL_CENTRE,
          DirectionTypes directions_in = {YDirectionType::Standard,
                                          ZDirectionType::Average});

  /*!
   * Destructor
   */
  ~Field2D() override;

  /// Data type
  using value_type = BoutReal;

  /// Ensure data is allocated
  Field2D& allocate();
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

  // these methods return Field2D to allow method chaining
  Field2D& setLocation(CELL_LOC new_location) {
    Field::setLocation(new_location);
    return *this;
  }
  Field2D& setDirectionY(YDirectionType d) {
    // This method included in case it is wanted in a templated function also dealing with
    // Field3D or FieldPerp - there is no difference between orthogonal and field-aligned
    // coordinates for Field2D, so should always have YDirectionType::Standard.
    ASSERT1(d == YDirectionType::Standard);
    directions.y = d;
    return *this;
  }

  /// Check if this field has yup and ydown fields
  bool hasParallelSlices() const {
    return true;
  }

  [[deprecated("Please use Field2D::hasParallelSlices instead")]]
  bool hasYupYdown() const {
    return hasParallelSlices();
  }
  
  Field2D& yup(std::vector<Field2D>::size_type UNUSED(index) = 0) {
    return *this;
  }
  const Field2D& yup(std::vector<Field2D>::size_type UNUSED(index) = 0) const {
    return *this;
  }

  Field2D& ydown(std::vector<Field2D>::size_type UNUSED(index) = 0) {
    return *this;
  }
  const Field2D& ydown(std::vector<Field2D>::size_type UNUSED(index) = 0) const {
    return *this;
  }

  Field2D& ynext(int UNUSED(dir)) { return *this; }
  const Field2D& ynext(int UNUSED(dir)) const { return *this; }

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

  /// Return a Region<Ind2D> reference to use to iterate over this field
  const Region<Ind2D>& getRegion(REGION region) const;  
  const Region<Ind2D>& getRegion(const std::string &region_name) const;

  Region<Ind2D>::RegionIndices::const_iterator begin() const {return std::begin(getRegion("RGN_ALL"));};
  Region<Ind2D>::RegionIndices::const_iterator end() const {return std::end(getRegion("RGN_ALL"));};
  
  BoutReal& operator[](const Ind2D &d) {
    return data[d.ind];
  }
  const BoutReal& operator[](const Ind2D &d) const {
    return data[d.ind];
  }
  BoutReal& operator[](const Ind3D &d);
  const BoutReal& operator[](const Ind3D &d) const;

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
      throw BoutException("Field2D: ({:d}, {:d}) index out of bounds ({:d} , {:d})\n", jx,
                          jy, nx, ny);
#endif
  
    return data[jx*ny + jy];
  }
  inline const BoutReal& operator()(int jx, int jy) const {
#if CHECK > 2
    if(!isAllocated())
      throw BoutException("Field2D: () operator on empty data");
    
    if((jx < 0) || (jx >= nx) || 
       (jy < 0) || (jy >= ny) )
      throw BoutException("Field2D: ({:d}, {:d}) index out of bounds ({:d} , {:d})\n", jx,
                          jy, nx, ny);
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

  bool is3D() const override { return false; }

#if CHECK > 0
  void doneComms() override { bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true; }
#else
  void doneComms() override {}
#endif

  friend class Vector2D;
  
  void applyBoundary(bool init=false) override;
  void applyBoundary(BoutReal time);
  void applyBoundary(const std::string &condition);
  void applyBoundary(const char* condition) { applyBoundary(std::string(condition)); }
  void applyBoundary(const std::string &region, const std::string &condition);
  void applyTDerivBoundary() override;
  void setBoundaryTo(const Field2D &f2d); ///< Copy the boundary region

  friend void swap(Field2D& first, Field2D& second) noexcept {
    using std::swap;

    // Swap base class members
    swap(static_cast<Field&>(first), static_cast<Field&>(second));

    swap(first.data, second.data);
    swap(first.nx, second.nx);
    swap(first.ny, second.ny);
    swap(first.deriv, second.deriv);
    swap(first.bndry_op, second.bndry_op);
    swap(first.boundaryIsCopy, second.boundaryIsCopy);
    swap(first.boundaryIsSet, second.boundaryIsSet);
    swap(first.bndry_op_par, second.bndry_op_par);
    swap(first.bndry_generator, second.bndry_generator);
  }

private:
  /// Array sizes (from fieldmesh). These are valid only if fieldmesh is not null
  int nx{-1}, ny{-1};

  /// Internal data array. Handles allocation/freeing of memory
  Array<BoutReal> data;

  /// Time-derivative, can be nullptr
  Field2D *deriv{nullptr};
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

inline Field2D toFieldAligned(const Field2D& f, const std::string& UNUSED(region) = "RGN_ALL") {
  return f;
}
[[deprecated("Please use toFieldAligned(const Field2D& f, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline Field2D toFieldAligned(const Field2D& f, REGION region) {
  return toFieldAligned(f, toString(region));
}

inline Field2D fromFieldAligned(const Field2D& f, const std::string& UNUSED(region) = "RGN_ALL") {
  return f;
}
[[deprecated("Please use fromFieldAligned(const Field2D& f, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline Field2D fromFieldAligned(const Field2D& f, REGION region) {
  return fromFieldAligned(f, toString(region));
}

#if CHECK > 0
/// Throw an exception if \p f is not allocated or if any
/// elements are non-finite (for CHECK > 2).
/// Loops over all points including the boundaries by
/// default (can be changed using the \p rgn argument
void checkData(const Field2D &f, const std::string& region = "RGN_NOBNDRY");
[[deprecated("Please use checkData(const Field2D& f, "
    "const std::string& region = \"RGN_NOBNDRY\") instead")]]
inline void checkData(const Field2D &f, REGION region) {
  return checkData(f, toString(region));
}
#else
inline void checkData(const Field2D &UNUSED(f), std::string UNUSED(region) = "RGN_NOBNDRY") {}
[[deprecated("Please use checkData(const Field2D& f, "
    "const std::string& region = \"RGN_NOBNDRY\") instead")]]
inline void checkData(const Field2D &UNUSED(f), REGION UNUSED(region)) {}
#endif

/// Force guard cells of passed field \p var to NaN
#if CHECK > 2
void invalidateGuards(Field2D &var);
#else
inline void invalidateGuards(Field2D &UNUSED(var)) {}
#endif

/// Returns a reference to the time-derivative of a field \p f
///
/// Wrapper around member function f.timeDeriv()
inline Field2D& ddt(Field2D &f) {
  return *(f.timeDeriv());
}

/// toString template specialisation
/// Defined in utils.hxx
template <>
inline std::string toString<>(const Field2D& UNUSED(val)) {
  return "<Field2D>";
}

/// Test if two fields are the same, by calculating
/// the minimum absolute difference between them
bool operator==(const Field2D &a, const Field2D &b);

std::ostream& operator<<(std::ostream &out, const Field2D &value);

#endif /* __FIELD2D_H__ */
