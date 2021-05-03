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

#include "bout/array.hxx"
#include "bout/region.hxx"

#include "bout/assert.hxx"

#include "bout/field_visitor.hxx"

#include "utils.hxx"

#include <vector>


/// Class for 3D X-Y-Z scalar fields
/*!
  This class represents a scalar field defined over the mesh.
  It handles memory management, and provides overloaded operators
  for operations on the data, iterators and access methods.

  Initialisation
  --------------

  Fields can be declared in any scope (even global),
  but cannot be accessed by index or used until the data
  is allocated.

      Field3D f;   // Declare variable, no data allocated
      f(0,0,0) = 1.0; // Error !

      f = 0.0;  // Allocates memory, fills with value (0.0)

      Field3D g(1.0); // Declares, allocates memory, fills with value (1.0)

      Field3D h;   // not allocated
      h.allocate();  // Data array allocated, values undefined
      f(0,0,0) = 1.0; // ok

  Copy-on-Write
  -------------

  A field is a reference to the underlying data array, so
  setting one field equal to another has the effect of making
  both fields share the same underlying data

      Field3D f(0.0);
      Field3D g = f; // f and g now share data
      f(0,0,0) = 1.0; // g is also modified

  Setting the entire field equal to a new value changes the reference:

      Field3D f(0.0);
      Field3D g = f; // f and g now share data
      g = 1.0;   // g and f are now separate

  To ensure that a field is unique, call allocate() which
  will make a copy of the underlying data if it is shared.

      Field3D f(0.0);
      Field3D g = f; // f and g now share data
      g.allocate();  // Data copied so g and f don't share data
      f(0,0,0) = 1.0; // ok

  Data access
  -----------

  Individual data indices can be accessed by index using
  round brackets:

      Field3D f;
      f(0,1,2) = 1.0;  // Set value
      BoutReal val = f(2,1,3);  // Get value

  If CHECK is greater than 2, this function will perform
  bounds checking. This will significantly slow calculations.

  Some methods, such as FFT routines, need access to
  a pointer to memory. For the Z dimension this can be done
  by passing only the X and Y indices

      BoutReal *data = f(0,1);

  `data` now points to `f(0,1,0)` and can be incremented to move in Z.

  Iteration
  ---------

  To loop over all points in a field, a for loop can be used
  to get the indices:

      Field3D f(0.0); // Allocate, set to zero

      for( const auto &i : f ) {  // Loop over all points, with index i
        f[i] = 1.0;
      }

  There is also more explicit looping over regions:

      for( const auto &i : f.region(RGN_ALL) ) {  // Loop over all points, with index i
        f[i] = 1.0;
      }

  Parallel (y) derivatives
  ------------------------

  In several numerical schemes the mapping along magnetic fields
  (default y direction) is a relatively complex map. To accommodate
  this, the values of a field in the positive (up) and negative (down)
  directions can be stored in separate fields.

      Field3D f(0.0); // f allocated, set to zero

      f.yup() // error; f.yup not allocated

      f.clearParallelSlices(); // f.yup_fields and f.ydown_fields are now empty
      f.yup() // error; f.yup not allocated

      To have separate fields for yup and ydown, first call

      f.splitParallelSlices(); // f.yup() and f.ydown() separate

      f.yup(); // ok
      f.yup()(0,1,0) // error; f.yup not allocated

      f.yup() = 1.0; // Set f.yup() field to 1.0

      f.yup()(0,1,0) // ok

 */
class Field3D : public Field, public FieldData {
 public:
  using ind_type = Ind3D;
  
  /*!
   * Constructor
   *
   * Note: the global "mesh" can't be passed here because
   * fields may be created before the mesh is.
   */
  Field3D(Mesh *localmesh = nullptr, CELL_LOC location_in=CELL_CENTRE,
          DirectionTypes directions_in =
            {YDirectionType::Standard, ZDirectionType::Standard});

  /*!
   * Copy constructor
   */
  Field3D(const Field3D& f);

  /// Move constructor
  Field3D(Field3D&& f) noexcept { swap(*this, f); }
  /// Constructor from 2D field
  Field3D(const Field2D& f);
  /// Constructor from value
  Field3D(BoutReal val, Mesh *localmesh = nullptr);
  /// Constructor from Array and Mesh
  Field3D(Array<BoutReal> data, Mesh* localmesh, CELL_LOC location = CELL_CENTRE,
          DirectionTypes directions_in = {YDirectionType::Standard,
                                          ZDirectionType::Standard});
  /// Destructor
  ~Field3D() override;

  /// Data type stored in this field
  using value_type = BoutReal;

  /*!
   * Ensures that memory is allocated and unique
   */
  Field3D& allocate();
  
  /*!
   * Test if data is allocated
   */
  bool isAllocated() const { return !data.empty(); } 
  
  /*!
   * Return a pointer to the time-derivative field
   *
   * The first time this is called, a new field will be
   * allocated. Subsequent calls return the same field
   */
  Field3D* timeDeriv();

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
  int getNz() const override {return nz;};

  // these methods return Field3D to allow method chaining
  Field3D& setLocation(CELL_LOC new_location) {
    Field::setLocation(new_location);
    return *this;
  }
  Field3D& setDirectionY(YDirectionType d) {
    directions.y = d;
    return *this;
  }

  /*!
   * Ensure that this field has separate fields
   * for yup and ydown.
   */
  void splitParallelSlices();

  [[deprecated("Please use Field3D::splitParallelSlices instead")]]
  void splitYupYdown() {
    splitParallelSlices();
  }

  /*!
   * Clear the parallel slices, yup and ydown
   */
  void clearParallelSlices();
  
  [[deprecated("Please use Field3D::clearParallelSlices instead")]]
  void mergeYupYdown() {
    clearParallelSlices();
  }

  /// Check if this field has yup and ydown fields
  bool hasParallelSlices() const {
    return !yup_fields.empty() and !ydown_fields.empty();
  }

  [[deprecated("Please use Field3D::hasParallelSlices instead")]]
  bool hasYupYdown() const {
    return hasParallelSlices();
  }

  /// Check if this field has yup and ydown fields
  /// Return reference to yup field
  Field3D &yup(std::vector<Field3D>::size_type index = 0) {
    ASSERT2(index < yup_fields.size());
    return yup_fields[index];
  }
  /// Return const reference to yup field
  const Field3D &yup(std::vector<Field3D>::size_type index = 0) const {
    ASSERT2(index < yup_fields.size());
    return yup_fields[index];
  }

  /// Return reference to ydown field
  Field3D &ydown(std::vector<Field3D>::size_type index = 0) {
    ASSERT2(index < ydown_fields.size());
    return ydown_fields[index];
  }

  /// Return const reference to ydown field
  const Field3D &ydown(std::vector<Field3D>::size_type index = 0) const {
    ASSERT2(index < ydown_fields.size());
    return ydown_fields[index];
  }

  /// Return the parallel slice at \p offset
  ///
  /// \p offset of 0 returns the main field itself
  Field3D& ynext(int offset);
  const Field3D& ynext(int offset) const;

  /// If \p twist_shift_enabled is true, does this Field3D require a twist-shift at branch
  /// cuts on closed field lines?
  bool requiresTwistShift(bool twist_shift_enabled);

  /////////////////////////////////////////////////////////
  // Data access

  /// Return a Region<Ind3D> reference to use to iterate over this field
  ///
  /// Example
  /// -------
  /// 
  /// This loops over the interior points, not the boundary
  /// and inside the loop the index is used to calculate the difference
  /// between the point one index up in x (i.xp()) and one index down
  /// in x (i.xm()), putting the result into a different field 'g'
  /// 
  /// for(const auto &i : f.getRegion(RGN_NOBNDRY)) {
  ///   g[i] = f[i.xp()] - f[i.xm()];
  /// }
  /// 
  const Region<Ind3D>& getRegion(REGION region) const;  
  const Region<Ind3D>& getRegion(const std::string &region_name) const;
  /// Use region provided by the default, and if none is set, use the provided one
  const Region<Ind3D>& getDefaultRegion(const std::string& region_name) const;
  void setRegion(const std::string& region_name);
  void resetRegion() { regionID = -1; };
  void setRegion(int id) { regionID = id; };
  int getRegionID() const { return regionID; };

  /// Return a Region<Ind2D> reference to use to iterate over the x- and
  /// y-indices of this field
  const Region<Ind2D>& getRegion2D(REGION region) const;
  const Region<Ind2D>& getRegion2D(const std::string &region_name) const;
  
  Region<Ind3D>::RegionIndices::const_iterator begin() const {return std::begin(getRegion("RGN_ALL"));};
  Region<Ind3D>::RegionIndices::const_iterator end() const {return std::end(getRegion("RGN_ALL"));};
  
  BoutReal& operator[](const Ind3D &d) {
    return data[d.ind];
  }
  const BoutReal& operator[](const Ind3D &d) const {
    return data[d.ind];
  }

  BoutReal& operator()(const IndPerp &d, int jy);
  const BoutReal& operator()(const IndPerp &d, int jy) const;

  BoutReal& operator()(const Ind2D &d, int jz);
  const BoutReal& operator()(const Ind2D &d, int jz) const;
  
  /*!
   * Direct access to the underlying data array
   *
   * If CHECK > 2 then bounds checking is performed
   * 
   * If CHECK <= 2 then no checks are performed, to
   * allow inlining and optimisation of inner loops
   */
  inline BoutReal& operator()(int jx, int jy, int jz) {
#if CHECK > 2
    // Perform bounds checking
    if(data.empty())
      throw BoutException("Field3D: () operator on empty data");
    
    if((jx < 0) || (jx >= nx) || 
       (jy < 0) || (jy >= ny) || 
       (jz < 0) || (jz >= nz))
      throw BoutException(
          "Field3D: ({:d}, {:d}, {:d}) operator out of bounds ({:d}, {:d}, {:d})", jx, jy,
          jz, nx, ny, nz);
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
      throw BoutException(
          "Field3D: ({:d}, {:d}, {:d}) operator out of bounds ({:d}, {:d}, {:d})", jx, jy,
          jz, nx, ny, nz);
#endif
    return data[(jx*ny +jy)*nz + jz];
  }
  
  /*!
   * Direct access to the underlying data array
   *
   * This version returns a pointer to a data array,
   * and is intended for use with FFT routines. The data
   * is guaranteed to be contiguous in Z index
   */
  inline const BoutReal* operator()(int jx, int jy) const {
#if CHECK > 2
    if(data.empty())
      throw BoutException("Field3D: () operator on empty data");

    if((jx < 0) || (jx >= nx) ||
       (jy < 0) || (jy >= ny))
      throw BoutException("Field3D: ({:d}, {:d}) operator out of bounds ({:d}, {:d})", jx,
                          jy, nx, ny);
#endif
    return &data[(jx*ny +jy)*nz];
  }

  inline BoutReal* operator()(int jx, int jy) {
#if CHECK > 2
    if(data.empty())
      throw BoutException("Field3D: () operator on empty data");

    if((jx < 0) || (jx >= nx) ||
       (jy < 0) || (jy >= ny))
      throw BoutException("Field3D: ({:d}, {:d}) operator out of bounds ({:d}, {:d})", jx,
                          jy, nx, ny);
#endif
    return &data[(jx*ny +jy)*nz];
  }
  
  /////////////////////////////////////////////////////////
  // Operators
  
  /// Assignment operators
  ///@{
  Field3D & operator=(const Field3D &rhs);
  Field3D & operator=(Field3D&& rhs);
  Field3D & operator=(const Field2D &rhs);
  /// return void, as only part initialised
  void      operator=(const FieldPerp &rhs);
  Field3D & operator=(BoutReal val);
  ///@}

  /// Addition operators
  ///@{
  Field3D & operator+=(const Field3D &rhs);
  Field3D & operator+=(const Field2D &rhs);
  Field3D & operator+=(BoutReal rhs);
  ///@}
  
  /// Subtraction operators
  ///@{
  Field3D & operator-=(const Field3D &rhs);
  Field3D & operator-=(const Field2D &rhs);
  Field3D & operator-=(BoutReal rhs);
  ///@}

  /// Multiplication operators
  ///@{
  Field3D & operator*=(const Field3D &rhs);
  Field3D & operator*=(const Field2D &rhs);
  Field3D & operator*=(BoutReal rhs);
  ///@}

  /// Division operators
  ///@{
  Field3D & operator/=(const Field3D &rhs);
  Field3D & operator/=(const Field2D &rhs);
  Field3D & operator/=(BoutReal rhs);
  ///@}
  
  // FieldData virtual functions
  
  bool isReal() const override   { return true; }         // Consists of BoutReal values
  bool is3D() const override     { return true; }         // Field is 3D
  int  byteSize() const override { return sizeof(BoutReal); } // Just one BoutReal
  int  BoutRealSize() const override { return 1; }

  /// Visitor pattern support
  void accept(FieldVisitor &v) override { v.accept(*this); }
  
#if CHECK > 0
  void doneComms() override { bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true; }
#else
  void doneComms() override {}
#endif

  friend class Vector3D;
  friend class Vector2D;

  Field3D& calcParallelSlices();

  void applyBoundary(bool init=false) override;
  void applyBoundary(BoutReal t);
  void applyBoundary(const std::string &condition);
  void applyBoundary(const char* condition) { applyBoundary(std::string(condition)); }
  void applyBoundary(const std::string &region, const std::string &condition);
  void applyTDerivBoundary() override;
  
  /// Copy the boundary values half-way between cells
  /// This uses 2nd order central differences to set the value
  /// on the boundary to the value on the boundary in field \p f3d.
  /// Note: does not just copy values in boundary region.
  void setBoundaryTo(const Field3D &f3d); 

  void applyParallelBoundary();
  void applyParallelBoundary(BoutReal t);
  void applyParallelBoundary(const std::string &condition);
  void applyParallelBoundary(const char* condition) { applyParallelBoundary(std::string(condition)); }
  void applyParallelBoundary(const std::string &region, const std::string &condition);
  void applyParallelBoundary(const std::string &region, const std::string &condition, Field3D *f);

  friend void swap(Field3D& first, Field3D& second) noexcept {
    using std::swap;

    // Swap base class members
    swap(static_cast<Field&>(first), static_cast<Field&>(second));

    swap(first.data, second.data);
    swap(first.background, second.background);
    swap(first.nx, second.nx);
    swap(first.ny, second.ny);
    swap(first.nz, second.nz);
    swap(first.deriv, second.deriv);
    swap(first.yup_fields, second.yup_fields);
    swap(first.ydown_fields, second.ydown_fields);
    swap(first.bndry_op, second.bndry_op);
    swap(first.boundaryIsCopy, second.boundaryIsCopy);
    swap(first.boundaryIsSet, second.boundaryIsSet);
    swap(first.bndry_op_par, second.bndry_op_par);
    swap(first.bndry_generator, second.bndry_generator);
  }
  
private:
  /// Boundary - add a 2D field
  const Field2D *background{nullptr};

  /// Array sizes (from fieldmesh). These are valid only if fieldmesh is not null
  int nx{-1}, ny{-1}, nz{-1};

  /// Internal data array. Handles allocation/freeing of memory
  Array<BoutReal> data;
  
  /// Time derivative (may be nullptr)
  Field3D *deriv{nullptr};

  /// Fields containing values along Y
  std::vector<Field3D> yup_fields{}, ydown_fields{};

  /// RegionID over which the field is valid
  int regionID{-1};
};

// Non-member overloaded operators

// Binary operators
FieldPerp operator+(const Field3D &lhs, const FieldPerp &rhs);
FieldPerp operator-(const Field3D &lhs, const FieldPerp &rhs);
FieldPerp operator*(const Field3D &lhs, const FieldPerp &rhs);
FieldPerp operator/(const Field3D &lhs, const FieldPerp &rhs);

Field3D operator+(const Field3D &lhs, const Field3D &rhs);
Field3D operator-(const Field3D &lhs, const Field3D &rhs);
Field3D operator*(const Field3D &lhs, const Field3D &rhs);
Field3D operator/(const Field3D &lhs, const Field3D &rhs);

Field3D operator+(const Field3D &lhs, const Field2D &rhs);
Field3D operator-(const Field3D &lhs, const Field2D &rhs);
Field3D operator*(const Field3D &lhs, const Field2D &rhs);
Field3D operator/(const Field3D &lhs, const Field2D &rhs);

Field3D operator+(const Field3D &lhs, BoutReal rhs);
Field3D operator-(const Field3D &lhs, BoutReal rhs);
Field3D operator*(const Field3D &lhs, BoutReal rhs);
Field3D operator/(const Field3D &lhs, BoutReal rhs);

Field3D operator+(BoutReal lhs, const Field3D &rhs);
Field3D operator-(BoutReal lhs, const Field3D &rhs);
Field3D operator*(BoutReal lhs, const Field3D &rhs);
Field3D operator/(BoutReal lhs, const Field3D &rhs);

/*!
 * Unary minus. Returns the negative of given field,
 * iterates over whole domain including guard/boundary cells.
 */
Field3D operator-(const Field3D &f);

// Non-member functions

/// Exponent: pow(lhs, lhs) is \p lhs raised to the power of \p rhs
///
/// Extra overloads not provided by the templates in field.hxx
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
Field3D pow(const Field3D& lhs, const Field2D& rhs, const std::string& rgn = "RGN_ALL");
[[deprecated("Please use pow(const Field3D& lhs, const Field2D& rhs"
    "const std::string& region = \"RGN_ALL\") instead")]]
inline Field3D pow(const Field3D &lhs, const Field2D &rhs, REGION rgn) {
  return pow(lhs, rhs, toString(rgn));
}
FieldPerp pow(const Field3D& lhs, const FieldPerp& rhs, const std::string& rgn = "RGN_ALL");
[[deprecated("Please use pow(const Field3D& lhs, const FieldPerp& rhs"
    "const std::string& region = \"RGN_ALL\") instead")]]
inline FieldPerp pow(const Field3D& lhs, const FieldPerp& rhs, REGION rgn) {
  return pow(lhs, rhs, toString(rgn));
}

#if CHECK > 0
/// Throw an exception if \p f is not allocated or if any
/// elements are non-finite (for CHECK > 2).
/// Loops over all points including the boundaries by
/// default (can be changed using the \p rgn argument
void checkData(const Field3D& f, const std::string& region = "RGN_NOBNDRY");
[[deprecated("Please use checkData(const Field3D& f, "
    "const std::string& region = \"RGN_NOBNDRY\") instead")]]
inline void checkData(const Field3D &f, REGION region) {
  return checkData(f, toString(region));
}
#else
/// Ignored with disabled CHECK; Throw an exception if \p f is not
/// allocated or if any elements are non-finite (for CHECK > 2)
inline void checkData(const Field3D& UNUSED(f), const std::string& UNUSED(region) = "RGN_NOBNDRY") {};
[[deprecated("Please use checkData(const Field3D& f, "
    "const std::string& region = \"RGN_NOBNDRY\") instead")]]
inline void checkData(const Field3D &UNUSED(f), REGION UNUSED(region)) {}
#endif

/// Fourier filtering, removes all except one mode
///
/// @param[in] var Variable to apply filter to
/// @param[in] N0  The component to keep
/// @param[in] rgn The region to calculate the result over
Field3D filter(const Field3D& var, int N0, const std::string& rgn = "RGN_ALL");
[[deprecated("Please use filter(const Field3D& var, int N0, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline Field3D filter(const Field3D& var, int N0, REGION rgn) {
  return filter(var, N0, toString(rgn));
}

/// Fourier low pass filtering. Removes modes
/// higher than \p zmax and optionally the zonal component
///
/// @param[in] var   Variable to apply filter to
/// @param[in] zmax  Maximum mode in Z
/// @param[in] keep_zonal  Keep the zonal component if true
/// @param[in] rgn   The region to calculate the result over
Field3D lowPass(const Field3D& var, int zmax, bool keep_zonal,
    const std::string& rgn = "RGN_ALL");
[[deprecated("Please use lowpass(const Field3D& var, int zmax, bool keep_zonal, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline Field3D lowPass(const Field3D& var, int zmax, bool keep_zonal, REGION rgn) {
  return lowPass(var, zmax, keep_zonal, toString(rgn));
}

/// The argument \p keep_zonal used to be integer "zmin" -- this was a
/// misnomer. Please use the version above which uses a bool instead
DEPRECATED(inline Field3D lowPass(const Field3D& var, int zmax, int keep_zonal,
                                  REGION rgn = RGN_ALL)) {
  ASSERT0(static_cast<bool>(keep_zonal) == keep_zonal);
  return lowPass(var, zmax, static_cast<bool>(keep_zonal), toString(rgn));
}

/// Fourier low pass filtering. Removes modes higher than \p zmax
///
/// @param[in] var   Variable to apply filter to
/// @param[in] zmax  Maximum mode in Z
/// @param[in] rgn   The region to calculate the result over
inline Field3D lowPass(const Field3D &var, int zmax, const std::string rgn = "RGN_ALL") {
  return lowPass(var, zmax, true, rgn);
}
[[deprecated("Please use lowpass(const Field3D& var, int zmax, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline Field3D lowPass(const Field3D &var, int zmax, REGION rgn) {
  return lowPass(var, zmax, toString(rgn));
}

/// Perform a shift by a given angle in Z
///
/// @param[inout] var  The variable to be modified in-place
/// @param[in] jx      X index
/// @param[in] jy      Y index
/// @param[in] zangle  The Z angle to apply
void shiftZ(Field3D &var, int jx, int jy, double zangle);

/// Apply a phase shift by a given angle \p zangle in Z to all points
///
/// @param[inout] var  The variable to modify in-place
/// @param[in] zangle  The angle to shift by in Z
/// @param[in] rgn     The region to calculate the result over
void shiftZ(Field3D &var, BoutReal zangle, const std::string& rgn="RGN_ALL");
[[deprecated("Please use shiftZ(const Field3D& var, BoutReal zangle, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline void shiftZ(Field3D &var, BoutReal zangle, REGION rgn) {
  return shiftZ(var, zangle, toString(rgn));
}

/// Average in the Z direction
///
/// @param[in] f     Variable to average
/// @param[in] rgn   The region to calculate the result over
Field2D DC(const Field3D &f, const std::string& rgn = "RGN_ALL");
[[deprecated("Please use DC(const Field3D& f, "
    "const std::string& region = \"RGN_ALL\") instead")]]
inline Field2D DC(const Field3D &f, REGION rgn) {
  return DC(f, toString(rgn));
}

/// Force guard cells of passed field \p var to NaN
#if CHECK > 2
void invalidateGuards(Field3D &var);
#else
inline void invalidateGuards(Field3D &UNUSED(var)) {}
#endif

/// Returns a reference to the time-derivative of a field \p f
///
/// Wrapper around member function f.timeDeriv()
inline Field3D& ddt(Field3D &f) {
  return *(f.timeDeriv());
}

/// toString template specialisation
/// Defined in utils.hxx
template <>
inline std::string toString<>(const Field3D& UNUSED(val)) {
  return "<Field3D>";
}

/// Test if two fields are the same, by calculating
/// the minimum absolute difference between them
bool operator==(const Field3D &a, const Field3D &b);

/// Output a string describing a Field3D to a stream
std::ostream& operator<<(std::ostream &out, const Field3D &value);

#endif /* __FIELD3D_H__ */
