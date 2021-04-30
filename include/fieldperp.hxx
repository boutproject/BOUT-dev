/**************************************************************************
 * Class for 2D X-Z slices
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
 **************************************************************************/

class FieldPerp;

#ifndef __FIELDPERP_H__
#define __FIELDPERP_H__

#include "field.hxx"

#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/region.hxx"

#include "unused.hxx"

#include <ostream>
#include <string>

class Field2D; // #include "field2d.hxx"
class Field3D; // #include "field3d.hxx"

/*!
 * Represents a 2D field perpendicular to the magnetic field
 * at a particular index in Y, which only varies in X-Z. 
 * 
 * Primarily used inside field solvers
 */ 
class FieldPerp : public Field {
 public:
  using ind_type = IndPerp;
    
  /*!
   * Constructor
   */
  FieldPerp(Mesh * fieldmesh = nullptr, CELL_LOC location_in=CELL_CENTRE,
            int yindex_in=-1,
            DirectionTypes directions_in =
              {YDirectionType::Standard, ZDirectionType::Standard});

  /*!
   * Copy constructor. After this the data
   * will be shared (non unique)
   */
  FieldPerp(const FieldPerp& f) = default;

  /*!
   * Move constructor
   */
  FieldPerp(FieldPerp &&rhs) = default;

  /*!
   * Constructor. This creates a FieldPerp using the global Mesh pointer (mesh)
   * allocates data, and assigns the value \p val to all points including
   * boundary cells.
   */ 
  FieldPerp(BoutReal val, Mesh *localmesh = nullptr);

  /*!
   * Constructor from Array and Mesh
   */
  FieldPerp(Array<BoutReal> data, Mesh* fieldmesh, CELL_LOC location_in = CELL_CENTRE,
            int yindex_in = -1,
            DirectionTypes directions_in = {YDirectionType::Standard,
                                            ZDirectionType::Standard});

  ~FieldPerp() override = default;

  /*!
   * Assignment operators
   */
  FieldPerp &operator=(const FieldPerp &rhs);
  FieldPerp &operator=(FieldPerp &&rhs) = default;
  FieldPerp &operator=(BoutReal rhs);

  /// Return a Region<IndPerp> reference to use to iterate over this field
  const Region<IndPerp>& getRegion(REGION region) const;  
  const Region<IndPerp>& getRegion(const std::string &region_name) const;
  const Region<IndPerp>& getDefaultRegion(const std::string& region_name) const {
    return getRegion(region_name);
  }

  Region<IndPerp>::RegionIndices::const_iterator begin() const {return std::begin(getRegion("RGN_ALL"));};
  Region<IndPerp>::RegionIndices::const_iterator end() const {return std::end(getRegion("RGN_ALL"));};
  
  inline BoutReal& operator[](const IndPerp &d) {
    return data[d.ind];
  }
  inline const BoutReal& operator[](const IndPerp &d) const {
    return data[d.ind];
  }  

  inline BoutReal& operator[](const Ind3D &d) {
    ASSERT3(d.y() == yindex);
    return operator()(d.x(), d.z()); //Could use mesh->ind3DtoPerp if we had access to mesh here
  }
  inline const BoutReal& operator[](const Ind3D &d) const {
    ASSERT3(d.y() == yindex);
    return operator()(d.x(), d.z());
  }

  /// Return the y index at which this field is defined. This value is
  /// local to each processor
  int getIndex() const { return yindex; }

  /// Return the globally defined y index if it's either an interior
  /// (grid) point, or a boundary point. Otherwise, return -1 to
  /// indicate a guard cell or an invalid value
  int getGlobalIndex() const;

  /// Set the (local) y index at which this field is defined
  ///
  /// This is used in arithmetic operations
  FieldPerp& setIndex(int y) {
    yindex = y;
    return *this;
  }

  /// Set the (local) y index at which this field is defined from a
  /// globally defined y index
  ///
  /// Only use the global y index if it's either an interior (grid)
  /// point, or a boundary point. Otherwise, sets yindex to -1 to
  /// indicate a guard cell or an invalid value
  FieldPerp& setIndexFromGlobal(int y_global);

  // these methods return FieldPerp to allow method chaining
  FieldPerp& setLocation(CELL_LOC new_location) override {
    Field::setLocation(new_location);
    return *this;
  }
  FieldPerp& setDirectionY(YDirectionType d) override {
    Field::setDirectionY(d);
    return *this;
  }

  /*!
   * Ensure that data array is allocated and unique
   */
  FieldPerp& allocate();

  /*!
   * True if the underlying data array is allocated.
   *
   *
   */
  bool isAllocated() const { return !data.empty(); }
  
  // operators
  
  const BoutReal* operator[](int jx) const {
    ASSERT2(!data.empty());
    ASSERT2( (jx >= 0) && (jx < nx) );
  
    return &data[jx*nz];
  }
  
  /*!
   * Returns a C-style array (pointer to first element) in Z
   * at a given X index. Used mainly for FFT routines
   */
  BoutReal* operator[](int jx) {
    ASSERT2(!data.empty());
    ASSERT2( (jx >= 0) && (jx < nx) );
    
    return &data[jx*nz];
  }

  /*!
   * Access to the underlying data array at a given x,z index
   * 
   * If CHECK > 2 then bounds checking is performed, otherwise
   * no checks are performed
   */ 
  BoutReal& operator()(int jx, int jz) {
#if CHECK > 2
    // Bounds check both indices
    if(data.empty())
      throw BoutException("FieldPerp: () operator on empty data");
    if((jx < 0) || (jx >= nx) || 
       (jz < 0) || (jz >= nz))
      throw BoutException("FieldPerp: ({:d}, {:d}) operator out of bounds ({:d}, {:d})",
                          jx, jz, nx, nz);
#endif
    return data[jx*nz + jz];
  }
  
  /*!
   * Const (read-only) access to the underlying data array.
   */ 
  const BoutReal& operator()(int jx, int jz) const {
#if CHECK > 2
    // Bounds check both indices
    if(data.empty())
      throw BoutException("FieldPerp: () operator on empty data");
    if((jx < 0) || (jx >= nx) || 
       (jz < 0) || (jz >= nz))
      throw BoutException("FieldPerp: ({:d}, {:d}) operator out of bounds ({:d}, {:d})",
                          jx, jz, nx, nz);
#endif
    return data[jx*nz + jz];
  }
  
  /*!
   * Access to the underlying data array. (X,Y,Z) indices for consistency with 
   * other field types
   * 
   */ 
  BoutReal& operator()(int jx, int UNUSED(jy), int jz) { return (*this)(jx, jz); }
  
  const BoutReal& operator()(int jx, int UNUSED(jy), int jz) const { return (*this)(jx, jz); }

  /*!
   * Addition, modifying in-place. 
   * This loops over the entire domain, including guard/boundary cells
   */
  FieldPerp & operator+=(const FieldPerp &rhs);
  FieldPerp & operator+=(const Field3D &rhs);
  FieldPerp & operator+=(const Field2D &rhs);
  FieldPerp & operator+=(BoutReal rhs);

  /*!
   * Subtraction, modifying in place. 
   * This loops over the entire domain, including guard/boundary cells
   */
  FieldPerp & operator-=(const FieldPerp &rhs);
  FieldPerp & operator-=(const Field3D &rhs);
  FieldPerp & operator-=(const Field2D &rhs);
  FieldPerp & operator-=(BoutReal rhs);

  /*!
   * Multiplication, modifying in place. 
   * This loops over the entire domain, including guard/boundary cells
   */
  FieldPerp & operator*=(const FieldPerp &rhs);
  FieldPerp & operator*=(const Field3D &rhs);
  FieldPerp & operator*=(const Field2D &rhs);
  FieldPerp & operator*=(BoutReal rhs);

  /*!
   * Division, modifying in place. 
   * This loops over the entire domain, including guard/boundary cells
   */
  FieldPerp & operator/=(const FieldPerp &rhs);
  FieldPerp & operator/=(const Field3D &rhs);
  FieldPerp & operator/=(const Field2D &rhs);
  FieldPerp & operator/=(BoutReal rhs);

  /*!
   * Return the number of nx points
   */
  int getNx() const override {return nx;};
  /*!
   * Return the number of ny points
   */
  int getNy() const override { return 1; };
  /*!
   * Return the number of nz points
   */
  int getNz() const override {return nz;};

  bool is3D() const override { return false; }

  friend void swap(FieldPerp& first, FieldPerp& second) noexcept;

private:
  /// The Y index at which this FieldPerp is defined
  int yindex{-1};

  /// The size of the data array
  int nx{-1}, nz{-1};

  /// The underlying data array
  Array<BoutReal> data;
};
  
// Non-member functions

// Non-member overloaded operators
  
FieldPerp operator+(const FieldPerp &lhs, const FieldPerp &rhs);
FieldPerp operator+(const FieldPerp &lhs, const Field3D &rhs);
FieldPerp operator+(const FieldPerp &lhs, const Field2D &rhs);
FieldPerp operator+(const FieldPerp &lhs, BoutReal rhs);
FieldPerp operator+(BoutReal lhs, const FieldPerp &rhs);

FieldPerp operator-(const FieldPerp &lhs, const FieldPerp &rhs);
FieldPerp operator-(const FieldPerp &lhs, const Field3D &rhs);
FieldPerp operator-(const FieldPerp &lhs, const Field2D &rhs);
FieldPerp operator-(const FieldPerp &lhs, BoutReal rhs);
FieldPerp operator-(BoutReal lhs, const FieldPerp &rhs);

FieldPerp operator*(const FieldPerp &lhs, const FieldPerp &rhs);
FieldPerp operator*(const FieldPerp &lhs, const Field3D &rhs);
FieldPerp operator*(const FieldPerp &lhs, const Field2D &rhs);
FieldPerp operator*(const FieldPerp &lhs, BoutReal rhs);
FieldPerp operator*(BoutReal lhs, const FieldPerp &rhs);

FieldPerp operator/(const FieldPerp &lhs, const FieldPerp &rhs);
FieldPerp operator/(const FieldPerp &lhs, const Field3D &rhs);
FieldPerp operator/(const FieldPerp &lhs, const Field2D &rhs);
FieldPerp operator/(const FieldPerp &lhs, BoutReal rhs);
FieldPerp operator/(BoutReal lhs, const FieldPerp &rhs);

/*!
 * Unary minus. Returns the negative of given field,
 * iterates over whole domain including guard/boundary cells.
 */
FieldPerp operator-(const FieldPerp &f);

/// Create a FieldPerp by slicing a 3D field at a given y
const FieldPerp sliceXZ(const Field3D& f, int y);

// Specialize newEmptyField templates for FieldPerp
/// Return an empty shell field of some type derived from Field, with metadata
/// copied and a data array that is allocated but not initialised.
template<>
inline FieldPerp emptyFrom<FieldPerp>(const FieldPerp& f) {
  return FieldPerp(f.getMesh(), f.getLocation(), f.getIndex(), {f.getDirectionY(), f.getDirectionZ()}).allocate();
}

#if CHECK > 0
void checkData(const FieldPerp &f, const std::string& region = "RGN_NOX");
[[deprecated("Please use checkData(const FieldPerp& f, "
    "const std::string& region = \"RGN_NOBNDRY\") instead")]]
inline void checkData(const FieldPerp &f, REGION region) {
  return checkData(f, toString(region));
}
#else
inline void checkData(const FieldPerp &UNUSED(f), const std::string& UNUSED(region) = "RGN_NOX") {}
[[deprecated("Please use checkData(const FieldPerp& f, "
    "const std::string& region = \"RGN_NOBNDRY\") instead")]]
inline void checkData(const FieldPerp &UNUSED(f), REGION UNUSED(region)) {}
#endif

/// Force guard cells of passed field \p var to NaN
#if CHECK > 2
void invalidateGuards(FieldPerp &var);
#else
inline void invalidateGuards(FieldPerp &UNUSED(var)) {}
#endif

/// toString template specialisation
/// Defined in utils.hxx
template <>
inline std::string toString<>(const FieldPerp& UNUSED(val)) {
  return "<FieldPerp>";
}

/// Test if two fields are the same, by checking that they are defined
/// at the same y-index, and if the minimum absolute difference
/// between them is less than 1e-10
bool operator==(const FieldPerp& a, const FieldPerp& b);

/// Output a string describing a FieldPerp to a stream
std::ostream& operator<<(std::ostream& out, const FieldPerp& value);

#endif
