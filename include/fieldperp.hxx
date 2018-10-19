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

#include "bout/dataiterator.hxx"
#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/region.hxx"

#include "unused.hxx"

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
  FieldPerp(Mesh * fieldmesh = nullptr);

  /*!
   * Copy constructor. After this the data
   * will be shared (non unique)
   */
  FieldPerp(const FieldPerp &f)
      : Field(f.fieldmesh), yindex(f.yindex), nx(f.nx), nz(f.nz), data(f.data) {}

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

  ~FieldPerp() override {}

  /*!
   * Assignment operators
   */
  FieldPerp &operator=(const FieldPerp &rhs);
  FieldPerp &operator=(FieldPerp &&rhs) = default;
  FieldPerp &operator=(BoutReal rhs);

  /*!
   * Iterators and data access
   */
  const DataIterator DEPRECATED(begin()) const;
  const DataIterator DEPRECATED(end()) const;

  const IndexRange DEPRECATED(region(REGION rgn)) const override;

  /// Return a Region<IndPerp> reference to use to iterate over this field
  const Region<IndPerp>& getRegion(REGION region) const;  
  const Region<IndPerp>& getRegion(const std::string &region_name) const;

  /*!
   * Direct data access using DataIterator indexing
   */
  inline BoutReal& DEPRECATED(operator[](const DataIterator &d)) {
    return operator()(d.x, d.z);
  }
  inline const BoutReal& DEPRECATED(operator[](const DataIterator &d)) const {
    return operator()(d.x, d.z);
  }
  BoutReal& DEPRECATED(operator[](const Indices &i)) {
    return operator()(i.x, i.z);
  }
  const BoutReal& DEPRECATED(operator[](const Indices &i)) const override{
    return operator()(i.x, i.z);
  }

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

  /*!
   * Returns the y index at which this field is defined
   */ 
  int getIndex() const {return yindex;}
  
  /*!
   * Sets the y index at which this field is defined
   *
   * This is used in arithmetic operations
   */
  void setIndex(int y) { yindex = y; }

  /*!
   * Ensure that data array is allocated and unique
   */
  void allocate();

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
      throw BoutException("FieldPerp: (%d, %d) operator out of bounds (%d, %d)", 
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
      throw BoutException("FieldPerp: (%d, %d) operator out of bounds (%d, %d)", 
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
  
 private:
  int yindex = -1; ///< The Y index at which this FieldPerp is defined

  /// The size of the data array
  int nx, nz;

  /// The underlying data array
  Array<BoutReal> data;
};
  
// Non-member overloaded operators
  
const FieldPerp operator+(const FieldPerp &lhs, const FieldPerp &rhs);
const FieldPerp operator+(const FieldPerp &lhs, const Field3D &rhs);
const FieldPerp operator+(const FieldPerp &lhs, const Field2D &rhs);
const FieldPerp operator+(const FieldPerp &lhs, BoutReal rhs);
inline const FieldPerp operator+(BoutReal lhs, const FieldPerp &rhs) {
  return rhs + lhs;
}

const FieldPerp operator-(const FieldPerp &lhs, const FieldPerp &other);
const FieldPerp operator-(const FieldPerp &lhs, const Field3D &other);
const FieldPerp operator-(const FieldPerp &lhs, const Field2D &other);
const FieldPerp operator-(const FieldPerp &lhs, BoutReal rhs);
const FieldPerp operator-(BoutReal lhs, const FieldPerp &rhs);

const FieldPerp operator*(const FieldPerp &lhs, const FieldPerp &other);
const FieldPerp operator*(const FieldPerp &lhs, const Field3D &other);
const FieldPerp operator*(const FieldPerp &lhs, const Field2D &other);
const FieldPerp operator*(const FieldPerp &lhs, BoutReal rhs);
inline const FieldPerp operator*(BoutReal lhs, const FieldPerp &rhs) {
  return rhs * lhs;
}

const FieldPerp operator/(const FieldPerp &lhs, const FieldPerp &other);
const FieldPerp operator/(const FieldPerp &lhs, const Field3D &other);
const FieldPerp operator/(const FieldPerp &lhs, const Field2D &other);
const FieldPerp operator/(const FieldPerp &lhs, BoutReal rhs);
const FieldPerp operator/(BoutReal lhs, const FieldPerp &rhs);

/*!
 * Unary minus. Returns the negative of given field,
 * iterates over whole domain including guard/boundary cells.
 */
FieldPerp operator-(const FieldPerp &f);

/// Unary plus. This does nothing, but is sometimes useful
/// when commenting out code or for clarification.
FieldPerp operator+(const FieldPerp &f) { return f; }

/// Square root
const FieldPerp sqrt(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Absolute value
const FieldPerp abs(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Exponential
const FieldPerp exp(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Natural logarithm
const FieldPerp log(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Sine trigonometric function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
const FieldPerp sin(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Cosine trigonometric function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
const FieldPerp cos(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Tangent trigonometric function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
const FieldPerp tan(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Hyperbolic sine function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
const FieldPerp sinh(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Hyperbolic cosine function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
const FieldPerp cosh(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Hyperbolic tangent function.
///
/// @param[in] f    Angle in radians
/// @param[in] rgn  The region to calculate the result over
///
/// This loops over the entire domain, including guard/boundary cells by
/// default (can be changed using the \p rgn argument).
/// If CHECK >= 3 then the result will be checked for non-finite numbers
const FieldPerp tanh(const FieldPerp &f, REGION rgn=RGN_ALL);

/// Create a unique copy of a FieldPerp, ensuring
/// that they do not share an underlying data array
const FieldPerp copy(const FieldPerp &f);

/// Sets a floor on var, so minimum of the return value is >= f
const FieldPerp floor(const FieldPerp &var, BoutReal f, REGION rgn=RGN_ALL);

/// Power, lhs ** rhs
FieldPerp pow(const FieldPerp &lhs, const FieldPerp &rhs, REGION rgn=RGN_ALL);
FieldPerp pow(const FieldPerp &lhs, BoutReal rhs, REGION rgn=RGN_ALL);
FieldPerp pow(BoutReal lhs, const FieldPerp &rhs, REGION rgn=RGN_ALL);

/// Create a FieldPerp by slicing a 3D field at a given y
const FieldPerp sliceXZ(const Field3D& f, int y);

/// Calculates the minimum of a field, excluding
/// the boundary/guard cells by default (this can be
/// changed with the rgn argument).
/// By default this is only on the local processor,
/// but setting allpe=true does a collective Allreduce
/// over all processors.
///
/// @param[in] f      The field to loop over
/// @param[in] allpe  Minimum over all processors?
/// @param[in] rgn    The region to calculate the result over
BoutReal min(const FieldPerp &f, bool allpe=false, REGION rgn=RGN_NOX);

/// Calculates the maximum of a field, excluding
/// the boundary/guard cells by default (this can be
/// changed with the rgn argument).
/// By default this is only on the local processor,
/// but setting allpe=true does a collective Allreduce
/// over all processors.
///
/// @param[in] f      The field to loop over
/// @param[in] allpe  Minimum over all processors?
/// @param[in] rgn    The region to calculate the result over
BoutReal max(const FieldPerp &f, bool allpe=false, REGION rgn=RGN_NOX);

/// Test if all values of this field are finite
/// Loops over the entire domain including boundaries by
/// default (can be changed using the \p rgn argument)
bool finite(const FieldPerp &f, REGION rgn=RGN_ALL);

#if CHECK > 0
void checkData(const FieldPerp &f, REGION region = RGN_NOX);
#else
inline void checkData(const FieldPerp &UNUSED(f), REGION UNUSED(region) = RGN_NOX) {}
#endif

/// Force guard cells of passed field \p var to NaN
#if CHECK > 2
void invalidateGuards(FieldPerp &var);
#else
inline void invalidateGuards(FieldPerp &UNUSED(var)) {}
#endif

#endif
