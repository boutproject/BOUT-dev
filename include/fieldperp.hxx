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
  ~FieldPerp() {}

  /*!
   * Assignment operators
   */
  FieldPerp &operator=(const FieldPerp &rhs);
  FieldPerp &operator=(FieldPerp &&rhs) = default;
  FieldPerp &operator=(BoutReal rhs);

  /*!
   * Iterators and data access
   */
  const DataIterator begin() const;
  const DataIterator end() const;

  const IndexRange region(REGION rgn) const;

  /*!
   * Direct data access using DataIterator indexing
   */
  inline BoutReal& operator[](const DataIterator &d) {
    return operator()(d.x, d.z);
  }
  inline const BoutReal& operator[](const DataIterator &d) const {
    return operator()(d.x, d.z);
  }
  BoutReal& operator[](const Indices &i) {
    return operator()(i.x, i.z);
  }
  const BoutReal& operator[](const Indices &i) const override{
    return operator()(i.x, i.z);
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
  void allocate() {
    if(data.empty()) {
      data = Array<BoutReal>(nx*nz);
    }else
      data.ensureUnique();
  }

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

  virtual int getNy() const override{ return 1;};
  
 private:
  int yindex; ///< The Y index at which this FieldPerp is defined

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
 * Create a unique copy of a FieldPerp, ensuring 
 * that they do not share an underlying data array
 */
const FieldPerp copy(const FieldPerp &f);

/*!
 * Create a FieldPerp by slicing a 3D field at a given y
 */
const FieldPerp sliceXZ(const Field3D& f, int y);

#endif
