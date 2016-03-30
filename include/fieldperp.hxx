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

#include "bout/deprecated.hxx"

#include "bout/dataiterator.hxx"
#include "bout/array.hxx"
#include "bout/assert.hxx"

class Field2D; // #include "field2d.hxx"
class Field3D; // #include "field3d.hxx"


class FieldPerp : public Field {
 public:
  FieldPerp();

  /*!
   * Copy constructor. After this the data
   * will be shared (non unique)
   */
  FieldPerp(const FieldPerp& f) : yindex(f.yindex),
				  nx(f.nx), nz(f.nz),
				  data(f.data) { }
  ~FieldPerp() {}

  /*!
   * Assignment operators
   */
  FieldPerp & operator=(const FieldPerp &rhs);
  FieldPerp & operator=(const BoutReal rhs);
  
  /*!
   * Iterators and data access
   */
  const DataIterator begin() const;
  const DataIterator end() const;

  inline BoutReal& operator[](DataIterator &d) {
    return operator()(d.x, d.z);
  }
  inline const BoutReal& operator[](DataIterator &d) const {
    return operator()(d.x, d.z);
  }
  BoutReal& operator[](const Indices &i) {
    return operator()(i.x, i.z);
  }
  const BoutReal& operator[](const Indices &i) const {
    return operator()(i.x, i.z);
  }
  
  int getIndex() const {return yindex;}
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

  bool isAllocated() const { return !data.empty(); }
  
  // operators

  const BoutReal* operator[](int jx) const {
    ASSERT2(!data.empty());
    ASSERT2( (jx >= 0) && (jx < nx) );
    
    return &data[jx*nz];
  }

  BoutReal* operator[](int jx) {
    ASSERT2(!data.empty());
    ASSERT2( (jx >= 0) && (jx < nx) );
    
    return &data[jx*nz];
  }

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
  
  BoutReal& operator()(int jx, int jy, int jz) { return (*this)(jx, jz); }
  
  const BoutReal& operator()(int jx, int jy, int jz) const { return (*this)(jx, jz); }

  FieldPerp & operator+=(const FieldPerp &rhs);
  FieldPerp & operator+=(const Field3D &rhs);
  FieldPerp & operator+=(const Field2D &rhs);
  FieldPerp & operator+=(const BoutReal &rhs);
  
  FieldPerp & operator-=(const FieldPerp &rhs);
  FieldPerp & operator-=(const Field3D &rhs);
  FieldPerp & operator-=(const Field2D &rhs);
  FieldPerp & operator-=(const BoutReal &rhs);

  FieldPerp & operator*=(const FieldPerp &rhs);
  FieldPerp & operator*=(const Field3D &rhs);
  FieldPerp & operator*=(const Field2D &rhs);
  FieldPerp & operator*=(const BoutReal &rhs);

  FieldPerp & operator/=(const FieldPerp &rhs);
  FieldPerp & operator/=(const Field3D &rhs);
  FieldPerp & operator/=(const Field2D &rhs);
  FieldPerp & operator/=(const BoutReal &rhs);
  
  // Stencils

  //void setStencil(bstencil *fval, bindex *bx) const;
  void setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  
 private:
  int yindex;

  int nx, nz;

  Array<BoutReal> data;
};

// Non-member overloaded operators

const FieldPerp operator+(const FieldPerp &lhs, const FieldPerp &rhs);
const FieldPerp operator+(const FieldPerp &lhs, const Field3D &rhs);
const FieldPerp operator+(const FieldPerp &lhs, const Field2D &rhs);
const FieldPerp operator+(const FieldPerp &lhs, const BoutReal &rhs);
inline const FieldPerp operator+(const BoutReal &lhs, const FieldPerp &rhs) {
  return rhs + lhs;
}

const FieldPerp operator-(const FieldPerp &lhs, const FieldPerp &other);
const FieldPerp operator-(const FieldPerp &lhs, const Field3D &other);
const FieldPerp operator-(const FieldPerp &lhs, const Field2D &other);
const FieldPerp operator-(const FieldPerp &lhs, const BoutReal &rhs);
const FieldPerp operator-(const BoutReal &lhs, const FieldPerp &rhs);

const FieldPerp operator*(const FieldPerp &lhs, const FieldPerp &other);
const FieldPerp operator*(const FieldPerp &lhs, const Field3D &other);
const FieldPerp operator*(const FieldPerp &lhs, const Field2D &other);
const FieldPerp operator*(const FieldPerp &lhs, const BoutReal &rhs);
inline const FieldPerp operator*(const BoutReal lhs, const FieldPerp &rhs) {
  return rhs * lhs;
}

const FieldPerp operator/(const FieldPerp &lhs, const FieldPerp &other);
const FieldPerp operator/(const FieldPerp &lhs, const Field3D &other);
const FieldPerp operator/(const FieldPerp &lhs, const Field2D &other);
const FieldPerp operator/(const FieldPerp &lhs, const BoutReal &rhs);
const FieldPerp operator/(const BoutReal lhs, const FieldPerp &rhs);

const FieldPerp copy(const FieldPerp &f);

/*!
 * Square a FieldPerp
 */
const FieldPerp SQ(const FieldPerp &f);

/*!
 * Create a FieldPerp by slicing a 3D field at a given y
 */
const FieldPerp sliceXZ(const Field3D& f, int y);

#endif
