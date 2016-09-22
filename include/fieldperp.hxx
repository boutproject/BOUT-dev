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

class Field2D; // #include "field2d.hxx"
class Field3D; // #include "field3d.hxx"

class FieldPerp : public Field {
 public:
  FieldPerp();
  FieldPerp(const FieldPerp& f); // Copy constructor
  FieldPerp(BoutReal val);
  ~FieldPerp();

  DEPRECATED(FieldPerp* clone() const);
  
  void set(const Field3D &f, int y);

  void setData(BoutReal **d) {data = d;}
  BoutReal **getData() const { return data; }

  int getIndex() const;
  void setIndex(int y);

  void allocate();

  // operators

  BoutReal* operator[](int jx) const;
  BoutReal & operator()(int jx, int jz);
  const BoutReal & operator()(int jx, int jz) const;

  FieldPerp & operator=(const FieldPerp &rhs);
  FieldPerp & operator=(const BoutReal rhs);

  FieldPerp & operator+=(const FieldPerp &rhs);
  FieldPerp & operator+=(const Field3D &rhs);
  FieldPerp & operator+=(const Field2D &rhs);
  FieldPerp & operator+=(const BoutReal rhs);
  
  FieldPerp & operator-=(const FieldPerp &rhs);
  FieldPerp & operator-=(const Field3D &rhs);
  FieldPerp & operator-=(const Field2D &rhs);
  FieldPerp & operator-=(const BoutReal rhs);

  FieldPerp & operator*=(const FieldPerp &rhs);
  FieldPerp & operator*=(const Field3D &rhs);
  FieldPerp & operator*=(const Field2D &rhs);
  FieldPerp & operator*=(const BoutReal rhs);

  FieldPerp & operator/=(const FieldPerp &rhs);
  FieldPerp & operator/=(const Field3D &rhs);
  FieldPerp & operator/=(const Field2D &rhs);
  FieldPerp & operator/=(const BoutReal rhs);

  FieldPerp & operator^=(const FieldPerp &rhs);
  FieldPerp & operator^=(const Field3D &rhs);
  FieldPerp & operator^=(const Field2D &rhs);
  FieldPerp & operator^=(const BoutReal rhs);

  // Binary operators

  const FieldPerp operator+(const FieldPerp &other) const;
  const FieldPerp operator+(const Field3D &other) const;
  const FieldPerp operator+(const Field2D &other) const;
  
  const FieldPerp operator-(const FieldPerp &other) const;
  const FieldPerp operator-(const Field3D &other) const;
  const FieldPerp operator-(const Field2D &other) const;

  const FieldPerp operator*(const FieldPerp &other) const;
  const FieldPerp operator*(const Field3D &other) const;
  const FieldPerp operator*(const Field2D &other) const;
  const FieldPerp operator*(const BoutReal rhs) const;

  const FieldPerp operator/(const FieldPerp &other) const;
  const FieldPerp operator/(const Field3D &other) const;
  const FieldPerp operator/(const Field2D &other) const;
  const FieldPerp operator/(const BoutReal rhs) const;

  const FieldPerp operator^(const FieldPerp &other) const;
  const FieldPerp operator^(const Field3D &other) const;
  const FieldPerp operator^(const Field2D &other) const;
  const FieldPerp operator^(const BoutReal rhs) const;

  // Functions
  
  friend const FieldPerp exp(const FieldPerp &f);
  friend const FieldPerp log(const FieldPerp &f);
  
  friend const FieldPerp sin(const FieldPerp &f);
  friend const FieldPerp cos(const FieldPerp &f);
  friend const FieldPerp tan(const FieldPerp &f);

  friend const FieldPerp sinh(const FieldPerp &f);
  friend const FieldPerp cosh(const FieldPerp &f);
  friend const FieldPerp tanh(const FieldPerp &f);


  // Stencils

  void setStencil(bstencil *fval, bindex *bx) const;
  void setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  
 private:
  
  BoutReal interpZ(int jx, int jz0, BoutReal zoffset, int order) const;

  int yindex;
  
  BoutReal **data;

  // Data stack: Blocks of memory for this class
  static int nblocks, max_blocks;
  static BoutReal ***block; // Pointer to blocks of memory

  void allocData();
  void freeData();
};

// Non-member overloaded operators

const FieldPerp operator*(const BoutReal lhs, const FieldPerp &rhs);
const FieldPerp operator/(const BoutReal lhs, const FieldPerp &rhs);
const FieldPerp operator^(const BoutReal lhs, const FieldPerp &rhs);

const FieldPerp copy(const FieldPerp &f);

#endif
