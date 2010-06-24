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

#ifndef __FIELD3D_H__
#define __FIELD3D_H__

#include "field.h"
#include "field2d.h"
#include "fieldperp.h"
#include "stencils.h"
#include "bout_types.h"

/// Structure to store blocks of memory for Field3D class
struct memblock3d {
  /// memory block
  BoutReal ***data;

  /// Number of references
  int refs;
  
  /// Pointer to next block in linked-list structure
  memblock3d *next;
  
  /// Pointer in list of all blocks
  memblock3d *all_next;
}; 


/// Class for 3D X-Y-Z scalar fields
/*!
  Recycles memory using a global stack. Overloaded operators
  provide operations on data.

  July 2008: Added FieldData virtual functions
  May 2008: Added reference counting to reduce memory copying
 */
class Field3D : public Field, public FieldData {
 public:
  /// Constructor
  Field3D();
  /// copy constructor
  Field3D(const Field3D& f);
  /// Constructor from 2D field
  Field3D(const Field2D& f);
  /// Destructor
  ~Field3D();

  Field3D* clone() const;

  /// Ensures that memory is allocated
  void allocate() const;
  /// Returns a pointer to internal data (REMOVE THIS)
  BoutReal*** getData() const;
  bool isAllocated() const { return block !=  NULL; } ///< Test if data is allocated

  /// Return a pointer to the time-derivative field
  Field3D* timeDeriv();

  /// Returns DC component
  const Field2D DC();

  // Staggered grids
  void setLocation(CELL_LOC loc); // Set variable location
  CELL_LOC getLocation() const; // Variable location

  // Operators
  
  /// Allows access to internal data using square-brackets
  BoutReal** operator[](int jx) const;
  BoutReal& operator[](bindex &bx) const;

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
  Field3D & operator*=(const BoutReal rhs);
  
  /// Division
  Field3D & operator/=(const Field3D &rhs);
  Field3D & operator/=(const Field2D &rhs);
  Field3D & operator/=(const BoutReal rhs);

  /// Exponentiation (use pow() function)
  Field3D & operator^=(const Field3D &rhs);
  Field3D & operator^=(const Field2D &rhs);
  Field3D & operator^=(const BoutReal rhs);
  
  // Binary operators

  const Field3D operator+() const;
  const Field3D operator+(const Field3D &other) const;
  const Field3D operator+(const Field2D &other) const;
  const FieldPerp operator+(const FieldPerp &other) const;
  const Field3D operator+(const BoutReal &rhs) const;

  const Field3D operator-() const;
  const Field3D operator-(const Field3D &other) const;
  const Field3D operator-(const Field2D &other) const;
  const FieldPerp operator-(const FieldPerp &other) const;
  const Field3D operator-(const BoutReal &rhs) const;

  const Field3D operator*(const Field3D &other) const;
  const Field3D operator*(const Field2D &other) const;
  const FieldPerp operator*(const FieldPerp &other) const;
  const Field3D operator*(const BoutReal rhs) const;

  const Field3D operator/(const Field3D &other) const;
  const Field3D operator/(const Field2D &other) const;
  const FieldPerp operator/(const FieldPerp &other) const;
  const Field3D operator/(const BoutReal rhs) const;

  const Field3D operator^(const Field3D &other) const;
  const Field3D operator^(const Field2D &other) const;
  const FieldPerp operator^(const FieldPerp &other) const;
  const Field3D operator^(const BoutReal rhs) const;

  // Stencils for differencing

  /// Takes a location and fills values in the stencil
  void setStencil(bstencil *fval, bindex *bx) const {
    setStencil(fval, bx, true);
  }
  void setStencil(bstencil *fval, bindex *bx, bool need_x) const;

  void setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;
  void setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc = CELL_DEFAULT) const;

  /// Shifts specified points by angle
  void shiftZ(int jx, int jy, double zangle); 
  /// Shift all points in z by specified angle
  const Field3D shiftZ(const Field2D zangle) const; 
  const Field3D shiftZ(const BoutReal zangle) const;
  /// Shifts to/from BoutReal-space (using zShift global variable)
  const Field3D shiftZ(bool toBoutReal) const; 
  /// virtual function to shift between BoutReal and shifted space
  void shiftToReal(bool toBoutReal) {
    *this = shiftZ(toBoutReal);
  }
  
  // Slicing

  // NOTE: No shifting done in z for x array
  void getXArray(int y, int z, rvec &xv) const;
  void getYArray(int x, int z, rvec &yv) const;
  void getZArray(int x, int y, rvec &zv) const;

  void setXArray(int y, int z, const rvec &xv);
  void setYArray(int x, int z, const rvec &yv);
  void setZArray(int x, int y, const rvec &zv);

  /// Take a slice through the data at constant y
  const FieldPerp slice(int y) const;

  // Functions
  
  const Field3D sqrt() const;
  const Field3D abs() const;
  BoutReal min(bool allpe=false) const;
  BoutReal max(bool allpe=false) const;

  // Friend operators
  friend const Field3D operator-(const BoutReal &lhs, const Field3D &rhs);
  friend const Field3D operator+(const BoutReal &lhs, const Field3D &rhs);

  // Friend functions
  
  friend const Field3D sin(const Field3D &f);
  friend const Field3D cos(const Field3D &f);
  friend const Field3D tan(const Field3D &f);

  friend const Field3D sinh(const Field3D &f);
  friend const Field3D cosh(const Field3D &f);
  friend const Field3D tanh(const Field3D &f);

  friend const Field3D filter(const Field3D &var, int N0);
  friend const Field3D lowPass(const Field3D &var, int zmax);
  friend const Field3D lowPass(const Field3D &var, int zmax, int zmin);

  friend bool finite(const Field3D &var);

  // FieldData virtual functions
  
  bool isReal() const   { return true; }         // Consists of BoutReal values
  bool is3D() const     { return true; }         // Field is 3D
  int  byteSize() const { return sizeof(BoutReal); } // Just one BoutReal
  int  BoutRealSize() const { return 1; }
  int  getData(int x, int y, int z, void *vptr) const;
  int  getData(int x, int y, int z, BoutReal *rptr) const;
  int  setData(int x, int y, int z, void *vptr);
  int  setData(int x, int y, int z, BoutReal *rptr);

  bool ioSupport() { return true; } ///< This class supports I/O operations
  BoutReal *getData(int component) { 
    BoutReal ***d = getData();
    return **d;
  }
  void zeroComponent(int component){
    *this = 0.0;
  }

#ifdef CHECK
  bool checkData(bool vital = false) const; ///< Checks if the data is all valid. 

  void doneComms() { bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true; }
#endif

  friend class Vector3D;
  
  static void cleanup(); // Frees all memory

 private:
  /// Interpolates in z using up to 4 points
  BoutReal interpZ(int jx, int jy, int jz0, BoutReal zoffset, int order) const;

  // NOTE: Data structures mutable, though logically const

  /// Data block for this object
  mutable memblock3d *block;

  /// Number of blocks allocated
  static int nblocks;
  /// Linked list of all memory blocks
  static memblock3d *blocklist;
  /// Linked list of free blocks
  static memblock3d *free_block;

  /// Get a new block of data, either from free list or allocate
  memblock3d* newBlock() const;
  /// Makes sure data is allocated and only referenced by this object
  void allocData() const;
  /// Releases the data array, putting onto global stack
  void freeData();
  
  CELL_LOC location; // Location of the variable in the cell
  
  Field3D *ddt; ///< Time derivative (may be NULL)
};

// Non-member overloaded operators

const Field3D operator*(const BoutReal lhs, const Field3D &rhs);
const Field3D operator/(const BoutReal lhs, const Field3D &rhs);
const Field3D operator^(const BoutReal lhs, const Field3D &rhs);

// Non-member functions
const Field3D sqrt(const Field3D &f);
const Field3D abs(const Field3D &f);
BoutReal min(const Field3D &f, bool allpe=false);
BoutReal max(const Field3D &f, bool allpe=false);

#endif /* __FIELD3D_H__ */
