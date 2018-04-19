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

#include "bout/dataiterator.hxx"

#include "bout/array.hxx"
#include "bout/region.hxx"

#include "bout/deprecated.hxx"
#include "bout/assert.hxx"

#include "bout/field_visitor.hxx"

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

  Indexing can also be done using DataIterator or Indices objects,
  defined in bout/dataiterator.hxx:

      Indices i = {0,1,0};

      f[i] = 1.0;  // Equivalent to f(0,1,0)

  This is primarily used to allow convenient iteration over fields

  Iteration
  ---------

  To loop over all points in a field, a for loop can be used
  to get the indices:

      Field3D f(0.0); // Allocate, set to zero

      for( auto i : f ) {  // Loop over all points, with index i
        f[i] = 1.0;
      }

  There is also more explicit looping over regions:

      for( auto i : f.region(RGN_ALL) ) {  // Loop over all points, with index i
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

      f.mergeYupYdown(); // f.yup() and f.ydown() now point to f
      f.yup()(0,1,0)  // ok, gives value of f at (0,1,0)

      To have separate fields for yup and ydown, first call

      f.splitYupYdown(); // f.yup() and f.ydown() separate

      f.yup(); // ok
      f.yup()(0,1,0) // error; f.yup not allocated

      f.yup() = 1.0; // Set f.yup() field to 1.0

      f.yup()(0,1,0) // ok

 */
class Field3D : public Field, public FieldData {
 public:
  /*!
   * Constructor
   *
   * Note: the global "mesh" can't be passed here because
   * fields may be created before the mesh is.
   */
  Field3D(Mesh *localmesh = nullptr);

  /*!
   * Copy constructor
   */
  Field3D(const Field3D& f);
  
  /// Constructor from 2D field
  Field3D(const Field2D& f);
  /// Constructor from value
  Field3D(BoutReal val, Mesh *localmesh = nullptr);
  /// Destructor
  ~Field3D();

  /// Data type stored in this field
  using value_type = BoutReal;

  /*!
   * Ensures that memory is allocated and unique
   */
  void allocate();
  
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

  /*!
   * Ensure that this field has separate fields
   * for yup and ydown.
   */
  void splitYupYdown();

  /*!
   * Ensure that yup and ydown refer to this field
   */
  void mergeYupYdown();
  
  /// Check if this field has yup and ydown fields
  bool hasYupYdown() const {
    return (yup_field != nullptr) && (ydown_field != nullptr);
  }

  /// Return reference to yup field
  Field3D& yup() { 
    ASSERT2(yup_field != nullptr); // Check for communicate
    return *yup_field; 
  }
  /// Return const reference to yup field
  const Field3D& yup() const { 
    ASSERT2(yup_field != nullptr);
    return *yup_field; 
  }
  
  /// Return reference to ydown field
  Field3D& ydown() { 
    ASSERT2(ydown_field != nullptr);
    return *ydown_field;
  }
  
  /// Return const reference to ydown field
  const Field3D& ydown() const { 
    ASSERT2(ydown_field != nullptr);
    return *ydown_field; 
  }

  /// Return yup if dir=+1, and ydown if dir=-1
  Field3D& ynext(int dir);
  const Field3D& ynext(int dir) const;

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
  
  const DataIterator iterator() const;

  /*!
   * These begin and end functions are used to iterate over
   * the indices of a field. Indices are used rather than
   * values since this allows expressions involving multiple fields.
   *
   * Example
   * -------
   *
   * Field3D objects f and g can be modified by 
   * 
   * for(auto i : f) {
   *   f[i] = 2.*f[i] + g[i];
   * }
   * 
   */
  const DataIterator begin() const;
  const DataIterator end() const;
  
  /*!
   * Returns a range of indices which can be iterated over
   * Uses the REGION flags in bout_types.hxx
   * 
   * Example
   * -------
   * 
   * This loops over the interior points, not the boundary
   * and inside the loop the index is used to calculate the difference
   * between the point one index up in x (i.xp()) and one index down
   * in x (i.xm()), putting the result into a different field 'g'
   * 
   * for(auto i : f.region(RGN_NOBNDRY)) {
   *   g[i] = f[i.xp()] - f[i.xm()];
   * }
   * 
   */
  const IndexRange region(REGION rgn) const override;

  /*!
   * Like Field3D::region(REGION rgn), but returns range
   * to iterate over only x-y, not z.
   * This is useful in the Fourier transform functions
   * which need an explicit loop in z.
   *
   */
  const IndexRange region2D(REGION rgn) const;

  /*!
   * Direct data access using DataIterator object.
   * This uses operator(x,y,z) so checks will only be
   * performed if CHECK > 2.
   */
  BoutReal& operator[](const DataIterator &d) {
    return operator()(d.x, d.y, d.z);
  }
  const BoutReal& operator[](const DataIterator &d) const {
    return operator()(d.x, d.y, d.z);
  }
  BoutReal& operator[](const Indices &i) override {
    return operator()(i.x, i.y, i.z);
  }
  const BoutReal& operator[](const Indices &i) const override {
    return operator()(i.x, i.y, i.z);
  }
  

  BoutReal& operator[](const Ind3D &d) {
    return data[d.ind];
  }
  const BoutReal& operator[](const Ind3D &d) const {
    return data[d.ind];
  }

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
      throw BoutException("Field3D: (%d, %d, %d) operator out of bounds (%d, %d, %d)", 
			  jx, jy, jz, nx, ny, nz);
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
      throw BoutException("Field3D: (%d, %d, %d) operator out of bounds (%d, %d, %d)", 
			  jx, jy, jz, nx, ny, nz);
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
      throw BoutException("Field3D: (%d, %d) operator out of bounds (%d, %d)",
                          jx, jy, nx, ny);
#endif
    return &data[(jx*ny +jy)*nz];
  }

  inline BoutReal* operator()(int jx, int jy) {
#if CHECK > 2
    if(data.empty())
      throw BoutException("Field3D: () operator on empty data");

    if((jx < 0) || (jx >= nx) ||
       (jy < 0) || (jy >= ny))
      throw BoutException("Field3D: (%d, %d) operator out of bounds (%d, %d)",
                          jx, jy, nx, ny);
#endif
    return &data[(jx*ny +jy)*nz];
  }
  
  /////////////////////////////////////////////////////////
  // Operators
  
  const Field3D operator+() const {return *this;}
  
  /// Assignment operators
  ///@{
  Field3D & operator=(const Field3D &rhs);
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

  DEPRECATED(void setBackground(const Field2D &f2d)); ///< Boundary is applied to the total of this and f2d
  void applyBoundary(bool init=false) override;
  void applyBoundary(BoutReal t);
  void applyBoundary(const string &condition);
  void applyBoundary(const char* condition) { applyBoundary(string(condition)); }
  void applyBoundary(const string &region, const string &condition);
  void applyTDerivBoundary() override;
  void setBoundaryTo(const Field3D &f3d); ///< Copy the boundary region

  void applyParallelBoundary();
  void applyParallelBoundary(BoutReal t);
  void applyParallelBoundary(const string &condition);
  void applyParallelBoundary(const char* condition) { applyParallelBoundary(string(condition)); }
  void applyParallelBoundary(const string &region, const string &condition);
  void applyParallelBoundary(const string &region, const string &condition, Field3D *f);
  
private:
  /// Boundary - add a 2D field
  const Field2D *background;

  int nx, ny, nz;  ///< Array sizes (from fieldmesh). These are valid only if fieldmesh is not null
  
  /// Internal data array. Handles allocation/freeing of memory
  Array<BoutReal> data;

  CELL_LOC location; ///< Location of the variable in the cell
  
  Field3D *deriv; ///< Time derivative (may be NULL)

  /// Pointers to fields containing values along Y
  Field3D *yup_field, *ydown_field;
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

/*!
 * Calculates the minimum of a field, excluding
 * the boundary/guard cells by default (can be changed
 * with rgn argument).
 * By default this is only on the local processor,
 * but setting allpe=true does a collective Allreduce
 * over all processors.
 *
 * @param[in] f  The field to loop over
 * @param[in] allpe  Minimum over all processors?
 * @param[in] rgn  the boundaries that should be ignored
 * 
 */
BoutReal min(const Field3D &f, bool allpe=false, REGION rgn=RGN_NOBNDRY);

/*!
 * Calculates the maximum of a field, excluding
 * the boundary/guard cells by default (can be changed
 * with rgn argument).
 * By default this is only on the local processor,
 * but setting allpe=true does a collective Allreduce
 * over all processors.
 *
 * @param[in] f  The field to loop over
 * @param[in] allpe  Minimum over all processors?
 * @param[in] rgn  the boundaries that should be ignored
 * 
 */
BoutReal max(const Field3D &f, bool allpe=false, REGION rgn=RGN_NOBNDRY);

/*!
 * Exponent: pow(a, b) is a raised to the power of b
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
Field3D pow(const Field3D &lhs, const Field3D &rhs, REGION rgn = RGN_ALL);
Field3D pow(const Field3D &lhs, const Field2D &rhs, REGION rgn = RGN_ALL);
Field3D pow(const Field3D &lhs, const FieldPerp &rhs, REGION rgn = RGN_ALL);
Field3D pow(const Field3D &lhs, BoutReal rhs, REGION rgn = RGN_ALL);
Field3D pow(BoutReal lhs, const Field3D &rhs, REGION rgn = RGN_ALL);

/*!
 * Square root
 * 
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D sqrt(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Absolute value (modulus, |f|)
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D abs(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Exponential: exp(f) is e to the power of f
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D exp(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Natural logarithm, inverse of exponential
 * 
 *     log(exp(f)) = f
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D log(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Sine trigonometric function. 
 *
 * @param[in] f  Angle in radians
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D sin(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Cosine trigonometric function. 
 *
 * @param[in] f  Angle in radians
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D cos(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Tangent trigonometric function. 
 *
 * @param[in] f  Angle in radians
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D tan(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Hyperbolic sine function. 
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D sinh(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Hyperbolic cosine function. 
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D cosh(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Hyperbolic tangent function. 
 *
 * This loops over the entire domain, including guard/boundary cells by
 * default (can be changed using the rgn argument)
 * If CHECK >= 3 then the result will be checked for non-finite numbers
 */
const Field3D tanh(const Field3D &f, REGION rgn = RGN_ALL);

/*!
 * Check if all values of a field are finite.
 * Loops over all points including the boundaries
 */
bool finite(const Field3D &var, REGION rgn = RGN_ALL);


#if CHECK > 0
/// Throw an exception if \p f is not allocated or if any
/// elements are non-finite (for CHECK > 2)
void checkData(const Field3D &f, REGION region = RGN_NOBNDRY);
#else
/// Ignored with disabled CHECK; Throw an exception if \p f is not
/// allocated or if any elements are non-finite (for CHECK > 2)
inline void checkData(const Field3D &UNUSED(f), REGION UNUSED(region) = RGN_NOBNDRY){};
#endif
 
/*!
 * Makes a copy of a field, ensuring that the underlying
 * data is not shared.
 */ 
const Field3D copy(const Field3D &f);

/*!
 * Apply a floor value to a field. Any value lower than
 * the floor is set to the floor.
 * 
 * @param[in] var  Variable to apply floor to
 * @param[in] f    The floor value
 *
 */
const Field3D floor(const Field3D &var, BoutReal f, REGION rgn = RGN_ALL);

/*!
 * Fourier filtering, removes all except one mode
 * 
 * @param[in] var Variable to apply filter to
 * @param[in] N0 The component to keep
 */
const Field3D filter(const Field3D &var, int N0, REGION rgn=RGN_ALL);

/*!
 * Fourier low pass filtering. Removes modes higher than zmax
 */ 
const Field3D lowPass(const Field3D &var, int zmax, REGION rgn=RGN_ALL);

/*!
 * Fourier low pass filtering. Removes modes
 * lower than zmin and higher than zmax
 */
const Field3D lowPass(const Field3D &var, int zmax, int zmin, REGION rgn=RGN_ALL);

/*!
 * Perform a shift by a given angle in Z
 *
 * @param[inout] var  The variable to be modified in-place
 * @param[in] jx   X index
 * @param[in] jy   Y index
 * @param[in] zangle   The Z angle to apply
 */
void shiftZ(Field3D &var, int jx, int jy, double zangle);

/*!
 * Apply a phase shift by a given angle in Z to all points
 * 
 * @param[inout] var  The variable to modify in-place
 * @param[in] zangle  The angle to shift by in Z
 */
void shiftZ(Field3D &var, double zangle, REGION rgn=RGN_ALL);

/*!
 * Average in the Z direction
 */ 
Field2D DC(const Field3D &f, REGION rgn=RGN_ALL);

/*!
 * Force guard cells of passed field to nan
 */ 
void invalidateGuards(Field3D &var);

/*!
 * @brief Returns a reference to the time-derivative of a field
 * 
 * Wrapper around member function f.timeDeriv()
 *
 */
inline Field3D& ddt(Field3D &f) {
  return *(f.timeDeriv());
}

#endif /* __FIELD3D_H__ */
