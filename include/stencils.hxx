/*!************************************************************************
 * \file stencils.hxx
 * 
 * Sets stencils for differencing, and handles indexing
 * 
 * NOTE: These methods, classes are all deprecated, and should be replaced
 *       with DataIterator iterations and indexing
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

#ifndef __STENCILS_H__
#define __STENCILS_H__

#include "bout_types.hxx"

/// Represents a value at a single index
///
/// 
class bvalue {
 public:
  int jx, jy, jz; ///< The index (jx,jy,jz)
  BoutReal val;   ///< The value at this index

  // Operators
  
  /// Copy value and index from RHS
  bvalue & operator=(const bvalue &rhs);
  
 
  bvalue & operator+=(const bvalue &rhs);  ///< Add rhs val. Doesn't check or modify indices
  bvalue & operator-=(const bvalue &rhs);  ///< Subtract rhs val. Doesn't check or modify indices
  bvalue & operator*=(const bvalue &rhs);  ///< Multiply by rhs val. Doesn't check or modify indices
  bvalue & operator*=(BoutReal rhs); ///< Multiply by rhs val  
  bvalue & operator/=(const bvalue &rhs);  ///< Divide by rhs val. Doesn't check or modify indices
  bvalue & operator/=(BoutReal rhs); ///< Divide by rhs val

  // Binary operators
  
  const bvalue operator+(const bvalue &other) const; ///< Add bvalues. Index is taken from the left operand
  const bvalue operator-(const bvalue &other) const; ///< Subtract bvalues. Index is taken from the left operand
  const bvalue operator*(const bvalue &other) const; ///< Multiply bvalues. Index is taken from the left operand
  const bvalue operator*(BoutReal rhs) const;  ///< Multiply by value 
  const bvalue operator/(const bvalue &other) const; ///< Divide bvalues. Index is taken from the left operand
  const bvalue operator/(BoutReal rhs) const; ///< Divide by value 
};

const bvalue operator*(BoutReal lhs, const bvalue &rhs); ///< Multiply a constant by a bvalue
const bvalue operator/(BoutReal lhs, const bvalue &rhs); ///< Divide a constant by a bvalue

/// Defines a set of values in 1D in the neighbourhood of an index
/// Used for calculating derivatives
class stencil {
 public:
  int jx, jy, jz; ///< central location
  
  BoutReal c, p, m, pp, mm; ///< stencil 2 each side of the centre
  
  /// constructor
  stencil();
  /// Constructor like pre r115 upwind methods (debugging)
  stencil(BoutReal fc);
  /// Constructor like pre r115 derivative methods
  stencil(BoutReal fc, BoutReal fm, BoutReal fp, BoutReal fmm, BoutReal fpp);

  /// Copy constructor
  stencil(const stencil &s);

  // operators
  
  stencil & operator=(const stencil &rhs);
  stencil & operator=(BoutReal rhs);

  stencil & operator+=(const stencil &rhs);
  stencil & operator+=(BoutReal rhs);
  stencil & operator-=(const stencil &rhs);
  stencil & operator-=(BoutReal rhs);
  stencil & operator*=(const stencil &rhs);
  stencil & operator*=(BoutReal rhs);
  stencil & operator/=(const stencil &rhs);
  stencil & operator/=(BoutReal rhs);

  const stencil operator+(const stencil &other) const;
  const stencil operator+(BoutReal other) const;
  const stencil operator-(const stencil &other) const;
  const stencil operator-(BoutReal other) const;
  const stencil operator*(const stencil &other) const;
  const stencil operator*(BoutReal other) const;
  const stencil operator/(const stencil &other) const;
  const stencil operator/(BoutReal other) const;

  friend const stencil operator+(BoutReal lhs, const stencil &rhs);
  friend const stencil operator-(BoutReal lhs, const stencil &rhs);
  friend const stencil operator*(BoutReal lhs, const stencil &rhs);
  friend const stencil operator/(BoutReal lhs, const stencil &rhs);
  
  BoutReal min() const; ///< Minimum value over stencil
  BoutReal max() const; ///< Maximum value over stencil
  const stencil abs() const; ///< Stencil with all values replaced by absolute values
};

class forward_stencil {
/// to be used only in lower/inner boundary guard cells
 public:
  int jx, jy, jz; ///< central location
  
  BoutReal m, c, p, p2, p3, p4; ///< stencil 4 forward of the centre
  
  /// constructor
  forward_stencil();

};

/// to be used only in upper/outer boundary guard cells
class backward_stencil {
 public:
  int jx, jy, jz; /// central location
  
  BoutReal p, c , m, m2, m3, m4; /// stencil 4 backward of the centre
  
  /// constructor
  backward_stencil();

};

BoutReal min(const stencil &s); ///< Minimum value over a stencil
BoutReal max(const stencil &s); ///< Maximum value over a stencil
const stencil abs(const stencil &s);  ///< Copy of stencil with all values replaced by absolute values

/// Represents indexes for stencils
/// This allowed for complicated connections between points
///
///
typedef struct {
  int jx,jy,jz, /* center   */
    
    /*1st order neigbors*/
    jxm, jxp,   /* (+/-) dx */
    jym, jyp,   /* (+/-) dy */
    jzm, jzp,   /* (+/-) dz */
    
    /*2nd order neigbors*/ 
    jx2m, jx2p,  /* (+/-) 2*dx */
    jy2m, jy2p,  /* (+/-) 2*dy */
    jz2m, jz2p;  /* (+/-) 2*dz */

  BoutReal yp_offset, ym_offset; ///< Toroidal offset (In index space)
  
  /// Offsets for sheared X derivatives (also in index space)
  BoutReal xp_offset, xm_offset;
  BoutReal x2p_offset, x2m_offset;

  /// What region is being looped over?
  REGION region;
  
} bindex;

void calc_index(bindex *bx); ///< Calculate all other indices in bindex, given central index

/// Start an iteration over indices
/// 
/// NOTE: This is done much better and more efficiently by DataIterator
/// 
/// @param[out] bx   The bindex which will be modified and used for iterating
/// @param[in]  region  The region over which to iterate. By default doesn't loop over boundary regions
///
void start_index(bindex *bx, REGION region = RGN_NOBNDRY);

/// move to next index in 3D
int next_index3(bindex *bx);

/// Move to next index in 2D
int next_index2(bindex *bx);

/// Move to next index in X-Z
int next_indexperp(bindex *bx);

void reverse_start_index(bindex *bx, REGION region = RGN_NOBNDRY); // Doesn't loop over boundary regions
void start_index_lasty(bindex *bx, REGION region = RGN_NOBNDRY);
int reverse_next_index3(bindex *bx);
int next_index_y(bindex *bx);
int previous_index_y(bindex *bx);

#endif /* __STENCILS_H__ */
