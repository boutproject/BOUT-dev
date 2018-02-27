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

BoutReal min(const stencil &s); ///< Minimum value over a stencil
BoutReal max(const stencil &s); ///< Maximum value over a stencil
const stencil abs(const stencil &s);  ///< Copy of stencil with all values replaced by absolute values

#endif /* __STENCILS_H__ */
