/**************************************************************************
 * Sets stencils for differencing, and handles indexing
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

#include "bout_types.h"

class bvalue {
 public:
  int jx, jy, jz;
  real val;

  // Operators
  
  bvalue & operator=(const bvalue &rhs);
  
  bvalue & operator+=(const bvalue &rhs);
  bvalue & operator-=(const bvalue &rhs);
  bvalue & operator*=(const bvalue &rhs);
  bvalue & operator*=(const real rhs);
  bvalue & operator/=(const bvalue &rhs);
  bvalue & operator/=(const real rhs);

  // Binary operators
  
  const bvalue operator+(const bvalue &other) const;
  const bvalue operator-(const bvalue &other) const;
  const bvalue operator*(const bvalue &other) const;
  const bvalue operator*(const real rhs) const;
  const bvalue operator/(const bvalue &other) const;
  const bvalue operator/(const real rhs) const;
};

const bvalue operator*(const real lhs, const bvalue &rhs);
const bvalue operator/(const real lhs, const bvalue &rhs);

class stencil {
 public:
  int jx, jy, jz; // central location
  
  real c, p, m, pp, mm; // stencil 2 each side of the centre
  
  /// constructor
  stencil();
  /// Constructor like pre r115 upwind methods (debugging)
  stencil(real fc);
  /// Constructor like pre r115 derivative methods
  stencil(real fc, real fm, real fp, real fmm, real fpp);

  /// Copy constructor
  stencil(const stencil &s);

  // operators
  
  stencil & operator=(const stencil &rhs);
  stencil & operator=(const real rhs);

  stencil & operator+=(const stencil &rhs);
  stencil & operator+=(const real &rhs);
  stencil & operator-=(const stencil &rhs);
  stencil & operator-=(const real &rhs);
  stencil & operator*=(const stencil &rhs);
  stencil & operator*=(const real &rhs);
  stencil & operator/=(const stencil &rhs);
  stencil & operator/=(const real &rhs);

  const stencil operator+(const stencil &other) const;
  const stencil operator+(const real &other) const;
  const stencil operator-(const stencil &other) const;
  const stencil operator-(const real &other) const;
  const stencil operator*(const stencil &other) const;
  const stencil operator*(const real &other) const;
  const stencil operator/(const stencil &other) const;
  const stencil operator/(const real &other) const;

  friend const stencil operator+(const real &lhs, const stencil &rhs);
  friend const stencil operator-(const real &lhs, const stencil &rhs);
  friend const stencil operator*(const real &lhs, const stencil &rhs);
  friend const stencil operator/(const real &lhs, const stencil &rhs);

  real min() const;
  real max() const;
  const stencil abs() const;
};

real min(const stencil &s);
real max(const stencil &s);
const stencil abs(const stencil &s);

class bstencil {
 public:
  
  int jx, jy, jz; // central location

  // Values

  real cc,  // center
    
    // 1st order neigbors
    xp, xm, /* (+/-) dx */
    yp, ym, /* (+/-) dy */
    zp, zm, /* (+/-) dz */
    
    // 2nd order neigbors  
    x2p, x2m, /* (+/-) 2*dx */
    y2p, y2m, /* (+/-) 2*dy */
    z2p, z2m; /* (+/-) 2*dz */
  
  // Operators

  bstencil & operator=(const bstencil &rhs);
  
  bstencil & operator+=(const bstencil &rhs);
  bstencil & operator-=(const bstencil &rhs);
  bstencil & operator*=(const bstencil &rhs);
  bstencil & operator*=(const real rhs);
  bstencil & operator/=(const bstencil &rhs);
  bstencil & operator/=(const real rhs);
  bstencil & operator^=(const bstencil &rhs);
  bstencil & operator^=(const real rhs);

  // Binary operators

  const bstencil operator+(const bstencil &other) const;
  const bstencil operator-(const bstencil &other) const;
  const bstencil operator*(const bstencil &other) const;
  const bstencil operator*(const real rhs) const;
  const bstencil operator/(const bstencil &other) const;
  const bstencil operator/(const real rhs) const;
  const bstencil operator^(const bstencil &other) const;
  const bstencil operator^(const real rhs) const;
  
 private:
  
  
};

const bstencil operator*(const real lhs, const bstencil &rhs);
const bstencil operator/(const real lhs, const bstencil &rhs);
const bstencil operator^(const real lhs, const bstencil &rhs);

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
  
  // Switches for twist-shift condition
  bool yp_shift, y2p_shift;
  bool ym_shift, y2m_shift;

  real yp_offset, ym_offset; // Toroidal offset (In index space)
  
  // Offsets for sheared X derivatives (also in index space)
  real xp_offset, xm_offset;
  real x2p_offset, x2m_offset;

  // What region is being looped over?
  REGION region;
  
} bindex;

void calc_index(bindex *bx);
void start_index(bindex *bx, REGION region = RGN_NOBNDRY); // Doesn't loop over boundary regions
int next_index3(bindex *bx);
int next_index2(bindex *bx);
int next_indexperp(bindex *bx);

#endif /* __STENCILS_H__ */
