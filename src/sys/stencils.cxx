/**************************************************************************
 * Stencil class for differencing, and handles indexing
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

#include <cmath>

typedef double BoutReal;

#include <globals.hxx>
#include <stencils.hxx>
#include <output.hxx>

/*******************************************************************************
 * stencil class
 *******************************************************************************/

stencil::stencil(const stencil &s)
{
  *this = s;
}

stencil::stencil(BoutReal fc)
{
  c = p = m = pp = mm = fc;
}

stencil::stencil(BoutReal fc, BoutReal fm, BoutReal fp, BoutReal fmm, BoutReal fpp)
{
  c  = fc;
  m  = fm;
  p  = fp;
  mm = fmm;
  pp = fpp;
}

stencil & stencil::operator=(const stencil &s)
{
  // Check for self-assignment
  if(this == &s)
    return(*this); // skip this assignment
  
  jx = s.jx;
  jy = s.jy;
  jz = s.jz;

  c = s.c;
  p = s.p;
  m = s.m;
  pp = s.pp;
  mm = s.mm;
  
  return *this;
}

stencil & stencil::operator=(const BoutReal rhs)
{
  jx = jy = jz = 0;

  c = p = m = pp = mm = rhs;

  return *this;
}

stencil & stencil::operator+=(const stencil &s)
{
#if CHECK > 0
  if((jx != s.jx) || (jy != s.jy) || (jz != s.jz)) {
    throw BoutException("Error: Stencil += not at same location\n");
  }
#endif

  c += s.c;
  p += s.p;
  m += s.m;
  pp += s.pp;
  mm += s.mm;

  return *this;
}

stencil & stencil::operator+=(BoutReal rhs)
{
  c += rhs;
  p += rhs;
  m += rhs;
  pp += rhs;
  mm += rhs;

  return *this;
}

stencil & stencil::operator-=(const stencil &s)
{
#if CHECK > 0
  if((jx != s.jx) || (jy != s.jy) || (jz != s.jz)) {
    throw BoutException("Error: Stencil -= not at same location\n");
  }
#endif

  c -= s.c;
  p -= s.p;
  m -= s.m;
  pp -= s.pp;
  mm -= s.mm;

  return *this;
}

stencil & stencil::operator-=(BoutReal rhs)
{
  c -= rhs;
  p -= rhs;
  m -= rhs;
  pp -= rhs;
  mm -= rhs;

  return *this;
}

stencil & stencil::operator*=(const stencil &s)
{
#if CHECK > 0
  if((jx != s.jx) || (jy != s.jy) || (jz != s.jz)) {
    throw BoutException("Error: Stencil *= not at same location\n");
  }
#endif

  c *= s.c;
  p *= s.p;
  m *= s.m;
  pp *= s.pp;
  mm *= s.mm;

  return *this;
}

stencil & stencil::operator*=(BoutReal rhs)
{
  c *= rhs;
  p *= rhs;
  m *= rhs;
  pp *= rhs;
  mm *= rhs;

  return *this;
}

stencil & stencil::operator/=(const stencil &s)
{
#if CHECK > 0
  if((jx != s.jx) || (jy != s.jy) || (jz != s.jz)) {
    throw BoutException("Error: Stencil /= not at same location\n");
  }
#endif

  c /= s.c;
  p /= s.p;
  m /= s.m;
  pp /= s.pp;
  mm /= s.mm;

  return *this;
}

stencil & stencil::operator/=(BoutReal rhs)
{
  c /= rhs;
  p /= rhs;
  m /= rhs;
  pp /= rhs;
  mm /= rhs;

  return *this;
}

const stencil stencil::operator+(const stencil &other) const
{
  stencil result = *this;
  result += other;
  return result;
}

const stencil stencil::operator+(BoutReal other) const
{
  stencil result = *this;
  result += other;
  return result;
}

const stencil stencil::operator-(const stencil &other) const
{
  stencil result = *this;
  result -= other;
  return result;
}

const stencil stencil::operator-(BoutReal other) const
{
  stencil result = *this;
  result -= other;
  return result;
}

const stencil stencil::operator*(const stencil &other) const
{
  stencil result = *this;
  result *= other;
  return result;
}

const stencil stencil::operator*(BoutReal other) const
{
  stencil result = *this;
  result *= other;
  return result;
}

const stencil stencil::operator/(const stencil &other) const
{
  stencil result = *this;
  result /= other;
  return result;
}

const stencil stencil::operator/(BoutReal other) const
{
  stencil result = *this;
  result /= other;
  return result;
}

const stencil operator+(BoutReal lhs, const stencil &rhs)
{
  return rhs + lhs;
}

const stencil operator-(BoutReal lhs, const stencil &rhs)
{
  stencil result;
  
  result.c = lhs - rhs.c;
  result.p = lhs - rhs.p;
  result.m = lhs - rhs.m;
  result.pp = lhs - rhs.pp;
  result.mm = lhs - rhs.mm;

  return result;
}

const stencil operator*(BoutReal lhs, const stencil &rhs)
{
  return rhs * lhs;
}

const stencil operator/(BoutReal lhs, const stencil &rhs)
{
  stencil result;
  
  result.c = lhs / rhs.c;
  result.p = lhs / rhs.p;
  result.m = lhs / rhs.m;
  result.pp = lhs / rhs.pp;
  result.mm = lhs / rhs.mm;

  return result;
}

BoutReal stencil::min() const
{
  BoutReal r;

  r = c;
  if(p < r)
    r = p;
  if(m < r)
    r = m;
  if(pp < r)
    r = pp;
  if(mm < r)
    r = mm;
  return r;
}

BoutReal stencil::max() const
{
  BoutReal r;

  r = c;
  if(p > r)
    r = p;
  if(m > r)
    r = m;
  if(pp > r)
    r = pp;
  if(mm > r)
    r = mm;
  return r;
}

const stencil stencil::abs() const
{
  stencil result;
  
  result.c = fabs(c);
  result.p = fabs(p);
  result.m = fabs(m);
  result.pp = fabs(pp);
  result.mm = fabs(mm);

  return result;
}

BoutReal min(const stencil &s)
{
  return s.min();
}

BoutReal max(const stencil &s)
{
  return s.max();
}

const stencil abs(const stencil &s)
{
  return s.abs();
}
