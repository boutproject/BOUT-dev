/**************************************************************************
 * Complex number class
 * On some machines the standard C++ complex class cannot be used because
 * there is a conflict between PVODE and the library's definition of BoutReal
 *
 * Changelog:
 * 
 * 2007-11 Ben Dudson <bd512@york.ac.uk>
 *    * Initial version
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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
#include <iostream>

#include <dcomplex.hxx>
#include <utils.hxx>

dcomplex::dcomplex()
{
}

dcomplex::dcomplex(const dcomplex &rhs)
{
  *this = rhs;
}

dcomplex::~dcomplex()
{
}

dcomplex & dcomplex::operator=(const dcomplex &rhs)
{
  r = rhs.r; i = rhs.i;
  return *this;
}

dcomplex & dcomplex::operator=(const BoutReal &rhs)
{
  r = rhs; i = 0.0;
  return *this;
}

dcomplex & dcomplex::operator+=(const dcomplex &rhs)
{
  r += rhs.r;
  i += rhs.i;
  return *this;
}

dcomplex & dcomplex::operator+=(const BoutReal &rhs)
{
  r += rhs;
  return *this;
}

dcomplex & dcomplex::operator-=(const dcomplex &rhs)
{
  r -= rhs.r;
  i -= rhs.i;
  return *this;
}

dcomplex & dcomplex::operator-=(const BoutReal &rhs)
{
  r -= rhs;
  return *this;
}

dcomplex & dcomplex::operator*=(const dcomplex &rhs)
{
  BoutReal rt;

  rt = r * rhs.r - i * rhs.i;

  i = i*rhs.r + r*rhs.i;
  r = rt;
  return *this;
}

dcomplex & dcomplex::operator*=(const BoutReal &rhs)
{
  r *= rhs;
  i *= rhs;
  return *this;
}

dcomplex & dcomplex::operator/=(const dcomplex &rhs)
{
  BoutReal c, rt;

  c = rhs.r*rhs.r + rhs.i*rhs.i;
  
  rt = (r*rhs.r + i*rhs.i)/c;
  i = (i*rhs.r - r*rhs.i)/c;
  r = rt;
  return *this;
}

dcomplex & dcomplex::operator/=(const BoutReal &rhs) {
  r /= rhs;
  i /= rhs;
  return *this;
}

const dcomplex dcomplex::operator-() const
{
  dcomplex result;

  result.r = -r;
  result.i = -i;
  
  return result;
}

const dcomplex dcomplex::operator+(const dcomplex &rhs) const
{
  dcomplex result;
  
  result.r = r + rhs.r;
  result.i = i + rhs.i;

  return result;
}

const dcomplex dcomplex::operator+(const BoutReal &rhs) const
{
  dcomplex result;

  result.r = r + rhs;
  result.i = i;

  return result;
}

const dcomplex dcomplex::operator-(const dcomplex &rhs) const
{
  dcomplex result;
  
  result.r = r - rhs.r;
  result.i = i - rhs.i;

  return result;
}

const dcomplex dcomplex::operator-(const BoutReal &rhs) const
{
  dcomplex result;

  result.r = r - rhs;
  result.i = i;

  return result;
}

const dcomplex dcomplex::operator*(const dcomplex &rhs) const
{
  dcomplex result;

  result.r = r*rhs.r - i*rhs.i;
  result.i = i*rhs.r + r*rhs.i;
  
  return result;
}

const dcomplex dcomplex::operator*(const BoutReal &rhs) const
{
  dcomplex result;

  result.r = r * rhs;
  result.i = i * rhs;

  return result;
}

const dcomplex dcomplex::operator/(const dcomplex &rhs) const
{
  dcomplex result;
  BoutReal c;

  c = rhs.r*rhs.r + rhs.i*rhs.i;
  result.r = (r*rhs.r + i*rhs.i)/c;
  result.i = (i*rhs.r - r*rhs.i)/c;

  return result;
}

const dcomplex dcomplex::operator/(const BoutReal &rhs) const
{
  dcomplex result;
  
  result.r = r / rhs;
  result.i = i / rhs;

  return result;
}

// Non-member overloaded operators

const dcomplex operator+(const BoutReal &lhs, const dcomplex &rhs)
{
  return rhs+lhs;
}

const dcomplex operator-(const BoutReal &lhs, const dcomplex &rhs)
{
  dcomplex result;

  result.r = lhs - rhs.r;
  result.i = -rhs.i;

  return result;
}

const dcomplex operator*(const BoutReal &lhs, const dcomplex &rhs)
{
  return rhs*lhs;
}

const dcomplex operator/(const BoutReal &lhs, const dcomplex &rhs)
{
  dcomplex result;
  BoutReal c;
  
  c = rhs.r*rhs.r + rhs.i*rhs.i;

  result.r = lhs*rhs.r / c;
  result.i = -lhs*rhs.i / c;

  return result;
}

// Boolean operators

bool dcomplex::operator==(const dcomplex &rhs) const
{
  return (r == rhs.r) && (i == rhs.i);
}

bool dcomplex::operator==(const BoutReal &rhs) const
{
  return (r == rhs) && (i == 0.0);
}

bool operator==(const BoutReal &lhs, const dcomplex &rhs)
{
  return (lhs == rhs.r) && (rhs.i == 0);
}

BoutReal abs(const dcomplex &c)
{
  return sqrt(SQ(c.r) + SQ(c.i));
}

const dcomplex conj(const dcomplex &c)
{
  dcomplex result;
  
  result.r = c.r;
  result.i = -c.i;

  return result;
}

BoutReal dcomplex::Real() const
{
  return r;
}

BoutReal dcomplex::Imag() const
{
  return i;
}

/****************************************************************
 * Functions
 ****************************************************************/

const dcomplex exp(const dcomplex &c)
{
  return exp(c.Real()) * dcomplex(cos(c.Imag()), sin(c.Imag()));
}

std::ostream &operator<<(std::ostream &stream, dcomplex c)
{
  stream << "(" << c.r << ", " << c.i << ")";
  
  return stream;
}

