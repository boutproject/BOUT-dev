/*!
 * \file dcomplex.h
 * \brief Complex number class definition
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
 */
#ifndef __DCOMPLEX_H__
#define __DCOMPLEX_H__

#include "bout_types.h"

#include <iostream>
#include <fstream>

//using std::ostream;

/*!
 * \brief Complex number class
 *
 * On some machines the standard C++ complex class cannot be used because
 * there is a conflict between PVODE and the library's definition of real
 *
 * \author B.Dudson
 * \date November 2007
 */
class dcomplex {
 public:
  dcomplex();
  dcomplex(const dcomplex &rhs);
  dcomplex(real rval, real ival);
  ~dcomplex();

  dcomplex & operator=(const dcomplex &rhs);
  dcomplex & operator=(const real &rval);

  dcomplex & operator+=(const dcomplex &rhs);
  dcomplex & operator+=(const real &rhs);

  dcomplex & operator-=(const dcomplex &rhs);
  dcomplex & operator-=(const real &rhs);

  dcomplex & operator*=(const dcomplex &rhs);
  dcomplex & operator*=(const real &rhs);

  dcomplex & operator/=(const dcomplex &rhs);
  dcomplex & operator/=(const real &rhs);

  const dcomplex operator-() const; // negate

  // Binary operators

  const dcomplex operator+(const dcomplex &rhs) const;
  const dcomplex operator+(const real &rhs) const;
  
  const dcomplex operator-(const dcomplex &rhs) const;
  const dcomplex operator-(const real &rhs) const;

  const dcomplex operator*(const dcomplex &rhs) const;
  const dcomplex operator*(const real &rhs) const;

  const dcomplex operator/(const dcomplex &rhs) const;
  const dcomplex operator/(const real &rhs) const;


  friend const dcomplex operator+(const real &lhs, const dcomplex &rhs);
  friend const dcomplex operator-(const real &lhs, const dcomplex &rhs);
  friend const dcomplex operator*(const real &lhs, const dcomplex &rhs);
  friend const dcomplex operator/(const real &lhs, const dcomplex &rhs);

  // Boolean operators

  bool operator==(const dcomplex &rhs) const;
  bool operator==(const real &rhs) const;
  friend bool operator==(const real &lhs, const dcomplex &rhs);

  friend real abs(const dcomplex &c);
  friend const dcomplex conj(const dcomplex &c); // Complex conjugate 

  real Real() const;
  real Imag() const;

  // Stream operators
  friend std::ostream &operator<<(std::ostream &stream, dcomplex c);
 private:
  real r, i;

};

const dcomplex exp(const dcomplex &c);

/// imaginary i
const dcomplex Im = dcomplex(0.0, 1.0);

#endif // __DCOMPLEX_H__
