/*!*************************************************************************
 * \file invert3x3.hxx
 *
 * A mix of short utilities for memory management, strings, and some
 * simple but common calculations
 *
 **************************************************************************
 * Copyright 2010-2024 B.D.Dudson, BOUT++ Team
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#pragma once

#include <bout/utils.hxx>

/// Explicit inversion of a 3x3 matrix \p a
///
/// The input \p small determines how small the determinant must be for
/// us to throw due to the matrix being singular (ill conditioned);
/// If small is less than zero then instead of throwing we return 1.
/// This is ugly but can be used to support some use cases.
template <typename T>
int invert3x3(Matrix<T>& a, BoutReal small = 1.0e-15) {
  TRACE("invert3x3");

  // Calculate the first co-factors
  T A = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1);
  T B = a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2);
  T C = a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0);

  // Calculate the determinant
  T det = a(0, 0) * A + a(0, 1) * B + a(0, 2) * C;

  if (std::abs(det) < std::abs(small)) {
    if (small >= 0) {
      throw BoutException("Determinant of matrix < {:e} --> Poorly conditioned", small);
    } else {
      return 1;
    }
  }

  // Calculate the rest of the co-factors
  T D = a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2);
  T E = a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0);
  T F = a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1);
  T G = a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1);
  T H = a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2);
  T I = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);

  // Now construct the output, overwrites input
  T detinv = 1.0 / det;

  a(0, 0) = A * detinv;
  a(0, 1) = D * detinv;
  a(0, 2) = G * detinv;
  a(1, 0) = B * detinv;
  a(1, 1) = E * detinv;
  a(1, 2) = H * detinv;
  a(2, 0) = C * detinv;
  a(2, 1) = F * detinv;
  a(2, 2) = I * detinv;

  return 0;
}
