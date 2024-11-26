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
#include <optional>

/// Explicit inversion of a 3x3 matrix \p a
///
/// If the matrix is singular (ill conditioned), the determinant is
/// return. Otherwise, an empty `std::optional` is return
namespace bout {
inline std::optional<BoutReal> invert3x3(Matrix<BoutReal>& a) {
  TRACE("invert3x3");

  // Calculate the first co-factors
  const BoutReal A = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1);
  const BoutReal B = a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2);
  const BoutReal C = a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0);

  // Calculate the determinant
  const BoutReal det = a(0, 0) * A + a(0, 1) * B + a(0, 2) * C;
  constexpr BoutReal small = 1.0e-15;
  if (std::abs(det) < std::abs(small)) {
    return det;
  }

  // Calculate the rest of the co-factors
  const BoutReal D = a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2);
  const BoutReal E = a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0);
  const BoutReal F = a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1);
  const BoutReal G = a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1);
  const BoutReal H = a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2);
  const BoutReal I = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);

  // Now construct the output, overwrites input
  const BoutReal detinv = 1.0 / det;

  a(0, 0) = A * detinv;
  a(0, 1) = D * detinv;
  a(0, 2) = G * detinv;
  a(1, 0) = B * detinv;
  a(1, 1) = E * detinv;
  a(1, 2) = H * detinv;
  a(2, 0) = C * detinv;
  a(2, 1) = F * detinv;
  a(2, 2) = I * detinv;

  return std::nullopt;
}
} // namespace bout
