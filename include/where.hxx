/*!*************************************************************************
 * \file where.hxx
 * 
 * A set of functions which choose between two values
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

#ifndef __WHERE_H__
#define __WHERE_H__

#include "field3d.hxx"
#include "field2d.hxx"

/// For each point, choose between two inputs based on a third input
///
/// @param[in] test   The value which determines which input to use
/// @param[in] gt0    Uses this value if test > 0.0
/// @param[in] le0    Uses this value if test <= 0.0
template <class T, class U, class V,
          class ResultType = typename std::common_type<T, U, V>::type>
auto where(const T& test, const U& gt0, const V& le0) -> ResultType {
  ASSERT1(areFieldsCompatible(test, gt0));
  ASSERT1(areFieldsCompatible(test, le0));

  ResultType result{emptyFrom(test)};

  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    if (test[i] > 0.0) {
      result[i] = gt0[i];
    } else {
      result[i] = le0[i];
    }
  }
  return result;
}

template <class T, class U,
          class ResultType = typename std::common_type<T, U>::type>
auto where(const T& test, const U& gt0, BoutReal le0) -> ResultType {
  ASSERT1(areFieldsCompatible(test, gt0));

  ResultType result{emptyFrom(test)};

  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    if (test[i] > 0.0) {
      result[i] = gt0[i];
    } else {
      result[i] = le0;
    }
  }
  return result;
}

template <class T, class V,
          class ResultType = typename std::common_type<T, V>::type>
auto where(const T& test, BoutReal gt0, const V& le0) -> ResultType {
  ASSERT1(areFieldsCompatible(test, le0));

  ResultType result{emptyFrom(test)};

  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    if (test[i] > 0.0) {
      result[i] = gt0;
    } else {
      result[i] = le0[i];
    }
  }
  return result;
}

template <class T, class ResultType = T>
auto where(const T& test, BoutReal gt0, BoutReal le0) -> ResultType {

  ResultType result{emptyFrom(test)};

  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    if (test[i] > 0.0) {
      result[i] = gt0;
    } else {
      result[i] = le0;
    }
  }
  return result;
}

#endif // __WHERE_H__
