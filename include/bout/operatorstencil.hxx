/*!***********************************************************************
 * \file operatorstencil.hxx 
 * Classes describing the geometry of stencils used for
 * differentiation operators. These can be used to determine how much
 * memory to preallocate when constructing a sparse matrix to
 * represent the operator.
 *
 **************************************************************************
 * Copyright 2019 C. MacMackin
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

#ifndef __OPERATORSTENCIL_H__
#define __OPERATORSTENCIL_H__

#include <vector>
#include <functional>
#include <type_traits>
#include <tuple>
#include <utility>
#include <iterator>
#include <algorithm>

#include <bout/mesh.hxx>
#include <bout/region.hxx>

/// A representation of offsets for indices, which can be added and
/// subtracted from them.
template<class T>
struct IndexOffset {
  static_assert(std::is_same<T, Ind3D>::value || std::is_same<T, Ind2D>::value
		|| std::is_same<T, IndPerp>::value,
		"IndexOffset only works with SpecificInd types");
  int dx = 0, dy = 0, dz = 0;

  const inline IndexOffset xp(int delta_x = 1) const {return {dx + delta_x, dy, dz}; }
  const inline IndexOffset xm(int delta_x = 1) const {return xp(-delta_x); }
  const inline IndexOffset yp(int delta_y = 1) const {return {dx, dy + delta_y, dz}; }
  const inline IndexOffset ym(int delta_y = 1) const {return yp(-delta_y); }
  const inline IndexOffset zp(int delta_z = 1) const {return {dx, dy, dz + delta_z}; }
  const inline IndexOffset zm(int delta_z = 1) const {return zp(-delta_z); }

  IndexOffset &operator+=(const IndexOffset& n) {
    dx += n.dx;
    dy += n.dy;
    dz += n.dz;
    return *this;
  }
  IndexOffset &operator-=(const IndexOffset& n) {
    dx -= n.dx;
    dy -= n.dy;
    dz -= n.dz;
    return *this;
  }
};

template<class T>
inline bool operator==(const IndexOffset<T> &lhs, const IndexOffset<T> &rhs) {
  return lhs.dx == rhs.dx && lhs.dy == rhs.dy && lhs.dz == rhs.dz;
}
template<class T>
inline bool operator!=(const IndexOffset<T> &lhs, const IndexOffset<T> &rhs) {
  return !operator==(lhs, rhs);
}
template<class T>
inline bool operator<(const IndexOffset<T> &lhs, const IndexOffset<T> &rhs) {
  if (lhs.dx != rhs.dx) {
    return lhs.dx < rhs.dx;
  } else if (lhs.dy != rhs.dy) {
    return lhs.dy < rhs.dy;
  } else {
    return lhs.dz < rhs.dz;
  }
}

template<class T>
const inline IndexOffset<T> operator+(IndexOffset<T> lhs, const IndexOffset<T> &rhs) {
  return lhs += rhs;
}
template<class T>
const inline IndexOffset<T> operator-(IndexOffset<T> lhs, const IndexOffset<T> &rhs) {
  return lhs -= rhs;
}

template<class T>
const inline T operator+(const T &lhs, const IndexOffset<T> &rhs) {
  return lhs.offset(rhs.dx, rhs.dy, rhs.dz);
}
template<class T>
const inline T operator+(const IndexOffset<T> &lhs, const T &rhs) {
  return operator+(rhs, lhs);
}
template<class T>
const inline T operator-(const T &lhs, const IndexOffset<T> &rhs) {
  // If CHECKLEVEL >= 3 then SpecificInd<N>.zm() complains about
  // negative values.
  return lhs.offset(-rhs.dx, -rhs.dy, -rhs.dz);
}

using OffsetInd3D = IndexOffset<Ind3D>;
using OffsetInd2D = IndexOffset<Ind2D>;
using OffsetIndPerp = IndexOffset<IndPerp>;

/// A class which can be used to represent the shape of a stencil used
/// to perform some operation. A stencil is made up of pairs of
/// stencil-parts and stencil-tests. A stencil-part is a vector of
/// OffsetIndices. When added to another index, the result is part of
/// the stencil. A stencil-test indicates whether a particular
/// stencil-part should be applied at a given index.
///
/// When trying to get the stencil part for an index, this class will
/// iterate through the part/test pairs in the order which they were
/// added, returning the first stencil-part with a test that
/// passes. If no stencil-part is found with a passing test, then an
/// error is thrown.
template <class T>
class OperatorStencil {
public:
  static_assert(std::is_same<T, Ind3D>::value || std::is_same<T, Ind2D>::value
		|| std::is_same<T, IndPerp>::value,
		"OperatorStencil only works with SpecificInd types");
  using offset = IndexOffset<T>;
  using stencil_part = std::vector<offset>;
  using stencil_test = std::function<bool(T)>;

  struct Stencil {
    stencil_test test;
    stencil_part part;
  };

  /// Add a stencil test/part pair
  void add(stencil_test test, stencil_part stencil) {
    stencils.push_back({test, stencil});
  }

  /// Get the ith stencil-part to have been added
  const stencil_part& getStencilPart(int i) const { return stencils[i].part; }

  /// Get the stencil-part to be used on this index. The method will
  /// iterate through the part/test pairs in the order which they were
  /// added, returning the first stencil-part with a test that
  /// passes. If no stencil-part is found with a passing test, then an
  /// error is thrown.
  const stencil_part& getStencilPart(const T &i) const {
    const auto& result = std::find_if(std::begin(stencils), std::end(stencils),
				      [&i](const auto& stencil) -> bool {
					return stencil.test(i); });
    if (result == std::end(stencils)) {
      throw BoutException("No stencil was specified for element " + toString(i));
    }
    return result->part;
  }

  /// Get the number of elements in the ith stencil-part
  int getStencilSize(const int i) const { return getStencilPart(i).size(); }

  /// Get the number of elements in the stencil part to be used at
  /// this index.
  int getStencilSize(const T i) const { return getStencilPart(i).size(); }

  /// Get the number of stencil-parts to have been added
  int getNumParts() const { return stencils.size(); }

  /// Returns a list of indices for which the stencils contain the
  /// argument
  const std::vector<T> getIndicesWithStencilIncluding(const T &i) const {
    std::vector<T> indices;
    int count = 0;
    for (const auto& item : stencils) {
      for (const auto& j : item.part) {
	T ind = i - j;
	if (getStencilNumber(ind) == count) {
	  indices.push_back(ind);
	}
      }
      count++;
    }
    return indices;
  }

  /// Iterators for the underlying vector data type
  using iterator = typename std::vector<Stencil>::iterator;
  using const_iterator =
      typename std::vector<Stencil>::const_iterator;
  using reverse_iterator =
      typename std::vector<Stencil>::reverse_iterator;
  using const_reverse_iterator =
      typename std::vector<Stencil>::const_reverse_iterator;
  
  iterator begin() { return std::begin(stencils); }
  const_iterator begin() const { return std::begin(stencils); }
  const_iterator cbegin() const { return std::cbegin(stencils); }
  iterator end() { return std::end(stencils); }
  const_iterator end() const { return std::end(stencils); }
  const_iterator cend() const { return std::cend(stencils); }
  reverse_iterator rbegin() { return std::rbegin(stencils); }
  const_reverse_iterator rbegin() const { return std::rbegin(stencils); }
  const_reverse_iterator crbegin() const { return std::crbegin(stencils); }
  reverse_iterator rend() { return std::rend(stencils); }
  const_reverse_iterator rend() const { return std::rend(stencils); }
  const_reverse_iterator crend() const { return std::crend(stencils); }

private:
  std::vector<Stencil> stencils = {};

  /// Returns the position of the first passing stencil test for this
  /// index, or -1 if no test passes.
  int getStencilNumber(const T &i) const {
    auto result = std::find_if(std::begin(stencils), std::end(stencils),
			       [&i](const auto& stencil) { return stencil.test(i); });
    if (result == std::end(stencils)) {
      return -1;
    }
    return std::distance(std::begin(stencils), result);
  }
};  

#endif // __OPERATORSTENCIL_H__
