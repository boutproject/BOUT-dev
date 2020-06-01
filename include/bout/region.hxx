/**************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Region class by D. Dickinson 2018
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

/// Flexible iterator for Field2D/Field3D
///
/// The `Region` class keeps hold of a set of indices which are used
/// to index a Field. Here, a Field refers to either a Field2D or a
/// Field3D. These indices can either be used directly, or blocks of
/// contiguous indices may be used instead, which allows OpenMP
/// parallelisation.
///
/// The BOUT_FOR helper macro is provided as a replacement
/// for for-loops.
///
/// Separate index classes are available for Field2Ds (Ind2D) and
/// Field3Ds (Ind3D). Note that while an Ind3D can be used to index a
/// Field2D, an Ind2D cannot be used to index a Field3D. This is
/// because an Ind2D essentially doesn't keep track of the
/// z-dimension.

#ifndef __REGION_H__
#define __REGION_H__

#include <algorithm>
#include <ostream>
#include <type_traits>
#include <utility>
#include <vector>

#include "bout_types.hxx"
#include "bout/assert.hxx"
#include "bout/openmpwrap.hxx"

/// The MAXREGIONBLOCKSIZE value can be tuned to try to optimise
/// performance on specific hardware. It determines what the largest
/// contiguous block size can be. As we hope the compiler will vectorise
/// the access to these contiguous blocks, the optimal MAXREGIONBLOCKSIZE
/// is likely related to the vector size etc.
#ifndef MAXREGIONBLOCKSIZE
#define MAXREGIONBLOCKSIZE 64
#endif
#if MAXREGIONBLOCKSIZE <= 0
#error "In region.hxx must set MAXREGIONBLOCKSIZE > 0"
#endif

/// Helper macros for iterating over a Region making use of the
/// contiguous blocks of indices
///
/// @param[in] index  The name of the index variable to use in the loop
/// @param[in] region An already existing Region
///
/// These macros all have the same basic form: an outer loop over
/// blocks of contiguous indices, and an inner loop over the indices
/// themselves. This allows the outer loop to be parallelised with
/// OpenMP while the inner loop can be vectorised with the CPU's
/// native SIMD instructions.
///
/// Alternative forms are also provided for loops that must be done
/// serially, as well as for more control over the OpenMP directives
/// used.
///
/// The different macros:
///
/// - BOUT_FOR: OpenMP-aware version that allows speedup with both
///   OpenMP and native SIMD vectorisation. This should be the
///   preferred form for most loops
/// - BOUT_FOR_SERIAL: for use with inherently serial loops. If BOUT++
///   was not compiled with OpenMP, BOUT_FOR falls back to using this
///   form
/// - BOUT_FOR_INNER: for use on loops inside OpenMP parallel regions,
///   e.g. in order to declare thread private variables outside of the
///   loop
/// - BOUT_FOR_OMP: the most generic form, that takes arbitrary OpenMP
///   directives as an extra argument
///
/// Example
/// -------
/// The following for-loop:
///
///     for (auto index = begin(region); index < end(region); ++index) {
///        A[index] = B[index] + C[index];
///     }
///
/// can be converted to a block region loop like so:
///
///     BOUT_FOR(index, region) {
///        A[index] = B[index] + C[index];
///     }
//


#define BOUT_FOR_SERIAL(index, region)                                                   \
  for (auto block = region.getBlocks().cbegin(), end = region.getBlocks().cend();        \
       block < end; ++block)                                                             \
    for (auto index = block->first; index < block->second; ++index)

#if BOUT_USE_OPENMP 
#define BOUT_FOR_OMP(index, region, omp_pragmas)                                         \
  BOUT_OMP(omp_pragmas)                                                                  \
  for (auto block = region.getBlocks().cbegin(); block < region.getBlocks().cend();      \
       ++block)                                                                          \
    for (auto index = block->first; index < block->second; ++index)
#else
// No OpenMP, so fall back to slightly more efficient serial form
#define BOUT_FOR_OMP(index, region, omp_pragmas) BOUT_FOR_SERIAL(index, region)
#endif

#define BOUT_FOR(index, region)                                                          \
  BOUT_FOR_OMP(index, region, parallel for schedule(OPENMP_SCHEDULE))

#define BOUT_FOR_INNER(index, region)                                                    \
  BOUT_FOR_OMP(index, region, for schedule(OPENMP_SCHEDULE) nowait)


enum class IND_TYPE { IND_3D = 0, IND_2D = 1, IND_PERP = 2 };

/// Indices base class for Fields -- Regions are dereferenced into these
///
/// Provides methods for offsetting by fixed amounts in x, y, z, as
/// well as a generic method for offsetting by any amount in multiple
/// directions.
///
/// Assumes that the offset is less than the grid size in that
/// direction. This assumption is checked for at CHECK=3. This
/// assumption implies that a `FieldPerp` cannot be offset in y, and a
/// `Field2D` cannot be offset in z. A stronger, more expensive check
/// that the resulting offset index doesn't go out of bounds can be
/// enabled at CHECK=4.
///
/// Also provides helper methods for converting Ind2D/Ind3D/IndPerp to x, y, z
/// indices
///
/// Examples
/// --------
///
///     Field3D field, result;
///     auto index = std::begin(region);
///
///     result = field[index->yp()] - field[index->ym()];
template<IND_TYPE N>
class SpecificInd {
public:
  int ind = -1; //< 1D index into Field
private:
  int ny = -1, nz = -1; //< Sizes of y and z dimensions

public:
  SpecificInd() = default;
  SpecificInd(int i, int ny, int nz) : ind(i), ny(ny), nz(nz){};
  explicit SpecificInd(int i) : ind(i) {};

  /// Pre-increment operator
  SpecificInd &operator++() {
    ++ind;
    return *this;
  }

  /// Post-increment operator
  SpecificInd operator++(int) {
    SpecificInd original(*this);
    ++ind;
    return original;
  }

  /// Pre-decrement operator
  SpecificInd &operator--() {
    --ind;
    return *this;
  }

  /// Post-decrement operator
  SpecificInd operator--(int) {
    SpecificInd original(*this);
    --ind;
    return original;
  }

  /// In-place addition
  SpecificInd &operator+=(SpecificInd n) {
    ind += n.ind;
    return *this;
  }

  SpecificInd &operator+=(int n) {
    ind += n;
    return *this;
  }

  /// In-place subtraction
  SpecificInd &operator-=(SpecificInd n) {
    ind -= n.ind;
    return *this;
  }

  SpecificInd &operator-=(int n) {
    ind -= n;
    return *this;
  }

  /// Modulus operator
  SpecificInd operator%(int n) {
    SpecificInd new_ind{ind % n, ny, nz};
    return new_ind;
  }

  /// Convenience functions for converting to (x, y, z)
  int x() const { return (ind / nz) / ny; }
  int y() const { return (ind / nz) % ny; }
  int z() const { return (ind % nz); }

  /// Templated routine to return index.?p(offset), where `?` is one of {x,y,z}
  /// and is determined by the `dir` template argument. The offset corresponds
  /// to the `dd` template argument.
  template<int dd, DIRECTION dir>
  const inline SpecificInd plus() const{
    static_assert(dir == DIRECTION::X || dir == DIRECTION::Y || dir == DIRECTION::Z
                      || dir == DIRECTION::YAligned || dir == DIRECTION::YOrthogonal,
                  "Unhandled DIRECTION in SpecificInd::plus");
    switch(dir) {
    case(DIRECTION::X):
      return xp(dd);
    case(DIRECTION::Y):
    case(DIRECTION::YAligned):
    case(DIRECTION::YOrthogonal):
      return yp(dd);
    case(DIRECTION::Z):
      return zp(dd);
    }
  }

  /// Templated routine to return index.?m(offset), where `?` is one of {x,y,z}
  /// and is determined by the `dir` template argument. The offset corresponds
  /// to the `dd` template argument.
  template<int dd, DIRECTION dir>
  const inline SpecificInd minus() const{
    static_assert(dir == DIRECTION::X || dir == DIRECTION::Y || dir == DIRECTION::Z
                      || dir == DIRECTION::YAligned || dir == DIRECTION::YOrthogonal,
                  "Unhandled DIRECTION in SpecificInd::minus");
    switch(dir) {
    case(DIRECTION::X):
      return xm(dd);
    case(DIRECTION::Y):
    case(DIRECTION::YAligned):
    case(DIRECTION::YOrthogonal):
      return ym(dd);
    case(DIRECTION::Z):
      return zm(dd);
    }
  }

  const inline SpecificInd xp(int dx = 1) const { return {ind + (dx * ny * nz), ny, nz}; }
  /// The index one point -1 in x
  const inline SpecificInd xm(int dx = 1) const { return xp(-dx); }
  /// The index one point +1 in y
  const inline SpecificInd yp(int dy = 1) const {
#if CHECK >= 4
    if (y() + dy < 0 or y() + dy >= ny) {
      throw BoutException("Offset in y ({:d}) would go out of bounds at {:d}", dy, ind);
    }
#endif
    ASSERT3(std::abs(dy) < ny);
    return {ind + (dy * nz), ny, nz};
  }
  /// The index one point -1 in y
  const inline SpecificInd ym(int dy = 1) const { return yp(-dy); }
  /// The index one point +1 in z. Wraps around zend to zstart
  /// An alternative, non-branching calculation is :
  /// ind + dz - nz * ((ind + dz) / nz  - ind / nz)
  /// but this appears no faster (and perhaps slower).  
  const inline SpecificInd zp(int dz = 1) const {
    ASSERT3(dz >= 0);
    dz = dz <= nz ? dz : dz % nz; //Fix in case dz > nz, if not force it to be in range
    return {(ind + dz) % nz < dz ? ind - nz + dz : ind + dz, ny, nz};
  }
  /// The index one point -1 in z. Wraps around zstart to zend
  /// An alternative, non-branching calculation is :
  /// ind - dz + nz * ( (nz + ind) / nz - (nz + ind - dz) / nz)
  /// but this appears no faster (and perhaps slower).
  const inline SpecificInd zm(int dz = 1) const {
    dz = dz <= nz ? dz : dz % nz; //Fix in case dz > nz, if not force it to be in range
    ASSERT3(dz >= 0);
    return {(ind) % nz < dz ? ind + nz - dz : ind - dz, ny, nz};
  }

  // and for 2 cells
  const inline SpecificInd xpp() const { return xp(2); }
  const inline SpecificInd xmm() const { return xm(2); }
  const inline SpecificInd ypp() const { return yp(2); }
  const inline SpecificInd ymm() const { return ym(2); }
  const inline SpecificInd zpp() const { return zp(2); }
  const inline SpecificInd zmm() const { return zm(2); }

  /// Generic offset of \p index in multiple directions simultaneously
  const inline SpecificInd offset(int dx, int dy, int dz) const {
    auto temp = (dz > 0) ? zp(dz) : zm(-dz);
    return temp.yp(dy).xp(dx);
  }
};

/// Relational operators
template<IND_TYPE N>
inline bool operator==(const SpecificInd<N> &lhs, const SpecificInd<N> &rhs) {
  return lhs.ind == rhs.ind;
}

template<IND_TYPE N>
inline bool operator!=(const SpecificInd<N> &lhs, const SpecificInd<N> &rhs) {
  return !operator==(lhs, rhs);
}

template<IND_TYPE N>
inline bool operator<(const SpecificInd<N> &lhs, const SpecificInd<N> &rhs) {
  return lhs.ind < rhs.ind;
}

template<IND_TYPE N>
inline bool operator>(const SpecificInd<N> &lhs, const SpecificInd<N> &rhs) {
  return operator<(rhs, lhs);
}

template<IND_TYPE N>
inline bool operator>=(const SpecificInd<N> &lhs, const SpecificInd<N> &rhs) {
  return !operator<(lhs, rhs);
}

template<IND_TYPE N>
inline bool operator<=(const SpecificInd<N> &lhs, const SpecificInd<N> &rhs) {
  return !operator>(lhs, rhs);
}

/// Arithmetic operators with integers
template<IND_TYPE N>
inline SpecificInd<N> operator+(SpecificInd<N> lhs, const SpecificInd<N> &rhs) { return lhs += rhs; }

template<IND_TYPE N>
inline SpecificInd<N> operator+(SpecificInd<N> lhs, int n) { return lhs += SpecificInd<N>(n); }

template<IND_TYPE N>
inline SpecificInd<N> operator+(int n, SpecificInd<N> rhs) { return rhs += SpecificInd<N>(n); }

template<IND_TYPE N>
inline SpecificInd<N> operator-(SpecificInd<N> lhs, int n) { return lhs -= SpecificInd<N>(n); }

template<IND_TYPE N>
inline SpecificInd<N> operator-(SpecificInd<N> lhs, const SpecificInd<N> &rhs) { return lhs -= rhs; }

/// Define aliases for global indices in 3D and 2D 
using Ind3D = SpecificInd<IND_TYPE::IND_3D>;
using Ind2D = SpecificInd<IND_TYPE::IND_2D>;
using IndPerp = SpecificInd<IND_TYPE::IND_PERP>;

/// Get string representation of Ind3D
inline const std::string toString(const Ind3D& i) {
  return "(" + std::to_string(i.x()) + ", "
             + std::to_string(i.y()) + ", "
             + std::to_string(i.z()) + ")";
}
/// Get string representation of Ind2D
inline const std::string toString(const Ind2D& i) {
  return "(" + std::to_string(i.x()) + ", "
             + std::to_string(i.y()) + ")";
}
/// Get string representation of IndPerp
inline const std::string toString(const IndPerp& i) {
  return "(" + std::to_string(i.x()) + ", "
             + std::to_string(i.z()) + ")";
}

/// Structure to hold various derived "statistics" from a particular region
struct RegionStats {
  int numBlocks = 0;           ///< How many blocks
  int minBlockSize = 0;        ///< Size of smallest block
  int numMinBlocks = 0;        ///< Number of blocks with min size
  int maxBlockSize = 0;        ///< Size of largest block
  int numMaxBlocks = 0;        ///< Number of blocks with max size
  int numSmallBlocks = 0;      ///< Number of "small" blocks, for definition see Region::getStats
  BoutReal maxImbalance = 0;   ///< Ratio of largest block to smallest
};

/// Provide an easy way to report a Region's statistics
inline std::ostream &operator<<(std::ostream &out, const RegionStats &stats){
  if ( stats.numBlocks == 0 ) {
    out << "Empty";
    return out;
  }
  out << "Total blocks : "<< stats.numBlocks;
  out << ", " << "min(count)/max(count) :";
  out << " " << stats.minBlockSize << " (" << stats.numMinBlocks << ")/";
  out << " " << stats.maxBlockSize << " (" << stats.numMaxBlocks << ")";
  out << ", " << "Max imbalance : " << stats.maxImbalance;
  out << ", " << "Small block count : " << stats.numSmallBlocks;
  return out;
}


/// Specifies a set of indices which can be iterated over and begin()
/// and end() methods for range-based for loops.
///
/// `Region` is templated on either `Ind2D` or `Ind3D` for `Field2D`s
/// or `Field3D`s, respectively. Trying to create a `Region` using any
/// other type is a compile time error.
///
/// The set of indices is also broken down into sets of contiguous
/// blocks of at most MAXREGIONBLOCKSIZE indices. This allows loops to
/// be parallelised with OpenMP. Iterating using a "block region" may
/// be more efficient, although it requires a bit more set up. The
/// helper macro BOUT_FOR is provided to simplify things.
///
/// Example
/// -------
///
/// The indices that form a region can be defined manually:
///
///     Region<Ind3D>::RegionIndices indices {0, 2, 4, 8, 3};
///     Region<Ind3D> region(indices);
///
/// then iterated over using begin() and end()
///
///     Field3D f(0.0);
///     for (auto i = region.begin(); i < region.end(); i++ ) {
///       f[i] = 1.0;
///     }
///
/// For the region constructed above the following would display
/// 0, 2, 4, 8, 3
///
///     for (auto i = region.begin(); i < region.end(); i++ ) {
///       output << i.ind << ",";
///     }
///
/// or the more convenient region for loop:
///
///     for (const auto &i : r) {
///       f[i] = a[i] + b[i];
///     }
///
/// For performance the BOUT_FOR macro should
/// allow OpenMP parallelisation and hardware vectorisation.
///
///     BOUT_FOR(i, region) {
///       f[i] = a[i] + b[i];
///     }
///
/// If you wish to vectorise but can't use OpenMP then
/// there is a serial verion of the macro:
///
///     BoutReal max=0.;
///     BOUT_FOR_SERIAL(i, region) {
///       max = f[i] > max ? f[i] : max;
///     }
template <typename T = Ind3D> class Region {
  // Following prevents a Region being created with anything other
  // than Ind2D, Ind3D or IndPerp as template type
  static_assert(std::is_base_of<Ind2D, T>::value || std::is_base_of<Ind3D, T>::value || std::is_base_of<IndPerp, T>::value,
                "Region must be templated with one of IndPerp, Ind2D or Ind3D");

public:
  using data_type = T;

  /// Indices to iterate over
  using RegionIndices = std::vector<T>;
  /// Start and end of contiguous region. This describes a range [block.first,block.second)
  using ContiguousBlock = std::pair<T, T>;
  /// Collection of contiguous regions
  using ContiguousBlocks = std::vector<ContiguousBlock>;

  // Type aliases for STL-container compatibility
  using value_type = T;
  using reference = value_type&;
  using const_reference = const value_type&;
  using size_type = typename RegionIndices::size_type;
  using iterator = typename RegionIndices::iterator;
  using const_iterator = typename RegionIndices::const_iterator;

  // NOTE::
  // Probably want to require a mesh in constructor, both to know nx/ny/nz
  // but also to ensure consistency etc.

  // Want to construct from:
  // * vector of ints
  // * existing region
  // * start/stop of each dimension.

  // Want to make this private to disable but think it may be needed as we put Regions
  // into maps which seems to need to be able to make "empty" objects.
  Region<T>() = default;

  Region<T>(int xstart, int xend, int ystart, int yend, int zstart, int zend, int ny,
            int nz, int maxregionblocksize = MAXREGIONBLOCKSIZE)
      : ny(ny), nz(nz) {
#if CHECK > 1
    if (std::is_base_of<Ind2D, T>::value) {
      if (nz != 1)
        throw BoutException(
            "Trying to make Region<Ind2D> with nz = {:d}, but expected nz = 1", nz);
      if (zstart != 0)
        throw BoutException(
            "Trying to make Region<Ind2D> with zstart = {:d}, but expected zstart = 0",
            zstart);
      if (zstart != 0)
        throw BoutException(
            "Trying to make Region<Ind2D> with zend = {:d}, but expected zend = 0", zend);
    }

    if (std::is_base_of<IndPerp, T>::value) {
      if (ny != 1)
        throw BoutException(
            "Trying to make Region<IndPerp> with ny = {:d}, but expected ny = 1", ny);
      if (ystart != 0)
        throw BoutException(
            "Trying to make Region<IndPerp> with ystart = {:d}, but expected ystart = 0",
            ystart);
      if (ystart != 0)
        throw BoutException(
            "Trying to make Region<IndPerp> with yend = {:d}, but expected yend = 0",
            yend);
    }
#endif
    
    indices = createRegionIndices(xstart, xend, ystart, yend, zstart, zend, ny, nz);
    blocks = getContiguousBlocks(maxregionblocksize);
  };

  Region<T>(RegionIndices &indices, int maxregionblocksize = MAXREGIONBLOCKSIZE) : indices(indices) {
    blocks = getContiguousBlocks(maxregionblocksize);
  };

  Region<T>(ContiguousBlocks &blocks) : blocks(blocks) {
    indices = getRegionIndices();
  };

  /// Destructor
  ~Region() = default;

  /// Expose the iterator over indices for use in range-based
  /// for-loops or with STL algorithms, etc.
  ///
  /// Note that if the indices are altered using these iterators, the
  /// blocks may become out of sync and will need to manually updated
  typename RegionIndices::iterator begin() { return std::begin(indices); };
  typename RegionIndices::const_iterator begin() const { return std::begin(indices); };
  typename RegionIndices::const_iterator cbegin() const { return indices.cbegin(); };
  typename RegionIndices::iterator end() { return std::end(indices); };
  typename RegionIndices::const_iterator end() const { return std::end(indices); };
  typename RegionIndices::const_iterator cend() const { return indices.cend(); };

  const ContiguousBlocks &getBlocks() const { return blocks; };
  const RegionIndices &getIndices() const { return indices; };

  /// Set the indices and ensure blocks updated
  void setIndices (RegionIndices &indicesIn, int maxregionblocksize = MAXREGIONBLOCKSIZE) {
    indices = indicesIn;
    blocks = getContiguousBlocks(maxregionblocksize);
  };

  /// Set the blocks and ensure indices updated
  void setBlocks (ContiguousBlocks &blocksIn) {
    blocks = blocksIn;
    indices = getRegionIndices();
  };

  /// Return a new Region that has the same indices as this one but
  /// ensures the indices are sorted.
  Region<T> asSorted(){
    auto sortedIndices = getIndices();
    std::sort(std::begin(sortedIndices), std::end(sortedIndices));
    return Region<T>(sortedIndices);
  };

  /// Sort this Region in place
  Region<T>& sort() {
    *this = this->asSorted();
    return *this;
  }

  /// Return a new Region that has the same indices as this one but
  /// ensures the indices are sorted and unique (i.e. not duplicate
  /// indices). Note this sorts the input.
  Region<T> asUnique() {
    // As we don't really expect a lot of duplicates this approach should be
    // OK. An alternative is to make a std::set from the indices and then
    // convert back to a vector, but this is typically more expensive.
    auto sortedIndices = getIndices();
    std::sort(std::begin(sortedIndices), std::end(sortedIndices));
    // Get iterator pointing at the end of the unique points
    auto newEnd = std::unique(std::begin(sortedIndices), std::end(sortedIndices));
    // Remove non-unique points from end of vector
    sortedIndices.erase(newEnd, std::end(sortedIndices));
    return Region<T>(sortedIndices);
  }

  /// Make this Region unique in-place
  Region<T>& unique() {
    *this = this->asUnique();
    return *this;
  }

  /// Return a new region equivalent to *this but with indices contained
  /// in mask Region removed
  Region<T> mask(const Region<T> & maskRegion){
    // Get mask indices and sort as we're going to be searching through
    // this vector so if it's sorted we can be more efficient
    auto maskIndices = maskRegion.getIndices();
    std::sort(std::begin(maskIndices), std::end(maskIndices));

    // Get the current set of indices that we're going to mask and then
    // use to create the result region.
    auto currentIndices = getIndices();

    // Lambda that returns true/false depending if the passed value is in maskIndices
    // With C++14 T can be auto instead
    auto isInVector = [&](T val){
      return std::binary_search(std::begin(maskIndices), std::end(maskIndices), val);
    };

    // Erase elements of currentIndices that are in maskIndices
    currentIndices.erase(
             std::remove_if(std::begin(currentIndices), std::end(currentIndices),
                    isInVector),
             std::end(currentIndices)
             );

    // Update indices
    setIndices(currentIndices);

    return *this; // To allow command chaining
  };

  /// Returns a new region including only indices contained in both
  /// this region and the other.
  Region<T> getUnion(const Region<T>& otherRegion) {
    // Get other indices and sort as we're going to be searching through
    // this vector so if it's sorted we can be more efficient
    auto otherIndices = otherRegion.getIndices();
    std::sort(std::begin(otherIndices), std::end(otherIndices));

    // Get the current set of indices that we're going to get the
    // union with and then use to create the result region.
    auto currentIndices = getIndices();

    // Lambda that returns true/false depending if the passed value is in otherIndices
    // With C++14 T can be auto instead
    auto notInVector = [&](T val) {
      return !std::binary_search(std::begin(otherIndices), std::end(otherIndices), val);
    };

    // Erase elements of currentIndices that are in maskIndices
    currentIndices.erase(
        std::remove_if(std::begin(currentIndices), std::end(currentIndices), notInVector),
        std::end(currentIndices));

    // Update indices
    setIndices(currentIndices);

    return *this; // To allow command chaining
  }

  /// Accumulate operator
  Region<T> & operator+=(const Region<T> &rhs){
    (*this) = (*this) + rhs;
    return *this;
  }

  /// Offset all indices by fixed value
  Region<T> &offset(int offset) {
    // Exit early if possible
    if (offset == 0) {
      return *this;
    }

    auto oldInd = getIndices();
    RegionIndices newInd(oldInd.size());

    for (unsigned int i = 0; i < oldInd.size(); i++) {
      newInd[i] = T{oldInd[i].ind + offset, ny, nz};
    }

    setIndices(newInd);

    return *this;
  }

  /// Shift all indices by fixed value but wrap around on a given
  /// period. This is intended to act in a similar way as numpy's roll.
  /// It should be helpful to calculate offset arrays for periodic
  /// directions (e.g. z). For example for shift = 1, period = mesh->LocalNz
  /// we would find the zplus indices. For shift = mesh->LocalNy*mesh->LocalNz,
  /// period = mesh->LocalNx*mesh->LocalNy*mesh->LocalNz we find xplus indices.
  Region<T> & periodicShift(int shift, int period){
    // Exit early if possible
    if ( shift == 0 || period == 1 ){
      return *this;
    }
    // Handle -ve shifts as a +ve shift
    if ( shift < 0 ){
      return periodicShift(period+shift, period);
    }
    auto newInd = getIndices();

    // The calculation of the periodic shifted index is as follows
    //   localPos = index + shift % period;  // Find the shifted position within the period
    //   globalPos = (index/period) * period; // Find which period block we're in
    //   newIndex = globalPos + localPos;
    for (unsigned int i = 0; i < newInd.size(); i++){
      int index = newInd[i].ind;
      int whichBlock = index / period;
      newInd[i].ind = ((index + shift) % period) + period * whichBlock;
    };

    setIndices(newInd);

    return *this;
  }

  /// Number of indices (possibly repeated)
  unsigned int size() const {
    return indices.size();
  }

  /// Returns a RegionStats struct desribing the region
  RegionStats getStats() const {
    RegionStats result;

    result.numBlocks = blocks.size();
    if ( result.numBlocks == 0 ) return result;
    
    std::vector<int> blockSizes(result.numBlocks);
    
    // Get the size of each block using lambda to calculate size
    std::transform(std::begin(blocks), std::end(blocks), std::begin(blockSizes),
		   [](const ContiguousBlock &a) { return a.second.ind - a.first.ind;});

    auto minMaxSize = std::minmax_element(std::begin(blockSizes), std::end(blockSizes));

    result.minBlockSize = *(minMaxSize.first); //Note have to derefence to get actual value
    result.numMinBlocks = std::count(std::begin(blockSizes), std::end(blockSizes), result.minBlockSize);

    result.maxBlockSize = *(minMaxSize.second); //Note have to derefence to get actual value
    result.numMaxBlocks = std::count(std::begin(blockSizes), std::end(blockSizes), result.maxBlockSize);

    
    result.maxImbalance = static_cast<BoutReal>(result.maxBlockSize)/static_cast<BoutReal>(result.minBlockSize);

    // Count the number of small blocks, defined as blocks less than smallSizeFrac of maxBlockSize
    const BoutReal smallSizeFrac = 0.5;
    result.numSmallBlocks =
      std::count_if(std::begin(blockSizes), std::end(blockSizes),
		    [&result, smallSizeFrac](int theSize) {
		      return theSize < smallSizeFrac * result.maxBlockSize;
		    });

    return result;
  }

  // We could sort indices (either forcibly or on request) to try
  // to ensure most contiguous+ordered access. Probably ok in common
  // use but would cause problems if order is important, for example
  // if we could use regions to transpose/reorder data. If regions were
  // sorted this would prevent this usage.

private:
  RegionIndices indices;   //< Flattened indices
  ContiguousBlocks blocks; //< Contiguous sections of flattened indices
  int ny = -1;             //< Size of y dimension
  int nz = -1;             //< Size of z dimension

  /// Helper function to create a RegionIndices, given the start and end
  /// points in x, y, z, and the total y, z lengths
  inline RegionIndices createRegionIndices(int xstart, int xend, int ystart, int yend,
                                           int zstart, int zend, int ny, int nz) {

    if ((xend + 1 <= xstart) ||
        (yend + 1 <= ystart) ||
        (zend + 1 <= zstart)) {
      // Empty region
      return {};
    }

    ASSERT1(ny > 0);
    ASSERT1(nz > 0);

    int len = (xend - xstart + 1) * (yend - ystart + 1) * (zend - zstart + 1);
    ASSERT1(len > 0);
    // Guard against invalid length ranges
    if (len <= 0 ) return {};

    RegionIndices region(len, {-1, ny, nz});

    int x = xstart;
    int y = ystart;
    int z = zstart;

    bool done = false;
    int j = -1;
    while (!done) {
      j++;
      region[j].ind = (x * ny + y) * nz + z;
      if (x == xend && y == yend && z == zend) {
        done = true;
      }
      ++z;
      if (z > zend) {
        z = zstart;
        ++y;
        if (y > yend) {
          y = ystart;
          ++x;
        }
      }
    }
    return region;
  }


  /// Returns a vector of all contiguous blocks contained in the passed region.
  /// Limits the maximum size of any contiguous block to maxBlockSize.
  /// A contiguous block is described by the inclusive start and the exclusive end
  /// of the contiguous block.
  ContiguousBlocks getContiguousBlocks(int maxregionblocksize) const {
    ASSERT1(maxregionblocksize>0);
    const int npoints = indices.size();
    ContiguousBlocks result;
    int index = 0; // Index within vector of indices

    while (index < npoints) {
      const T startIndex = indices[index];
      int count =
          1; // We will always have at least startPair in the block so count starts at 1

      // Consider if the next point should be added to this block
      for (index++; count < maxregionblocksize; index++) {
        if (index >= npoints) {
          break;
        }
        if ((indices[index].ind - indices[index - 1].ind) == 1) {
          count++;
        } else { // Reached the end of this block so break
          break;
        }
      }

      // This contains the inclusive end currently
      T lastIndex = indices[index-1];
      // Increase the index stored by one to get exclusive end
      lastIndex++;
      // Add pair to output, denotes inclusive start and exclusive end
      result.push_back({startIndex, lastIndex});
    }

    return result;
  }


  /// Constructs the vector of indices from the stored blocks information
  RegionIndices getRegionIndices() {
    RegionIndices result;
    // This has to be serial unless we can make result large enough in advance
    // otherwise there will be a race between threads to extend the vector
    BOUT_FOR_SERIAL(curInd, (*this)) {
      result.push_back(curInd);
    }
    return result;
  }
};

/// Return a new region with sorted indices
template<typename T>
Region<T> sort(Region<T> &region) {
  return region.asSorted();
}

/// Return a new region with unique indices
template<typename T>
Region<T> unique(Region<T> &region) {
  return region.asUnique();
}

/// Return a masked version of a region
template<typename T>
Region<T> mask(const Region<T> &region, const Region<T> &mask) {
  auto result = region;
  return result.mask(mask);
}

/// Return the union of two regions
template <typename T>
Region<T> getUnion(const Region<T>& region, const Region<T>& otherRegion) {
  auto result = region;
  return result.getUnion(otherRegion);
}

/// Return a new region with combined indices from two Regions
/// This doesn't attempt to avoid duplicate elements or enforce
/// any sorting etc. but could be done if desired.
/// -
/// Addition is currently simple and just extends. Probably mostly ok 
/// but we could seek to remove duplicate points. Note we do
/// want to allow duplicate points (one reason we use vector and
/// not set) but what if we add a region that has some duplicates?
/// We could retain them but common usage would probably not want
/// the duplicates.
template<typename T>
Region<T> operator+(const Region<T> &lhs, const Region<T> &rhs){
  auto indices = lhs.getIndices(); // Indices is a copy of the indices
  auto indicesRhs = rhs.getIndices();
  indices.insert(std::end(indices), std::begin(indicesRhs), std::end(indicesRhs));
  return Region<T>(indices);
}

/// Returns a new region based on input but with indices offset by
/// a constant
template<typename T>
Region<T> offset(const Region<T> &region, int offset){
  auto result = region;
  return result.offset(offset);
}

/// Return the number of indices in a Region
template<typename T>
unsigned int size(const Region<T> &region){
  return region.size();
}

#endif /* __REGION_H__ */
