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
/// The BLOCK_REGION_LOOP helper macro is provided as a replacement
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
#include <type_traits>
#include <utility>
#include <vector>

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
/// @param[in] region An already existing Region
/// @param[in] index  The name of the index variable to use in the loop
/// The rest of the arguments form the loop body.
///
/// The following loops vectorise well when index is type int, needs testing
/// with current approach where index is Ind2D or ind3D
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
///     BLOCK_REGION_LOOP(region, index,
///        A[index] = B[index] + C[index];
///     )
#define BLOCK_REGION_LOOP_SERIAL(region, index, ...)                                     \
  {                                                                                      \
    const auto blocks = region.getBlocks();                                              \
    for (auto block = blocks.begin(); block < blocks.end(); ++block) {                   \
      for (auto index = block->first; index < block->second; ++index) {                  \
        __VA_ARGS__                                                                      \
      }                                                                                  \
    }                                                                                    \
  }

#define BLOCK_REGION_LOOP(region, index, ...)                                            \
  {                                                                                      \
    const auto blocks = region.getBlocks();                                              \
  BOUT_OMP(parallel for)                                                                 \
    for (auto block = blocks.begin(); block < blocks.end(); ++block) {                   \
      for (auto index = block->first; index < block->second; ++index) {                  \
        __VA_ARGS__                                                                      \
      }                                                                                  \
    }                                                                                    \
  }

/// Indices base class for Fields -- Regions are dereferenced into these
class SpecificInd {
public:
  int ind; //< 1D index into Field
  SpecificInd() : ind(-1){};
  SpecificInd(int i) : ind(i){};

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

  /// In-place subtraction
  SpecificInd &operator-=(SpecificInd n) {
    ind -= n.ind;
    return *this;
  }

  /// Modulus operator
  SpecificInd operator%(int n) {
    SpecificInd new_ind(ind % n);
    return new_ind;
  }
};

/// Relational operators
inline bool operator==(const SpecificInd &lhs, const SpecificInd &rhs) {
  return lhs.ind == rhs.ind;
}
inline bool operator!=(const SpecificInd &lhs, const SpecificInd &rhs) {
  return !operator==(lhs, rhs);
}
inline bool operator<(const SpecificInd &lhs, const SpecificInd &rhs) {
  return lhs.ind < rhs.ind;
}
inline bool operator>(const SpecificInd &lhs, const SpecificInd &rhs) {
  return operator<(rhs, lhs);
}
inline bool operator>=(const SpecificInd &lhs, const SpecificInd &rhs) {
  return !operator<(lhs, rhs);
}
inline bool operator<=(const SpecificInd &lhs, const SpecificInd &rhs) {
  return !operator>(lhs, rhs);
}

/// Index-type for `Field3D`s
class Ind3D : public SpecificInd {
public:
  Ind3D() : SpecificInd(){};
  Ind3D(int i) : SpecificInd(i){};
  Ind3D(SpecificInd baseIn) : SpecificInd(baseIn){};

  // Note operator= from base class is always hidden
  // by implicit method so have to be explicit
  Ind3D &operator=(int i) {
    ind = i;
    return *this;
  }

  Ind3D &operator=(SpecificInd i) {
    ind = i.ind;
    return *this;
  }
};

/// Arithmetic operators with integers
inline Ind3D operator+(Ind3D lhs, const Ind3D &rhs) { return lhs += rhs; }
inline Ind3D operator+(Ind3D lhs, int n) { return lhs += n; }
inline Ind3D operator+(int n, Ind3D rhs) { return rhs += n; }
inline Ind3D operator-(Ind3D lhs, int n) { return lhs -= n; }
inline Ind3D operator-(Ind3D lhs, const Ind3D &rhs) { return lhs -= rhs; }

/// Index-type for `Field2D`s
class Ind2D : public SpecificInd {
public:
  Ind2D() : SpecificInd(){};
  Ind2D(int i) : SpecificInd(i){};
  Ind2D(SpecificInd baseIn) : SpecificInd(baseIn){};

  Ind2D &operator=(int i) {
    ind = i;
    return *this;
  }

  Ind2D &operator=(SpecificInd i) {
    ind = i.ind;
    return *this;
  }
};

/// Arithmetic operators with integers
inline Ind2D operator+(Ind2D lhs, const Ind2D &rhs) { return lhs += rhs; }
inline Ind2D operator+(Ind2D lhs, int n) { return lhs += n; }
inline Ind2D operator+(int n, Ind2D rhs) { return rhs += n; }
inline Ind2D operator-(Ind2D lhs, int n) { return lhs -= n; }
inline Ind2D operator-(Ind2D lhs, const Ind2D &rhs) { return lhs -= rhs; }

/// Structure to hold various derived "statistics" from a particular region
struct RegionStats {
  int numBlocks = 0;           // How many blocks
  int minBlockSize = 0;        // Size of smallest block
  int numMinBlocks = 0;        // Number of blocks with min size
  int maxBlockSize = 0;        // Size of largest block
  int numMaxBlocks = 0;        // Number of blocks with max size
  int numSmallBlocks = 0;      // Number of "small" blocks, for definition see Region::getStats
  BoutReal maxImbalance = 0;   // Ratio of largest block to smallest
};

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
/// helper macro BLOCK_REGION_LOOP is provided to simplify things.
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
///     for (auto &i : r) {
///       f[i] = a[i] + b[i];
///     }
///
/// For performance the BLOCK_REGION_LOOP macro should
/// allow OpenMP parallelisation and hardware vectorisation.
///
///     BLOCK_REGION_LOOP(region, i,
///       f[i] = a[i] + b[i];
///     );
///
/// If you wish to vectorise but can't use OpenMP then
/// there is a serial verion of the macro:
///
///     BoutReal max=0.;
///     BLOCK_REGION_LOOP_SERIAL(region, i,
///       max = f[i] > max ? f[i] : max;
///     );
template <typename T = Ind3D> class Region {
  // Following prevents a Region being created with anything other
  // than Ind2D or Ind3D as template type
  static_assert(std::is_base_of<Ind2D, T>::value || std::is_base_of<Ind3D, T>::value,
                "Region must be templated with either Ind2D or Ind3D");

public:
  typedef T data_type;

  /// Indices to iterate over
  typedef std::vector<T> RegionIndices;
  /// Start and end of contiguous region. This describes a range [block.first,block.second)
  typedef std::pair<T, T> ContiguousBlock;
  /// Collection of contiguous regions
  typedef std::vector<ContiguousBlock> ContiguousBlocks;

  // NOTE::
  // Probably want to require a mesh in constructor, both to know nx/ny/nz
  // but also to ensure consistency etc.

  // Want to construct from:
  // * vector of ints
  // * existing region
  // * start/stop of each dimension.

  // Want to make this private to disable but think it may be needed as we put Regions
  // into maps which seems to need to be able to make "empty" objects.
  Region<T>(){};

  Region<T>(int xstart, int xend, int ystart, int yend, int zstart, int zend, int ny,
            int nz, int maxregionblocksize = MAXREGIONBLOCKSIZE) {
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
  ~Region(){};

  /// Expose the iterator over indices for use in range-based
  /// for-loops or with STL algorithms, etc.
  ///
  /// Note that if the indices are altered using these iterators, the
  /// blocks may become out of sync and will need to manually updated
  typename RegionIndices::iterator begin() { return std::begin(indices); };
  typename RegionIndices::const_iterator cbegin() const { return indices.cbegin(); };
  typename RegionIndices::iterator end() { return std::end(indices); };
  typename RegionIndices::const_iterator cend() const { return indices.cend(); };

  ContiguousBlocks getBlocks() const { return blocks; };
  RegionIndices getIndices() const { return indices; };

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
  Region<T> sort(){
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
  Region<T> unique(){
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
      newInd[i] = T(oldInd[i].ind + offset);
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
      newInd[i] = ((index + shift) % period) + period * whichBlock;
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
    result.numSmallBlocks = std::count_if(std::begin(blockSizes), std::end(blockSizes), [&blockSizes, &result, smallSizeFrac](int theSize) { return theSize < smallSizeFrac * result.maxBlockSize;}) ;

    return result;
  }

  // TODO: Should be able to add regions (would just require extending
  // indices and recalculating blocks). This raises question of should
  // we be able to subtract regions, and if so what does that mean.
  // Addition could be simple and just extend or we could seek to
  // remove duplicate points. Former probably mostly ok. Note we do
  // want to allow duplicate points (one reason we use vector and
  // not set) but what if we add a region that has some duplicates?
  // We could retain them but common usage would probably not want
  // the duplicates.

  // We could sort indices (either forcibly or on request) to try
  // to ensure most contiguous+ordered access. Probably ok in common
  // use but would cause problems if order is important, for example
  // if we could use regions to transpose/reorder data. If regions were
  // sorted this would prevent this usage.

private:
  RegionIndices indices;   //< Flattened indices
  ContiguousBlocks blocks; //< Contiguous sections of flattened indices

  /// Helper function to create a RegionIndices, given the start and end
  /// points in x, y, z, and the total y, z lengths
  inline RegionIndices createRegionIndices(int xstart, int xend, int ystart, int yend,
                                           int zstart, int zend, int ny, int nz) {

    if ( (xend + 1 <= xstart) ||
         (yend + 1 <= ystart) ||
         (zend + 1 <= zstart) ) {
      // Empty region
      return {};
    }
    
    ASSERT1(ny > 0);
    ASSERT1(nz > 0);

    int len = (xend - xstart + 1) * (yend - ystart + 1) * (zend - zstart + 1);
    ASSERT1(len > 0);
    RegionIndices region;
    // Guard against invalid length ranges
    if (len <= 0 ) return region;
    region.resize(len);

    int x = xstart;
    int y = ystart;
    int z = zstart;

    bool done = false;
    int j = -1;
    while (!done) {
      j++;
      region[j] = (x * ny + y) * nz + z;
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
    BLOCK_REGION_LOOP_SERIAL((*this), curInd, result.push_back(curInd););
    return result;
  }
};

/// Return a new region with sorted indices
template<typename T>
Region<T> sort(Region<T> &region) {
  return region.asSorted();
};

/// Return a new region with unique indices
template<typename T>
Region<T> unique(Region<T> &region) {
  return region.asUnique();
};

/// Return a masked version of a region
template<typename T>
Region<T> mask(const Region<T> &region, const Region<T> &mask) {
  auto result = region;
  return result.mask(mask);
};

/// Return a new region with combined indices from two Regions
/// This doesn't attempt to avoid duplicate elements or enforce
/// any sorting etc. but could be done if desired.
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
