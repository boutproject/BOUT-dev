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

#ifndef __REGION_H__
#define __REGION_H__

#include <type_traits>
#include <utility>
#include <vector>

#include "bout/openmpwrap.hxx"

/*
 * The MAXREGIONBLOCKSIZE value can be tuned to try to optimise
 * performance on specific hardware. It determines what the largest
 * contiguous block size can be. As we hope the compiler will vectorise
 * the access to these contiguous blocks, the optimal MAXREGIONBLOCKSIZE
 * is likely related to the vector size etc.
 */
#ifndef MAXREGIONBLOCKSIZE
#define MAXREGIONBLOCKSIZE 64
#endif
#if MAXREGIONBLOCKSIZE <= 0
#error "In region.hxx must set MAXREGIONBLOCKSIZE > 0"
#endif

// The following loops vectorise well well index is type int, needs testing
// with current approach where index is Ind2D or ind3D
#define BLOCK_REGION_LOOP_SERIAL(region, index, ...)                                     \
  {                                                                                      \
    const auto blocks = region.getBlocks();                                              \
    for (auto block = blocks.begin(); block < blocks.end(); ++block) {                   \
      for (auto index = block->first; index < block->second; ++index) {                 \
        __VA_ARGS__                                                                      \
      }                                                                                  \
    }                                                                                    \
  }

#define BLOCK_REGION_LOOP(region, index, ...)                                            \
  {                                                                                      \
    const auto blocks = region.getBlocks();                                              \
  BOUT_OMP(parallel for)                                                                 \
    for (auto block = blocks.begin(); block < blocks.end(); ++block) {                   \
      for (auto index = block->first; index < block->second; ++index) {                 \
        __VA_ARGS__                                                                      \
      }                                                                                  \
    }                                                                                    \
  }

class SpecificInd {
public:
  int ind;
  SpecificInd() : ind(-1){};
  SpecificInd(int i) : ind(i){};
  SpecificInd &operator++() {
    ++ind;
    return *this;
  }
  SpecificInd operator++(int) {
    SpecificInd original(*this);
    ++ind;
    return original;
  }
  bool operator==(const SpecificInd &x) const { return ind == x.ind; }
  bool operator!=(const SpecificInd &x) const { return ind != x.ind; }
  bool operator>(const SpecificInd &x) const { return ind > x.ind; }
  bool operator>=(const SpecificInd &x) const { return ind >= x.ind; }
  bool operator<(const SpecificInd &x) const { return ind < x.ind; }
  bool operator<=(const SpecificInd &x) const { return ind <= x.ind; }
};
// Make identical subclasses with different names
class Ind3D : public SpecificInd {
public:
  Ind3D() : SpecificInd(){};
  Ind3D(int i) : SpecificInd(i){};
  // Note operator= from base class is always hidden
  // by implicit method so have to be explicit
  Ind3D &operator=(int i) {
    ind = i;
    return *this;
  }
};
class Ind2D : public SpecificInd {
public:
  Ind2D() : SpecificInd(){};
  Ind2D(int i) : SpecificInd(i){};
  Ind2D &operator=(int i) {
    ind = i;
    return *this;
  }
};

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
  // Probably want to require a mesh in constructor, both to know nx/ny/nz o
  // but also to ensure consistency etc.

  // Want to construct from:
  // * vector of ints
  // * existing region
  // * start/stop of each dimension.

  // Want to make this private to disable but think it may be needed as we put Regions
  // into
  // maps which seems to need to be able to make "empty" objects.
  Region<T>(){};

  Region<T>(int xstart, int xend, int ystart, int yend, int zstart, int zend, int ny,
            int nz) {
    indices = createRegionIndices(xstart, xend, ystart, yend, zstart, zend, ny, nz);
    blocks = getContiguousBlocks();
  };

  Region<T>(RegionIndices &indices) : indices(indices) {
    blocks = getContiguousBlocks();
  };

  /// Destructor
  ~Region(){};

  ContiguousBlocks getBlocks() { return blocks; }; // Setter not appropriate
  RegionIndices getIndices() { return indices; };  // Setter could be ok

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
  RegionIndices indices;   // Flattened indices
  ContiguousBlocks blocks; // Contiguous sections of flattened indices

  /// Helper function to create a RegionIndices, given the start and end
  /// points in x, y, z, and the total y, z lengths
  inline RegionIndices createRegionIndices(int xstart, int xend, int ystart, int yend,
                                           int zstart, int zend, int ny, int nz) {

    ASSERT1(xend + 1 > xstart);
    ASSERT1(yend + 1 > ystart);
    ASSERT1(zend + 1 > zstart);
    ASSERT1(ny > 0);
    ASSERT1(nz > 0);

    int len = (xend - xstart + 1) * (yend - ystart + 1) * (zend - zstart + 1);
    RegionIndices region(len);
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

  /*
   * Returns a vector of all contiguous blocks contained in the passed region.
   * Limits the maximum size of any contiguous block to maxBlockSize.
   * A contiguous block is described by the inclusive start and the exclusive end
   * of the contiguous block.
   */
  ContiguousBlocks getContiguousBlocks() const {
    const int npoints = indices.size();
    ContiguousBlocks result;
    int index = 0; // Index within vector of indices

    while (index < npoints) {
      const T startIndex = indices[index];
      int count =
          1; // We will always have at least startPair in the block so count starts at 1

      // Consider if the next point should be added to this block
      for (index++; count < MAXREGIONBLOCKSIZE; index++) {
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
};

#endif /* __REGION_H__ */
