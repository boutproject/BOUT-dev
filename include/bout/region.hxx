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

template<typename T>
class Region;

#ifndef __REGION_H__
#define __REGION_H__

#include <type_traits>
#include <vector>
#include <utility>

#include "bout/array.hxx"

#ifndef MAXREGIONBLOCKSIZE
#define MAXREGIONBLOCKSIZE 64
#endif

//The following loops assume blocks are made up of ints but regions are now
//ind3D or ind2D to allow care to be taken in indexing fields
//
//Currently set up to only use regions stored in a region map, however
//may want to just pass in a Region instance instead?

#define BLOCK_REGION_LOOP_SERIAL(regionStr, index, ...)			\
  {const auto blocks = mesh->getRegion(regionStr).blocks;		\
  for(auto block = blocks.begin(); block < blocks.end() ; ++block){ \
  for(auto index = (*block).first ; index <= (*block).second; ++index){ \
  __VA_ARGS__ \
    }}}

#define BLOCK_REGION_LOOP(regionStr, index, ...)			\
  {const auto blocks = mesh->getRegion(regionStr).blocks;		\
  BOUT_OMP(parallel for) \
  for(auto block = blocks.begin(); block < blocks.end() ; ++block){ \
  for(auto index = (*block).first ; index <= (*block).second; ++index){ \
  __VA_ARGS__ \
    }}}


class specificInd{
public:
  int ind;
  specificInd& operator++(){++ind;return *this;}
  specificInd& operator++(int){++ind;return *this;}
  bool operator==(const specificInd &x) const {
    return ind == x.ind;
  }
  bool operator!=(const specificInd &x) const {
    return ind != x.ind;
  }
  bool operator>(const specificInd &x) const {
    return ind > x.ind;
  }
  bool operator>=(const specificInd &x) const {
    return ind >= x.ind;
  }
  bool operator<(const specificInd &x) const {
    return ind < x.ind;
  }
  bool operator<=(const specificInd &x) const {
    return ind <= x.ind;
  }
};
//Make identical subclasses with different names
class ind3D: public specificInd{};
class ind2D: public specificInd{};

//Simple structs to act as "strict typedef" allowing us to overload
//on ind type and protect against indexing a field with the wrong type
//of index.
//Need to define a ++() operator 
// struct ind3D{
//   int ind;
//   ind3D operator++(int)
// };

// struct ind2D{
//   int ind;
// };


template<typename T>
class Region {
  //Following prevents a Region being created with anything other
  //than ind2D or ind3D as template type
  static_assert(std::is_base_of<ind2D,T>::value || std::is_base_of<ind3D,T>::value, "Region must be templated with either ind2D or ind3D");

  
 public:
  typedef T data_type;

  /// Indices to iterate over
  typedef std::vector<T> regionIndices;
  /// Start and end of contiguous region
  typedef std::pair<T, T> contiguousBlock;
  /// Collection of contiguous regions
  typedef std::vector<contiguousBlock> contiguousBlocks;

  //NOTE::
  // Probably want to require a mesh in constructor, both to know nx/ny/nz o
  // but also to ensure consistency etc.
  
  //Want to construct from:
  // * vector of ints
  // * existing region
  // * start/stop of each dimension.

  //Want to make this private to disable but think it may be needed as we put Regions into
  //maps which seems to need to be able to make "empty" objects.
  Region(){};
  
  Region(int xstart, int xend, int ystart, int yend,
	 int zstart, int zend, int ny, int nz);
  Region(regionIndices& indices);

  /// Destructor
  ~Region(){};
private:
  regionIndices indices;  //Flattened indices
  contiguousBlocks blocks;//Contiguous sections of flattened indices

  /// Helper function to create a regionIndices, given the start and end
  /// points in x, y, z, and the total y, z lengths
  inline regionIndices createRegionIndices(int xstart, int xend, int ystart, int yend,
  					   int zstart, int zend, int ny, int nz);
    
  /* 
   * Returns a vector of all contiguous blocks contained in the passed region.
   * Limits the maximum size of any contiguous block to maxBlockSize.
   * A contiguous block is described by the inclusive start and the exclusive end
   * of the contiguous block.
   */
  contiguousBlocks getContiguousBlocks() const;
};

#endif /* __REGION_H__ */
