/*


 */

#ifndef __SINGLEDATAITERATOR_H__
#define __SINGLEDATAITERATOR_H__

#include "unused.hxx"
#include <iostream>
#include <iterator>
#include <output.hxx>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
int SDI_spread_work(int num_work, int thread, int max_thread);
#endif

/*!
 * Set of indices - DataIterator is dereferenced into these
 */
struct SIndices {
  int i; // index of array
  int nx;
  int ny;
  int nz;
};

/*!
 * Region for the SingleDataIterator to iterate over
 */
typedef std::vector<int> RegionIndices;

#define SDI_GET_END ((void *)NULL)

/*!
 * Provides range-based iteration over indices.
 * If OpenMP is enabled, then this divides work between threads.
 *
 * This is used mainly to loop over the indices of fields,
 * and provides convenient ways to index
 *
 * Example
 * -------
 *
 * Start,end values for (x,y,z) can be specified directly:
 *
 *     for(d = DataIterator(xs, xe, ys, ye, zs, ze); !d.done(); ++d) {
 *       // print index
 *       output.write("%d,%d,%d\n", d.x, d.y, d.z);
 *       // Index into a Field3D variable 'f'
 *       output.write("Value = %e\n", f[d]);
 *     }
 *
 * Usually DataIterator is used to loop over fields. Field3D::begin()
 * and Field3D::end() return DataIterator objects:
 *
 *     Field3D f(0.0); // Initialise field
 *     for(auto i : f) { // Loop over all indices, including guard cells
 *       f[i] = i.x; // Indexing using DataIterator
 *     }
 *
 */
class SingleDataIterator
    : public std::iterator<std::bidirectional_iterator_tag, SIndices> {
private:
/*!
 * This initialises OpenMP threads if enabled, and
 * divides iteration index ranges between threads
 */
#ifdef _OPENMP
  void omp_init();
#endif
  void idx_to_xyz(int i);

public:
  /*!
   * Constructor. This sets index ranges.
   * If OpenMP is enabled, the index range is divided
   * between threads using the omp_init method.
   */
  SingleDataIterator(int nx, int ny, int nz, RegionIndices &region)
      : icount(0), icountstart(0), icountend(region.size()), region_iter(region.cbegin()),
        nx(nx), ny(ny), nz(nz), region(region), isEnd(false) {
#ifdef _OPENMP
    omp_init();
#endif
  }

  /*!
   * set end();
   * use as DataIterator(int,int,int,int,int,int,SDI_GET_END);
   */
  SingleDataIterator(int nx, int ny, int nz, RegionIndices &region, void *UNUSED(dummy))
      : icount(0), icountstart(0), icountend(region.size()), region_iter(region.cend()),
        nx(nx), ny(ny), nz(nz), region(region), isEnd(true) {
#ifdef _OPENMP
    omp_init();
#endif
    next();
  }

  /*!
   * The index variables, updated during loop
   * Should make these private and provide getters?
   */
  int icount;
  int icountstart, icountend;
  std::vector<int>::const_iterator region_iter;
  int nx, ny, nz;
  const RegionIndices &region;

  /*!
   * Resets DataIterator to the start of the range
   */
  SingleDataIterator begin() {
    region_iter = region.cbegin();
    std::advance(region_iter, icountstart);
    return *this;
  }

  /*!
   * Sets DataIterator to one index past the end of the range
   */
  SingleDataIterator end() {
    region_iter = region.cbegin();
    std::advance(region_iter, icountend);
    // ++region_iter;
    return *this;
  }

  /// Advance to the next index
  void next() { ++region_iter; }
  /// Rewind to the previous index
  void prev() { --region_iter; }

  /// Pre-increment operator. Use this rather than post-increment when possible
  SingleDataIterator &operator++() {
    next();
    return *this;
  }

  /// Post-increment operator
  SingleDataIterator operator++(int) {
    SingleDataIterator tmp(*this);
    next();
    return tmp;
  }

  /// Pre-decrement operator
  SingleDataIterator &operator--() {
    prev();
    return *this;
  }

  /// Post-decrement operator
  SingleDataIterator operator--(int) {
    SingleDataIterator tmp(*this);
    prev();
    return tmp;
  }

  /*!
   * Dereference operators
   * These are needed because the C++11 for loop
   * dereferences the iterator
   */
  SingleDataIterator &operator*() {
    return *this;
  }

  /*!
   * Const dereference operator.
   * Needed because C++11 for loop dereferences the iterator
   */
  // const SingleDataIterator &operator*() const { icount = *region_iter; return *this; }

  /*!
   * Add an offset to the index for general stencils
   */
  const SIndices offset(int dx, int dy, int dz) const {
    int z0 = *region_iter % nz;
    if (dz > 0) {
      int zp = z0;
      for (int j = 0; j < dz; ++j)
        zp = (zp == nz - 1 ? 0 : zp + 1);
      return {*region_iter + ny * nz * dx + nz * dy + zp - z0, nx, ny, nz};
    } else {
      int zm = z0;
      for (; dz != 0; ++dz)
        zm = (zm == 0 ? nz - 1 : zm - 1);
      return {*region_iter + ny * nz * dx + nz * dy + zm - z0, nx, ny, nz};
    }
  }

  /// Convert to a 2D index
  int i2d() const { return *region_iter / nz; }

  /// Convert to x, y, z indices
  int x() const { return (*region_iter / nz) / ny; }
  int y() const { return (*region_iter / nz) % ny; }
  int z() const { return (*region_iter % nz); }

  /*
   * Shortcuts for common offsets, one cell
   * in each direction.
   */

  /// The index one point +1 in x
  const SIndices xp() const { return {*region_iter + ny * nz, nx, ny, nz}; }
  const SIndices xpzp() const {
    return {(*region_iter + 1) % nz == 0 ? *region_iter + ny * nz - nz + 1 : *region_iter + ny * nz + 1, nx,
            ny, nz};
  }
  const SIndices xpzm() const {
    return {*region_iter % nz == 0 ? *region_iter + ny * nz + nz - 1 : *region_iter + ny * nz - 1, nx, ny,
            nz};
  }
  /// The index one point -1 in x
  const SIndices xm() const { return {*region_iter - ny * nz, nx, ny, nz}; }
  const SIndices xmzp() const {
    return {(*region_iter + 1) % nz == 0 ? *region_iter - ny * nz - nz + 1 : *region_iter - ny * nz + 1, nx,
            ny, nz};
  }
  const SIndices xmzm() const {
    return {*region_iter % nz == 0 ? *region_iter - ny * nz + nz - 1 : *region_iter - ny * nz - 1, nx, ny,
            nz};
  }
  /// The index one point +1 in y
  const SIndices yp() const { return {*region_iter + nz, nx, ny, nz}; }
  const SIndices ypp() const { return {*region_iter + 2 * nz, nx, ny, nz}; }
  /// The index one point -1 in y
  const SIndices ym() const { return {*region_iter - nz, nx, ny, nz}; }
  const SIndices ymm() const { return {*region_iter - 2 * nz, nx, ny, nz}; }
  /// The index one point +1 in z. Wraps around zend to zstart
  const SIndices zp() const {
    return {(*region_iter + 1) % nz == 0 ? *region_iter - nz + 1 : *region_iter + 1, nx, ny, nz};
  }
  /// The index one point -1 in z. Wraps around zstart to zend
  const SIndices zm() const {
    return {*region_iter % nz == 0 ? *region_iter + nz - 1 : *region_iter - 1, nx, ny, nz};
  }

private:
  SingleDataIterator(); // Disable null constructor

  const bool isEnd;
};

/*!
 * Specifies a range of indices which can be iterated over
 * and begin() and end() methods for range-based for loops
 *
 * Example
 * -------
 *
 * Index ranges can be defined manually:
 *
 *     IndexRange r(0, 10, 0, 20, 0, 30);
 *
 * then iterated over using begin() and end()
 *
 *     for( DataIterator i = r.begin(); i != r.end(); i++ ) {
 *       output.write("%d,%d,%d\n", i.x, i.y, i.z);
 *     }
 *
 * or the more convenient range for loop:
 *
 *     for( auto i : r ) {
 *       output.write("%d,%d,%d\n", i.x, i.y, i.z);
 *     }
 *
 * A common use for this class is to loop over
 * regions of a field:
 *
 *     Field3D f(0.0);
 *     for( auto i : f.region(REGION_NOBNDRY) ) {
 *       f[i] = 1.0;
 *     }
 *
 * where REGION_NOBNDRY specifies a region not including
 * boundary/guard cells.
 */
struct SIndexRange {
  SIndexRange(int nx, int ny, int nz, RegionIndices &region)
      : nx(nx), ny(ny), nz(nz), region(region) {}

  int nx, ny, nz;
  RegionIndices &region;

  const SingleDataIterator begin() const {
    return SingleDataIterator(nx, ny, nz, region);
  }
  const SingleDataIterator end() const {
    return SingleDataIterator(nx, ny, nz, region, SDI_GET_END);
  }
};

inline void SingleDataIterator::idx_to_xyz(int i) {
  // function for debugging
  // print x,y,z for a given icount

  // i = (x*ny+y)*nz+z
  output << "i = " << i << ", x = " << ((i / nz) / ny) << ", y = " << (i / nz) % ny
         << ", z = " << (i % nz) << "\n";
};

#ifdef _OPENMP
inline int SDI_spread_work(int work, int cp, int np) {
  // Spread work between threads. If number of points do not
  // spread evenly between threads, put the remaining "rest"
  // points on the threads 0 ... rest-1.
  int pp = work / np;
  int rest = work % np;
  int result = pp * cp;
  if (rest > cp) {
    result += cp;
  } else {
    result += rest;
  }
  return result;
};

inline void SingleDataIterator::omp_init() {
  // In the case of OPENMP we need to calculate the range
  int threads = omp_get_num_threads();
  if (threads > 1) {
    int work = region.size();
    int current_thread = omp_get_thread_num();
    icountstart = SDI_spread_work(work, current_thread, threads);
    icountend = SDI_spread_work(work, current_thread + 1, threads);
  }
};
#endif

/// Comparison operator
inline bool operator==(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return lhs.region_iter == rhs.region_iter;
}

inline bool operator!=(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return !operator==(lhs, rhs);
}

inline bool operator<(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return lhs.region_iter < rhs.region_iter;
}

inline bool operator>(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return operator<(rhs, lhs);
}

inline bool operator<=(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return !operator>(lhs, rhs);
}

inline bool operator>=(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return !operator<(lhs, rhs);
}

#endif // __SINGLEDATAITERATOR_H__
