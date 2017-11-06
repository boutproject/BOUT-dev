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
  SIndices(int i, int nx, int ny, int nz) : i(i), nx(nx), ny(ny), nz(nz) {}
  int i; // index of array
  int nx;
  int ny;
  int nz;

  // This is a little gross, but required so that dereferencing
  // SingleDataIterators works as expected
  SIndices *operator->() { return this; }
  const SIndices *operator->() const { return this; }

  /// Convert to a 2D index
  int i2d() const { return i / nz; }

  /// Convert to x, y, z indices
  int x() const { return (i / nz) / ny; }
  int y() const { return (i / nz) % ny; }
  int z() const { return (i % nz); }

  /*
   * Shortcuts for common offsets, one cell
   * in each direction.
   */

  /// The index one point +1 in x
  SIndices xp() const { return {i + ny * nz, nx, ny, nz}; }
  SIndices xpzp() const {
    return {(i + 1) % nz == 0 ? i + ny * nz - nz + 1 : i + ny * nz + 1, nx, ny, nz};
  }
  SIndices xpzm() const {
    return {i % nz == 0 ? i + ny * nz + nz - 1 : i + ny * nz - 1, nx, ny, nz};
  }
  /// The index one point -1 in x
  SIndices xm() const { return {i - ny * nz, nx, ny, nz}; }
  SIndices xmzp() const {
    return {(i + 1) % nz == 0 ? i - ny * nz - nz + 1 : i - ny * nz + 1, nx, ny, nz};
  }
  SIndices xmzm() const {
    return {i % nz == 0 ? i - ny * nz + nz - 1 : i - ny * nz - 1, nx, ny, nz};
  }
  /// The index one point +1 in y
  SIndices yp() const { return {i + nz, nx, ny, nz}; }
  SIndices ypp() const { return {i + 2 * nz, nx, ny, nz}; }
  /// The index one point -1 in y
  SIndices ym() const { return {i - nz, nx, ny, nz}; }
  SIndices ymm() const { return {i - 2 * nz, nx, ny, nz}; }
  /// The index one point +1 in z. Wraps around zend to zstart
  SIndices zp() const {
    return {(i + 1) % nz == 0 ? i - nz + 1 : i + 1, nx, ny, nz};
  }
  /// The index one point -1 in z. Wraps around zstart to zend
  SIndices zm() const { return {i % nz == 0 ? i + nz - 1 : i - 1, nx, ny, nz}; }

  /*!
   * Add an offset to the index for general stencils
   */
  SIndices offset(int dx, int dy, int dz) const {
    int z0 = i % nz;
    if (dz > 0) {
      int zp = z0;
      for (int j = 0; j < dz; ++j) {
        zp = (zp == nz - 1 ? 0 : zp + 1);
      }
      return {i + ny * nz * dx + nz * dy + zp - z0, nx, ny, nz};
    } else {
      int zm = z0;
      for (; dz != 0; ++dz) {
        zm = (zm == 0 ? nz - 1 : zm - 1);
      }
      return {i + ny * nz * dx + nz * dy + zm - z0, nx, ny, nz};
    }
  }
};

/*!
 * Region for the SingleDataIterator to iterate over
 */
typedef std::vector<int> RegionIndices;

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
class SingleDataIterator {
private:
/*!
 * This initialises OpenMP threads if enabled, and
 * divides iteration index ranges between threads
 */
#ifdef _OPENMP
  void omp_init();
#endif

  SingleDataIterator(); // Disable null constructor

  /// The internal start, end indices of region_iter, needed for OpenMP
  int icountstart, icountend;
  /// Vector of indices this SDI should iterate over
  std::vector<int>::const_iterator region_iter;

public:
  // iterator traits
  using difference_type = int;
  using value_type = SIndices;
  using pointer = const SIndices *;
  using reference = const SIndices &;
  using iterator_category = std::random_access_iterator_tag;

  // SIndexRange needs to be a friend so that it can modify the
  // private members to create the begin and end iterators
  friend class SIndexRange;
  /*!
   * Constructor. This sets index ranges.
   * If OpenMP is enabled, the index range is divided
   * between threads using the omp_init method.
   */
  SingleDataIterator(int nx, int ny, int nz, RegionIndices &region)
      : icountstart(0), icountend(region.size()), region_iter(region.cend()),
        nx(nx), ny(ny), nz(nz), region(region) {
#ifdef _OPENMP
    omp_init();
#endif
  }

  /// Size of the x, y, z dimensions
  const int nx, ny, nz;
  /// Set of indices to iterate over
  const RegionIndices &region;

  /// Pre-increment operator. Use this rather than post-increment when possible
  SingleDataIterator &operator++() {
    ++region_iter;
    return *this;
  }

  /// Post-increment operator
  SingleDataIterator operator++(int) {
    SingleDataIterator tmp(*this);
    ++(*this);
    return tmp;
  }

  /// Pre-decrement operator
  SingleDataIterator &operator--() {
    --region_iter;
    return *this;
  }

  /// Post-decrement operator
  SingleDataIterator operator--(int) {
    SingleDataIterator tmp(*this);
    --(*this);
    return tmp;
  }

  /*!
   * Dereference operators
   */
  SIndices operator*() const { return {*region_iter, nx, ny, nz}; }
  // This is a little gross, but required so that we don't need to do:
  //     (*iter).i
  // in order to access the elements of a SIndices, where iter is a
  // SingleDataIterator
  SIndices operator->() const { return {*region_iter, nx, ny, nz}; }

  /// Arithmetic operators
  SingleDataIterator &operator+=(int n) {
    this->region_iter += n;
    return *this;
  }

  SingleDataIterator &operator-=(int n) { return *this += -n; }

  /// Indexing operator
  SIndices operator[](int n) const { return {*(region_iter + n), nx, ny, nz}; }

  // The following three operators are friend functions so that they
  // have access to the private region_iter

  /// Equivalence operator
  friend bool operator==(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
    return lhs.region_iter == rhs.region_iter;
  }

  /// Less-than operator
  friend bool operator<(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
    return lhs.region_iter < rhs.region_iter;
  }

  /// Iterator difference operator
  friend difference_type operator-(const SingleDataIterator &lhs,
                                   const SingleDataIterator &rhs) {
    return lhs.region_iter - rhs.region_iter;
  }
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
class SIndexRange {
public:
  SIndexRange(int nx, int ny, int nz, RegionIndices &region)
      : nx(nx), ny(ny), nz(nz), region(region) {}

  int nx, ny, nz;
  RegionIndices &region;

  /*!
   * Resets DataIterator to the start of the range
   */
  SingleDataIterator begin() const {
    SingleDataIterator iter{nx, ny, nz, region};
    iter.region_iter = region.cbegin();
    std::advance(iter.region_iter, iter.icountstart);
    return iter;
  }

  /*!
   * Sets DataIterator to one index past the end of the range
   */
  SingleDataIterator end() const {
    SingleDataIterator iter{nx, ny, nz, region};
    iter.region_iter = region.cbegin();
    std::advance(iter.region_iter, iter.icountend);
    return iter;
  }
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

/// Relational operators
inline bool operator!=(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return !operator==(lhs, rhs);
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

/// Arithmetic operators with integers, needed for constant-time random access
inline SingleDataIterator operator+(SingleDataIterator lhs, int n) { return lhs += n; }
inline SingleDataIterator operator+(int n, SingleDataIterator rhs) { return rhs += n; }
inline SingleDataIterator operator-(SingleDataIterator lhs, int n) { return lhs -= n; }

#endif // __SINGLEDATAITERATOR_H__
