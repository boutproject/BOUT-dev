/*


 */

#ifndef __DATAITERATOR_H__
#define __DATAITERATOR_H__

#include <iterator>
#include <iostream>
#include "unused.hxx"
#include "assert.hxx"

#ifdef _OPENMP
#include <omp.h>
inline int DI_spread_work(int num_work, int thread, int max_thread);
#endif

/*!
 * Set of indices - DataIterator is dereferenced into these
 */
struct Indices {
  int x;
  int y;
  int z;
};

#define DI_GET_END ((void *) NULL)

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
class DataIterator
  : public std::iterator<std::bidirectional_iterator_tag, Indices> {
private:
  /*!
   * This initialises OpenMP threads if enabled, and
   * divides iteration index ranges between threads
   */
#ifdef _OPENMP
  inline void omp_init(bool end);
#endif
public:
  /*!
   * Constructor. This sets index ranges.
   * If OpenMP is enabled, the index range is divided
   * between threads using the omp_init method.
   */ 
  DataIterator(int xs, int xe,
	       int ys, int ye,
	       int zs, int ze) : 
#ifndef _OPENMP
    x(xs), y(ys), z(zs),
    // start / end : start and end point of the iterator
    xstart(xs),   ystart(ys),   zstart(zs),
    xend(xe),     yend(ye),     zend(ze),
    // Min / Max are the values of the domain we iterate over
    xmin(xstart), ymin(ystart), zmin(zstart),
    xmax(xend),   ymax(yend),   zmax(zend),
#else
    // In the case of OpenMP each processor has a subset.
    // Therefore we initiallise here only the common ones, i.e. the
    // total size of the domain, not what this processor is doing
    xmin(xs),     ymin(ys),     zmin(zs),
    xmax(xe),     ymax(ye),     zmax(ze),
#endif
    isEnd(false)
  {
#ifdef _OPENMP
    omp_init(false);
#endif
  }

  /*!
   * set end();
   * use as DataIterator(int,int,int,int,int,int,DI_GET_END);
   */
  DataIterator(int xs, int xe,
               int ys, int ye,
               int zs, int ze, void* UNUSED(dummy)) :
#ifndef _OPENMP
    x(xe), y(ye), z(ze),
    xstart(xs),   ystart(ys),   zstart(zs),
    xend(xe),     yend(ye),     zend(ze),
    xmin(xstart), ymin(ystart), zmin(zstart),
    xmax(xend),   ymax(yend),   zmax(zend),
#else
    xmin(xs), ymin(ys),   zmin(zs),
    xmax(xe), ymax(ye),   zmax(ze),
#endif
    isEnd(true)
  {
#ifdef _OPENMP
    omp_init(true);
#endif
    next();
  }
  
  /*!
   * The index variables, updated during loop
   * Should make these private and provide getters?
   */
  int x, y, z;

  /// Pre-increment operator. Use this rather than post-increment when possible
  DataIterator& operator++() { next(); return *this; }
  
  /// Post-increment operator
  DataIterator operator++(int) { DataIterator tmp(*this); next(); return tmp; }
  
  /// Pre-decrement operator
  DataIterator& operator--() { prev(); return *this; }
  
  /// Post-decrement operator
  DataIterator operator--(int) { DataIterator tmp(*this); prev(); return tmp; }

  /// Comparison operator. Most common use is in for loops
  inline bool operator!=(const DataIterator& rhs) const {
    if (rhs.isEnd){
      return !this->done();
    } else {
      return  !(x == rhs.x && y == rhs.y && z == rhs.z);
    }
  }
  
  /*!
   * Dereference operators
   * These are needed because the C++11 for loop
   * dereferences the iterator
   */
  DataIterator& operator*() {
    return *this;
  }
  
  /*!
   * Const dereference operator. 
   * Needed because C++11 for loop dereferences the iterator
   */
  const DataIterator& operator*() const {
    return *this;
  }

  /*!
   * Add an offset to the index for general stencils
   */
  const inline Indices offset(int dx, int dy, int dz) const {
    if (dz>0){
      int zp=z;
      for (int j=0;j<dz;++j)
        zp=(zp == zmax ? zmin : zp+1);
      return {x+dx, y+dy, zp };
    } else {
      int zm=z;
      for (;dz!= 0;++dz)
        zm = (zm == zmin ? zmax : zm-1);
      return {x+dx, y+dy, zm };
    }
  }
  
  /*
   * Shortcuts for common offsets, one cell
   * in each direction.
   */
  // general implementation for an abitrary number of cells
  // which defaults to one if used as .xp();
  /// The index one point +1 in x
  const inline Indices xp(int i=1) const { return {x+i, y, z}; }
  /// The index one point -1 in x
  const inline Indices xm(int i=1) const { return {x-i, y, z}; }
  /// The index one point +1 in y
  const inline Indices yp(int i=1) const { return {x, y+i, z}; }
  /// The index one point -1 in y
  const inline Indices ym(int i=1) const { return {x, y-i, z}; }
  /// The index one point +1 in z. Wraps around zend to zstart
  const inline Indices zp(int i=1) const {
    ASSERT3(i >= 0);
    int zp=z;
    for (int j=0;j<i;++j)
      zp=(zp == zmax ? zmin : zp+1);
    return {x, y, zp }; }
  /// The index one point -1 in z. Wraps around zstart to zend
  const inline Indices zm(int i=1) const {
    ASSERT3(i >= 0);
    int zm=z;
    for (;i!= 0;--i)
      zm = (zm == zmin ? zmax : zm-1);
    return {x, y, zm }; }
  // and for 2 cells
  const inline Indices xpp() const { return xp(2); }
  const inline Indices xmm() const { return xm(2); }
  const inline Indices ypp() const { return yp(2); }
  const inline Indices ymm() const { return ym(2); }
  const inline Indices zpp() const { return zp(2); }
  const inline Indices zmm() const { return zm(2); }
  
  /*!
   * Resets DataIterator to the start of the range
   */
  void start() {
    x = xstart; y = ystart; z = zstart;
  }

  /*!
   * Sets DataIterator to one index past the end of the range
   */ 
  void end() {
    x = xend; y = yend; z = zend;
    next();
  }

  /*!
   * Checks if finished looping. Is this more efficient than
   * using the more idiomatic it != DataIterator::end() ?
   */
  bool done() const {
#ifndef _OPENMP
    return (x > xend) || (x < xstart);
#else //_OPENMP
    return (x == xend && y == yend && z > zend)
      || (x == xend && y > yend)
      || x > xend ||
      (x <= xstart && y <= ystart && z < zstart)  ;
#endif //_OPENMP
  }
  
private:
  DataIterator(); // Disable null constructor

  /// start / end : start and end point of THIS iterator
#ifndef _OPENMP
  const int xstart, ystart, zstart;
  const int xend, yend, zend;
#else
  /// start / end : local to THIS processor
  int xstart, ystart, zstart;
  int xend, yend, zend;
#endif
  /// min / max : size of the domain
  /// same for all processors
  int xmin, ymin, zmin;
  int xmax, ymax, zmax;

  const bool isEnd;
  /// Advance to the next index
  void next() {
    ++z;
    if(z > zmax) {
      z = zmin;
      ++y;
      if(y > ymax) {
	y = ymin;
	++x;
      }
    }
  }

  /// Rewind to the previous index
  void prev() {
    --z;
    if(z < zmin) {
      z = zmax;
      --y;
      if(y < ymin) {
	y = ymax;
	--x;
      }
    }
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
 *     for( auto i : f.region(RGN_NOBNDRY) ) {
 *       f[i] = 1.0;
 *     }
 * 
 * where RGN_NOBNDRY specifies a region not including
 * boundary/guard cells. 
 */
struct IndexRange {
  int xstart, xend;
  int ystart, yend;
  int zstart, zend;
  
  const DataIterator begin() const {
    return DataIterator(xstart, xend, 
                        ystart, yend,
                        zstart, zend);
  }
  const DataIterator end() const {
    return DataIterator(xstart, xend, 
			ystart, yend,
			zstart, zend, DI_GET_END);
  }
};

#ifdef _OPENMP
inline int DI_spread_work(int work,int cp,int np){
  int pp=work/np;
  int rest=work%np;
  int result=pp*cp;
  if (rest > cp){
    result +=cp;
  } else {
    result +=rest;
  }
  return result;
};

inline void DataIterator::omp_init(bool end){
  // In the case of OPENMP we need to calculate the range
  int threads=omp_get_num_threads();
  if (threads > 1){
    int ny=ymax-ymin+1;
    int nz=zmax-zmin+1;
    int work  = (xmax-xmin+1)*ny*nz;
    int current_thread = omp_get_thread_num();
    int begin_index = DI_spread_work(work,current_thread,threads);
    int end_index   = DI_spread_work(work,current_thread+1,threads);
    --end_index;
    zend   = (end_index   % nz) + zmin;
    zstart = (begin_index % nz) + zmin;
    end_index   /= nz;
    begin_index /= nz;
    yend   = (end_index   % ny) + ymin;
    ystart = (begin_index % ny) + ymin;
    end_index   /= ny;
    begin_index /= ny;
    xend   = end_index   + xmin;
    xstart = begin_index + xmin;
  } else {
    zstart = zmin;
    zend   = zmax;
    ystart = ymin;
    yend   = ymax;
    xstart = xmin;
    xend   = xmax;
  }
  if (!end){
    x=xstart;
    y=ystart;
    z=zstart;
  } else {
    x=xend;
    y=yend;
    z=zend;
  }
};
#endif

#endif // __DATAITERATOR_H__
