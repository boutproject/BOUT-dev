/*


 */

#ifndef __DATAITERATOR_H__
#define __DATAITERATOR_H__

#include <iterator>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
int DI_spread_work(int num_work, int thread, int max_thread);
#endif

// Set of indices - DataIterator is dereferenced into these
struct Indices {
  int x;
  int y;
  int z;
};

#define DI_GET_END ((void *) NULL)

class DataIterator
  : public std::iterator<std::bidirectional_iterator_tag, Indices> {
private:
  void omp_init(int xs,int xe,bool end);
public:
  /// Constructor. This would set ranges. Could depend on thread number
  DataIterator(int xs, int xe,
	       int ys, int ye,
	       int zs, int ze) : 
#ifndef _OPENMP
    x(xs), y(ys), z(zs),
    xstart(xs),   ystart(ys),   zstart(zs),
    xmin(xstart), ymin(ystart), zmin(zstart),
    xend(xe),     yend(ye),     zend(ze),
    xmax(xend),   ymax(yend),   zmax(zend),
#else
    xmin(xstart), ymin(ys),     zmin(zs),
    xmax(xend),   ymax(ye),     zmax(ze),
#endif
    isEnd(false)
  {
    omp_init(xs,xe,false);
  }

  // set end();
  // use as DataIterator(int,int,int,int,int,int,DI_GET_END);
  DataIterator(int xs, int xe,
	       int ys, int ye,
	       int zs, int ze,void* dummy) :
#ifndef _OPENMP
    x(xe), y(ye), z(ze),
    xstart(xs),   ystart(ys),   zstart(zs),
    xmin(xstart), ymin(ystart), zmin(zstart),
    xend(xe),     yend(ye),     zend(ze),
    xmax(xend),   ymax(yend),   zmax(zend),
#else
    xmin(xstart), ymin(ys),   zmin(zs),
    xmax(xend),   ymax(ye),   zmax(ze),
#endif
    isEnd(true)
  {
    omp_init(xs,xe,true);
    next();
  }
  /// The index variables, updated during loop
  // Should make these private and provide getters?
  int x, y, z;

  /// Increment operators
  DataIterator& operator++() { next(); return *this; }
  DataIterator operator++(int) { DataIterator tmp(*this); next(); return tmp; }
  DataIterator& operator--() { prev(); return *this; }
  DataIterator operator--(int) { DataIterator tmp(*this); prev(); return tmp; }

  // Comparison operator
  inline bool operator!=(const DataIterator& rhs) const {
    //return  !(x == rhs.x && y == rhs.y && z == rhs.z);
    if (rhs.isEnd){
      return !this->done();
    } else {
      return  !(x == rhs.x && y == rhs.y && z == rhs.z);
    }
  }
  
  /*
   * Dereference operators
   * These are needed because the C++11 for loop
   * dereferences the iterator
   */
  DataIterator& operator*() {
    return *this;
  }
  const DataIterator& operator*() const {
    return *this;
  }

  /*
   * Add an offset to the index for general stencils
   */
  const Indices offset(int dx, int dy, int dz) const {
    int nz = zend-zstart+1;
    return {x+dx, y+dy, ((z+dz-zstart + nz) % nz) + zstart};
  }
  /*
   * Shortcuts for common offsets, one cell
   * in each direction.
   */
  const Indices xp() const { return {x+1, y, z}; }
  const Indices xm() const { return {x-1, y, z}; }
  const Indices yp() const { return {x, y+1, z}; }
  const Indices ym() const { return {x, y-1, z}; }
  // Z indices should wrap in zstart, zend range (?)
  const Indices zp() const { return {x, y, z == zend ? zstart : z+1}; }
  const Indices zm() const { return {x, y, z == zstart ? zend : z-1}; }
  
  void start() {
    x = xstart; y = ystart; z = zstart;
  }

  void end() {
    x = xend; y = yend; z = zend;
    next();
  }

  /// Checks if finished looping. Is this more efficient than
  /// using the more idiomatic it != MeshIterator::end() ?
  bool done() const {
#ifndef _OPENMP
    return (x > xend) || (x < xstart);
#else //_OPENMP
    return (x == xend && y == yend && z > zend)
      || x > xend ||
      (x <= xstart && y <= ystart && z < zstart)  ;
#endif //_OPENMP
  }
  
private:
  DataIterator(); // Disable null constructor

  int xstart, xend;
  int ystart, yend;
  int zstart, zend;
#ifndef _OPENMP
  int &xmin, &ymin, &zmin;
  int &xmax, &ymax, &zmax;
#else
  int xmin, ymin, zmin;
  int xmax, ymax, zmax;
#endif
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

/*
 * Specifies a range of indices which can be iterated over
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

inline void DataIterator::omp_init(int xs, int xe,bool end){
  xmin=xs;
  xmax=xe;
  // In the case of OPENMP we need to calculate the range
  int threads=omp_get_num_threads();
  if (threads > 1){
    int ny=ymax-ymin+1;
    int nz=zmax-zmin+1;
    int work  = (xmax-xmin+1)*ny*nz;
    int current_thread = omp_get_thread_num();
    int begin = DI_spread_work(work,current_thread,threads);
    int end   = DI_spread_work(work,current_thread+1,threads);
    --end;
    zend   = (end   % nz) + zmin;
    zstart = (begin % nz) + zmin;
    end   /= nz;
    begin /= nz;
    yend   = (end   % ny) + ymin;
    ystart = (begin % ny) + ymin;
    end   /= ny;
    begin /= ny;
    xend   = end;
    xstart = begin;
  } else {
    zstart = zmin;
    zend   = zmax;
    ystart = ymin;
    yend   = ymax;
    xstart = xs;
    xend   = xe;
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
#else
inline void DataIterator::omp_init(int xs, int xe,bool end){;};
#endif

#endif // __DATAITERATOR_H__
