/*


 */

#ifndef __DATAITERATOR_H__
#define __DATAITERATOR_H__

#include <iterator>
#include <iostream>

// Set of indices - DataIterator is dereferenced into these
struct Indices {
  int x;
  int y;
  int z;
};

class DataIterator
  : public std::iterator<std::forward_iterator_tag, Indices> {
public:
  /// Constructor. This would set ranges. Could depend on thread number
  DataIterator(int xs, int xe,
	       int ys, int ye,
	       int zs, int ze) : x(xs), y(ys), z(zs),
				 xstart(xs), xend(xe),
				 ystart(ys), yend(ye),
				 zstart(zs), zend(ze) {
  }
  DataIterator(int x, int xs, int xe,
	       int y, int ys, int ye,
	       int z, int zs, int ze) : x(x), y(y), z(z),
					xstart(xs), xend(xe),
					ystart(ys), yend(ye),
					zstart(zs), zend(ze) {
  }
  
  /// The index variables, updated during loop
  int x, y, z;

  /// Increment operators
  DataIterator& operator++() { next(); return *this; }
  DataIterator operator++(int) { DataIterator tmp(*this); next(); return tmp; }

  // Comparison operator
  inline bool operator!=(const DataIterator& rhs) const {
    return (x != rhs.x) || (y != rhs.y) || (z != rhs.z);
  }

  // Dereference operators
  Indices operator*() {
    return {x, y, z};
  }
  const Indices operator*() const {
    return {x, y, z};
  }

  void start() {
    x = xstart; y = ystart; z = zstart;
  }
  
  /// Checks if finished looping. Is this more efficient than
  /// using the more idiomatic it != MeshIterator::end() ?
  bool done() const {
    return x > xend;
  }
  
private:
  DataIterator(); // Disable null constructor
  
  int xstart, xend;
  int ystart, yend;
  int zstart, zend;
  
  /// Advance to the next index
  void next() {
    ++z;
    if(z > zend) {
      z = zstart;
      ++y;
      if(y > yend) {
	y = ystart;
	++x;
      }
    }
  }
};

/*
template<>
struct DataIterable {
  DataIterable();
}
*/

#endif // __DATAITERATOR_H__
