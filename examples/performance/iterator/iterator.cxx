/*
 * Testing performance of an iterator over the mesh
 *
 */

#include <bout.hxx>

#include <chrono>
#include <iostream>
#include <iterator>
#include <time.h>

// A simple iterator over a 3D set of indices
class MeshIterator
  : public std::iterator< std::forward_iterator_tag, Indices > {
public:
  /// Constructor. This would set ranges. Could depend on thread number
  MeshIterator() : x(0), y(0), z(0), xstart(0), ystart(0), zstart(0) {
    xend = mesh->LocalNx-1;
    yend = mesh->LocalNy-1;
    zend = mesh->LocalNz;
  }

  MeshIterator(int x, int y, int z) : x(x), y(y), z(z), xstart(0), ystart(0), zstart(0) {
    xend = mesh->LocalNx-1;
    yend = mesh->LocalNy-1;
    zend = mesh->LocalNz;
  }

  /// The index variables, updated during loop
  int x, y, z;

  /// Increment operators
  MeshIterator& operator++() { next(); return *this; }
  MeshIterator& operator++(int) { next(); return *this; }

  // Comparison operator
  bool operator!=(const MeshIterator& rhs) const {
    return (x != rhs.x) || (y != rhs.y) || (z != rhs.z);
  }

  // Dereference operator
  Indices operator*() {
    return {x, y, z};
  }

  /// Checks if finished looping. Is this more efficient than
  /// using the more idiomatic it != MeshIterator::end() ?
  bool isDone() const {
    return x > xend;
  }

private:
  int xstart, xend;
  int ystart, yend;
  int zstart, zend;

  /// Advance to the next index
  void next() {
    z++;
    if(z > zend) {
      z = zstart;
      y++;
      if(y > yend) {
        y = ystart;
        x++;
      }
    }
  }
};

// Begin/end iterators
MeshIterator begin(Mesh* mesh) {
  return MeshIterator(0, 0, 0);
}

MeshIterator end(Mesh* mesh) {
  return MeshIterator(mesh->LocalNx-1, mesh->LocalNy-1, mesh->LocalNz);
}

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  Field3D a = 1.0;
  Field3D b = 2.0;

  Field3D result;
  result.allocate();

  typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
  typedef std::chrono::duration<double> Duration;
  using namespace std::chrono;
  
  // A single loop over block data
  
  BoutReal *ad = &a(0,0,0);
  BoutReal *bd = &b(0,0,0);
  BoutReal *rd = &result(0,0,0);
  
  // Loop over data so first test doesn't have a disadvantage from caching
  for(int i=0;i<10;++i) {
    for(int j=0;j<mesh->LocalNx*mesh->LocalNy*mesh->LocalNz;++j) {
      rd[j] = ad[j] + bd[j];
    }
  }
  
  SteadyClock start1 = steady_clock::now();
  int len = mesh->LocalNx*mesh->LocalNy*mesh->LocalNz;
  for(int i=0;i<10;++i) {
    for(int j=0;j<len;++j) {
      rd[j] = ad[j] + bd[j];
    }
  }
  Duration elapsed1 = steady_clock::now() - start1;

  // Nested loops over block data
  SteadyClock start2 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(int i=0;i<mesh->LocalNx;++i) {
      for(int j=0;j<mesh->LocalNy;++j) {
        for(int k=0;k<mesh->LocalNz;++k) {
          result(i,j,k) = a(i,j,k) + b(i,j,k);
        }
      }
    }
  }
  Duration elapsed2 = steady_clock::now() - start2;

  // MeshIterator over block data
  SteadyClock start3 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(MeshIterator i; !i.isDone(); ++i){
      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y,i.z);
    }
  }
  Duration elapsed3 = steady_clock::now() - start3;

  // DataIterator using begin(), end()
  SteadyClock start4 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(DataIterator i = begin(result), rend=end(result); i != rend; ++i){
      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y,i.z);
    }
  }
  Duration elapsed4 = steady_clock::now() - start4;

  // DataIterator with done()
  SteadyClock start5 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(DataIterator i = begin(result); !i.done() ; ++i){
      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y,i.z);
    }
  }
  Duration elapsed5 = steady_clock::now() - start5;

  // Range based for DataIterator with indices
  SteadyClock start6 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(auto i : result){
      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y,i.z);
    }
  }
  Duration elapsed6 = steady_clock::now() - start6;
  
  // Range based DataIterator 
  SteadyClock start9 = steady_clock::now();
  for (int x=0;x<10;++x) {
    for (auto i : result) {
      result[i] = a[i] + b[i];
    }
  }
  Duration elapsed9 = steady_clock::now() - start9;

  // DataIterator over fields
  SteadyClock start10 = steady_clock::now();
  for(int x=0;x<10;x++)
    for(DataIterator d = result.iterator(); !d.done(); d++)
      result[d] = a[d] + b[d];
  Duration elapsed10 = steady_clock::now() - start10;
  
  output << "TIMING\n======\n";
  output << "C loop                     : " << elapsed1.count() << std::endl;
  output << "----- (x,y,z) indexing ----" << std::endl;
  output << "Nested loops               : " << elapsed2.count() << std::endl;
  output << "MeshIterator               : " << elapsed3.count() << std::endl;
  output << "DataIterator (begin/end)   : " << elapsed4.count() << std::endl;
  output << "DataIterator (begin/done)  : " << elapsed5.count() << std::endl;
  output << "C++11 range-based for      : " << elapsed6.count() << std::endl;
  output << "------ [i] indexing -------" << std::endl;
  output << "C++11 Range-based for      : " << elapsed9.count() << std::endl;
  output << "DataIterator (done)        : " << elapsed10.count() << std::endl;

  BoutFinalise();
  return 0;
}
