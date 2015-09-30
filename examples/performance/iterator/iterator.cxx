/*
 * Testing performance of an iterator over the mesh
 *
 */

#include <bout.hxx>

#include <algorithm>
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
    xend = mesh->ngx-1;
    yend = mesh->ngy-1;
    zend = mesh->ngz-1;
  }

  MeshIterator(int x, int y, int z) : x(x), y(y), z(z), xstart(0), ystart(0), zstart(0) {
    xend = mesh->ngx-1;
    yend = mesh->ngy-1;
    zend = mesh->ngz-1;
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
  return MeshIterator(mesh->ngx-1, mesh->ngy-1, mesh->ngz-1);
}

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  Field3D a = 1.0;
  Field3D b = 2.0;

  Field3D result;
  result.allocate();

  BoutReal ***ad = a.getData();
  BoutReal ***bd = b.getData();
  BoutReal ***rd = result.getData();

  typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
  typedef std::chrono::duration<double> Duration;
  using namespace std::chrono;
  
  // Nested loops over block data
  SteadyClock start1 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(int i=0;i<mesh->ngx;++i) {
      for(int j=0;j<mesh->ngy;++j) {
        for(int k=0;k<mesh->ngz;++k) {
          rd[i][j][k] = ad[i][j][k] + bd[i][j][k];
        }
      }
    }
  }
  Duration elapsed1 = steady_clock::now() - start1;

  // A single loop over block data
  SteadyClock start2 = steady_clock::now();
  for(int i=0;i<10;++i) {
    for(int j=0;j<mesh->ngx*mesh->ngy*mesh->ngz;++j) {
      rd[0][0][j] = ad[0][0][j] + bd[0][0][j];
    }
  }
  Duration elapsed2 = steady_clock::now() - start2;

  // MeshIterator over block data
  // Has default ctor
  SteadyClock start3 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(MeshIterator i; !i.isDone(); ++i){
      rd[i.x][i.y][i.z] = ad[i.x][i.y][i.z] + bd[i.x][i.y][i.z];
    }
  }
  Duration elapsed3 = steady_clock::now() - start3;

  // DataIterator over block data
  // No default ctor
  SteadyClock start4 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(DataIterator i = begin(result); i != end(result) ; ++i){
      rd[i.x][i.y][i.z] = ad[i.x][i.y][i.z] + bd[i.x][i.y][i.z];
    }
  }
  Duration elapsed4 = steady_clock::now() - start4;

  // DataIterator over block data with done()
  // No default ctor
  SteadyClock start5 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(DataIterator i = begin(result); !i.done() ; ++i){
      rd[i.x][i.y][i.z] = ad[i.x][i.y][i.z] + bd[i.x][i.y][i.z];
    }
  }
  Duration elapsed5 = steady_clock::now() - start5;

  // Range based for DataIterator over data
  SteadyClock start6 = steady_clock::now();
  for(int x=0;x<10;x++) {
    for(auto i : result){
      rd[i.x][i.y][i.z] = ad[i.x][i.y][i.z] + bd[i.x][i.y][i.z];
    }
  }
  Duration elapsed6 = steady_clock::now() - start6;

  // for_each DataIterator over data
  SteadyClock start7 = steady_clock::now();
  for(int x=0;x<10;x++) {
  std::for_each(begin(result), end(result), [&](Indices i) {
      rd[i.x][i.y][i.z] = ad[i.x][i.y][i.z] + bd[i.x][i.y][i.z];
    });
  }
  Duration elapsed7 = steady_clock::now() - start7;

  // Nested loop over fields
  SteadyClock start8 = steady_clock::now();
  for(int x=0;x<10;x++)
    for(int i=0;i<mesh->ngx;++i)
      for(int j=0;j<mesh->ngy;++j)
        for(int k=0;k<mesh->ngz;++k)
          result(i,j,k) = a(i,j,k) + b(i,j,k);
  Duration elapsed8 = steady_clock::now() - start8;

  // Range based DataIterator over fields
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

  // for_each DataIterator over field
  SteadyClock start11 = steady_clock::now();
  for(int x=0;x<10;x++) {
  std::for_each(begin(result), end(result), [&](Indices i) {
      result[i] = a[i] + b[i];
    });
  }
  Duration elapsed11 = steady_clock::now() - start11;

  output << "TIMING\n======\n";
  output << "Nested (data)              : " << elapsed1.count() << std::endl;
  output << "C loop (data)              : " << elapsed2.count() << std::endl;
  output << "MeshIterator (data)        : " << elapsed3.count() << std::endl;
  output << "DataIterator (data)        : " << elapsed4.count() << std::endl;
  output << "DataIterator (data/done()) : " << elapsed5.count() << std::endl;
  output << "Range-based (data)         : " << elapsed6.count() << std::endl;
  output << "For_each (data)            : " << elapsed7.count() << std::endl;
  output << "Nested (field)             : " << elapsed8.count() << std::endl;
  output << "Range-based (field)        : " << elapsed9.count() << std::endl;
  output << "DataIterator (field)       : " << elapsed10.count() << std::endl;
  output << "For_each (field)           : " << elapsed11.count() << std::endl;

  // //----------------------------------------
  // // Operating on a single field

  // // Nested loops
  // SteadyClock start9 = steady_clock::now();
  // for(int x=0;x<10;++x) {
  //   for(int i=0;i<mesh->ngx;++i) {
  //     for(int j=0;j<mesh->ngy;++j) {
  //       for(int k=0;k<mesh->ngz;++k) {
  //         rd[i][j][k] = 3;
  //       }
  //     }
  //   }
  // }
  // Duration elapsed9 = steady_clock::now() - start9;

  // // Single loop
  // SteadyClock start10 = steady_clock::now();
  // for(int x=0;x<10;x++) {
  //   for(int j=0;j<mesh->ngx*mesh->ngy*mesh->ngz;++j) {
  //     rd[0][0][j] = 3;
  //   }
  // }
  // Duration elapsed10 = steady_clock::now() - start10;

  // // MeshIterator
  // SteadyClock start11 = steady_clock::now();
  // for(int i=0;i<10;++i) {
  //   for(MeshIterator i; !i.isDone(); ++i){
  //     rd[i.x][i.y][i.z] = 3;
  //   }
  // }
  // Duration elapsed11 = steady_clock::now() - start11;

  // // Range-based for (FieldIterator)
  // // SteadyClock start12 = steady_clock::now();
  // // for(int x=0;x<10;x++) {
  // //   for (FieldIterator r_it : result) {
  // //     r_it = 3;
  // //   }
  // // }
  // // Duration elapsed12 = steady_clock::now() - start12;

  // // Range-based for (MeshIterator)
  // SteadyClock start13 = steady_clock::now();
  // for(int x=0;x<10;x++) {
  //   for (auto i : mesh) {
  //     rd[i.x][i.y][i.z] = 3;
  //   }
  // }
  // Duration elapsed13 = steady_clock::now() - start13;

  // output << "TIMING (single fields) \n======================\n";
  // output << "Nested (single)			: " << elapsed9.count() << std::endl;
  // output << "C loop (single)			: " << elapsed10.count() << std::endl;
  // output << "Iterator (single)		: " << elapsed11.count() << std::endl;
  // // output << "Range-based (field) (single)	: " << elapsed12.count() << std::endl;
  // output << "Range-based (mesh) (single)	: " << elapsed13.count() << std::endl;


  BoutFinalise();
  return 0;
}
