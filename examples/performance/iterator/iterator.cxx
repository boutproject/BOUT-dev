/*
 * Testing performance of an iterator over the mesh
 *
 */

#include <bout.hxx>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <time.h>

// Set of indices - MeshIterator is dereferenced into these
struct Indices {
  int x;
  int y;
  int z;
};

// A simple iterator over a 3D set of indices
class MeshIterator {
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

// Essentially the same as MeshIterator, except it iterates over Field3D directly
class FieldIterator {
public:
  /// Constructor. This would set ranges. Could depend on thread number
  FieldIterator(Field3D& field) :
    x(0), y(0), z(0), xstart(0), ystart(0), zstart(0) {
    xend = mesh->ngx-1;
    yend = mesh->ngy-1;
    zend = mesh->ngz-1;
    data = field.getData();
  }
  FieldIterator(Field3D& field, int x, int y, int z) :
    x(x), y(y), z(z), xstart(0), ystart(0), zstart(0) {
    xend = mesh->ngx-1;
    yend = mesh->ngy-1;
    zend = mesh->ngz-1;
    data = field.getData();
  }

  /// The index variables, updated during loop
  int x, y, z;

  /// Increment operators
  FieldIterator& operator++() { next(); return *this; }
  FieldIterator& operator++(int) { next(); return *this; }

  // Comparison operator
  bool operator!=(const FieldIterator& rhs) const {
    return (x != rhs.x) || (y != rhs.y) || (z != rhs.z);
  }

  // Dereference operator
  BoutReal& operator*() {
    return data[x][y][z];
  }

  /// Checks if finished looping. Is this more efficient than
  /// using the more idiomatic it != MeshIterator::end() ?
  bool isDone() const {
    return x > xend;
  }

  BoutReal*** data;

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

FieldIterator begin(Field3D& field) {
  return FieldIterator(field, 0, 0, 0);
}

FieldIterator end(Field3D& field) {
  return FieldIterator(field, mesh->ngx-1, mesh->ngy-1, mesh->ngz-1);
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
  
  // Nested loops
  clock_t time1 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start1, end1;
  start1 = std::chrono::steady_clock::now();
  for(int x=0;x<10;x++) {
    for(int i=0;i<mesh->ngx;i++) {
      for(int j=0;j<mesh->ngy;j++) {
        for(int k=0;k<mesh->ngz;k++) {
          rd[i][j][k] = ad[i][j][k] + bd[i][j][k];
        }
      }
    }
  }
  end1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed1 = end1 - start1;
  time1 = clock() - time1;

  // A single loop
  clock_t time2 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start2, end2;
  start2 = std::chrono::steady_clock::now();
  for(int i=0;i<10;i++) {
    for(int j=0;j<mesh->ngx*mesh->ngy*mesh->ngz;j++) {
      rd[0][0][j] = ad[0][0][j] + bd[0][0][j];
    }
  }
  end2 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed2 = end2 - start2;
  time2 = clock() - time2;

  // Mesh Iterator
  clock_t time3 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start3, end3;
  start3 = std::chrono::steady_clock::now();
  for(int x=0;x<10;x++) {
    for(MeshIterator i; !i.isDone(); i++){
      rd[i.x][i.y][i.z] = ad[i.x][i.y][i.z] + bd[i.x][i.y][i.z];
    }
  }
  end3 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed3 = end3 - start3;
  time3 = clock() - time3;


  // Triple threat Field3D iterators
  clock_t time4 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start4, end4;
  start4 = std::chrono::steady_clock::now();
  for (int x=0;x<10;++x) {
    for (auto a_it(begin(a)), b_it(begin(b)), r_it(begin(result));
	 a_it != end(a);
	 ++a_it, ++b_it, ++r_it) {
      *r_it = *a_it + *b_it;
    }
  }
  end4 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed4 = end4 - start4;
  time4 = clock() - time4;

  // Range based for MeshIterator
  clock_t time5 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start5, end5;
  start5 = std::chrono::steady_clock::now();
  for(int x=0;x<10;x++) {
    for(auto i : mesh){
      rd[i.x][i.y][i.z] = ad[i.x][i.y][i.z] + bd[i.x][i.y][i.z];
    }
  }
  end5 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed5 = end5 - start5;
  time5 = clock() - time5;


  // for_each MeshIterator
  clock_t time6 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start6, end6;
  start6 = std::chrono::steady_clock::now();
  for(int x=0;x<6;x++) {
  std::for_each(begin(mesh), end(mesh), [&](Indices i) {
      rd[i.x][i.y][i.z] = ad[i.x][i.y][i.z] + bd[i.x][i.y][i.z];
    });
  }
  end6 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed6 = end6 - start6;
  time6 = clock() - time6;

  //////////////

  clock_t time7 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start7, end7;
  start7 = std::chrono::steady_clock::now();
  for(int x=0;x<10;x++)
    for(int i=0;i<mesh->ngx;i++)
      for(int j=0;j<mesh->ngy;j++) 
        for(int k=0;k<mesh->ngz;k++)
          result(i,j,k) = a(i,j,k) + b(i,j,k);
  end7 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed7 = end7 - start7;
  time7 = clock() - time7;

  clock_t time8 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start8, end8;
  start8 = std::chrono::steady_clock::now();
  for(int x=0;x<10;x++)
    for(DataIterator d = result.iterator(); !d.done(); d++)
      result[d] = a[d] + b[d];
  end8 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed8 = end8 - start8;
  time8 = clock() - time8;

  output << "TIMING\n======\n";
  output << "Nested		: " << time1 << "\t:" << elapsed1.count() << std::endl;
  output << "C loop		: " << time2 << "\t:" << elapsed2.count() << std::endl;
  output << "Iterator		: " << time3 << "\t:" << elapsed3.count() << std::endl;
  output << "Field Iterator	: " << time4 << "\t:" << elapsed4.count() << std::endl;
  output << "Range-based	: " << time5 << "\t:" << elapsed5.count() << std::endl;
  output << "For_each		: " << time6 << "\t:" << elapsed6.count() << std::endl;
  output << "Field loops	: " << time7 << "\t:" << elapsed7.count() << std::endl;
  output << "Field iterator	: " << time8 << "\t:" << elapsed8.count() << std::endl;

  //----------------------------------------
  // Operating on a single field

  // Nested loops
  clock_t time9 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start9, end9;
  start9 = std::chrono::steady_clock::now();
  for(int x=0;x<10;++x) {
    for(int i=0;i<mesh->ngx;i++) {
      for(int j=0;j<mesh->ngy;j++) {
        for(int k=0;k<mesh->ngz;k++) {
          rd[i][j][k] = 3;
        }
      }
    }
  }
  end9 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed9 = end9 - start9;
  time9 = clock() - time9;

  // Single loop
  clock_t time10 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start10, end10;
  start10 = std::chrono::steady_clock::now();
  for(int x=0;x<10;x++) {
    for(int j=0;j<mesh->ngx*mesh->ngy*mesh->ngz;j++) {
      rd[0][0][j] = 3;
    }
  }
  end10 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed10 = end10 - start10;
  time10 = clock() - time10;

  // MeshIterator
  clock_t time11 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start11, end11;
  start11 = std::chrono::steady_clock::now();
  for(int i=0;i<10;i++) {
    for(MeshIterator i; !i.isDone(); i++){
      rd[i.x][i.y][i.z] = 3;
    }
  }
  end11 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed11 = end11 - start11;
  time11 = clock() - time11;

  // Range-based for (FieldIterator)
  // clock_t time12 = clock();
  // std::chrono::time_point<std::chrono::steady_clock> start12, end12;
  // start12 = std::chrono::steady_clock::now();
  // for(int x=0;x<10;x++) {
  //   for (FieldIterator r_it : result) {
  //     r_it = 3;
  //   }
  // }
  // end12 = std::chrono::steady_clock::now();
  // std::chrono::duration<double> elapsed12 = end12 - start12;
  // time12 = clock() - time12;

  // Range-based for (MeshIterator)
  clock_t time13 = clock();
  std::chrono::time_point<std::chrono::steady_clock> start13, end13;
  start13 = std::chrono::steady_clock::now();
  for(int x=0;x<10;x++) {
    for (auto i : mesh) {
      rd[i.x][i.y][i.z] = 3;
    }
  }
  end13 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed13 = end13 - start13;
  time13 = clock() - time13;

  output << "TIMING (single fields) \n======================\n";
  output << "Nested (single)			: " << time9 << "\t:" << elapsed9.count() << std::endl;
  output << "C loop (single)			: " << time10 << "\t:" << elapsed10.count() << std::endl;
  output << "Iterator (single)		: " << time11 << "\t:" << elapsed11.count() << std::endl;
  // output << "Range-based (field) (single)	: " << time12 << "\t:" << elapsed12.count() << std::endl;
  output << "Range-based (mesh) (single)	: " << time13 << "\t:" << elapsed13.count() << std::endl;


  BoutFinalise();
  return 0;
}
