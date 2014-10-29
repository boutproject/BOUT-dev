/*
 * Testing performance of an iterator over the mesh
 * 
 */

#include <bout.hxx>

#include <time.h>

// A simple iterator over a 3D set of indices
class MeshIterator {
public:
  /// Constructor. This would set ranges. Could depend on thread number
  MeshIterator() : x(0), y(0), z(0), xstart(0), ystart(0), zstart(0) {
    xend = mesh->ngx-1;
    yend = mesh->ngy-1;
    zend = mesh->ngz-1;
  }
  
  /// The index variables, updated during loop
  int x, y, z;

  /// Increment operators
  MeshIterator& operator++() { next(); return *this; }
  MeshIterator& operator++(int) { next(); return *this; }

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
  for(int x=0;x<10;x++) {
    for(int i=0;i<mesh->ngx;i++) {
      for(int j=0;j<mesh->ngy;j++) {
        for(int k=0;k<mesh->ngz;k++) {
          rd[i][j][k] = ad[i][j][k] + bd[i][j][k];
        }
      }
    }
  }
  time1 = clock() - time1;
  
  clock_t time2 = clock();
  // A single loop
  for(int i=0;i<10;i++) {
    for(int j=0;j<mesh->ngx*mesh->ngy*mesh->ngz;j++) {
      rd[0][0][j] = ad[0][0][j] + bd[0][0][j];
    }
  }
  time2 = clock() - time2;
  
  // Iterator
  clock_t time3 = clock();
  for(int x=0;x<10;x++) {
    for(MeshIterator i; !i.isDone(); i++)
      rd[i.x][i.y][i.z] = ad[i.x][i.y][i.z] + bd[i.x][i.y][i.z];
  }
  time3 = clock() - time3;
  
  output << "TIMING\n======\n";
  output << "Nested: " << time1 << endl;
  output << "C loop: " << time2 << endl;
  output << "Iterator: " << time3 << endl;
  
  BoutFinalise();
  return 0;
}
