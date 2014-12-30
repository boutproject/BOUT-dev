/*
 * Timing of arithmetic operations
 *
 */

#include <bout/physicsmodel.hxx>

#include <bout/expr.hxx>

#include <time.h>

class Arithmetic : public PhysicsModel {
protected: 
  int init(bool restarting) {
    
    Field3D a = 1.0;
    Field3D b = 2.0;
    Field3D c = 3.0;

    Field3D result1, result2, result3, result4;

    // Using Field methods (classic operator overloading)
    
    result1 = 2.*a + b * c;

    clock_t c1 = clock();
    result1 = 2.*a + b * c;
    clock_t c2 = clock();
    
    // Using C loops
    result2.allocate();
    BoutReal *rd = result2[0][0];
    BoutReal *ad = a[0][0];
    BoutReal *bd = b[0][0];
    BoutReal *cd = c[0][0];
    for(int i=0, iend=(mesh->ngx*mesh->ngy*mesh->ngz)-1; i != iend; i++) {
      *rd = 2.*(*ad) + (*bd)*(*cd);
      rd++;
      ad++;
      bd++;
      cd++;
    }
    clock_t c3 = clock();
    
    // Template expressions
    
    result3 = eval3D(add(mul(2,a), mul(b,c)));

    clock_t c4 = clock();

    // Range iterator
    result4.allocate();
    for(auto i : result4)
      result4[i] = 2.*a[i] + b[i] * c[i];
    
    clock_t c5 = clock();

    output << "TIMING\n======\n";
    output << "Fields: " << c2 - c1 << endl;
    output << "C loop: " << c3 - c2 << endl;
    output << "Templates: " << c4 - c3 << endl;
    output << "Range For: " << c5 - c4 << endl;
    output << "\nUnits: " << 1./((double)CLOCKS_PER_SEC) << " seconds" << endl << endl;
    
    return 1;
  }
};

BOUTMAIN(Arithmetic);
