/*
 * Timing of arithmetic operations
 *
 */

#include <bout/physicsmodel.hxx>

#include <bout/expr.hxx>

#include <chrono>

typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
typedef std::chrono::duration<double> Duration;
using namespace std::chrono;

class Arithmetic : public PhysicsModel {
protected: 
  int init(bool restarting) {
    
    Field3D a = 1.0;
    Field3D b = 2.0;
    Field3D c = 3.0;

    Field3D result1, result2, result3, result4;

    // Using Field methods (classic operator overloading)
    
    result1 = 2.*a + b * c;

    SteadyClock start1 = steady_clock::now();
    result1 = 2.*a + b * c;
    Duration elapsed1 = steady_clock::now() - start1;
    
    // Using C loops
    result2.allocate();
    BoutReal *rd = &result2(0,0,0);
    BoutReal *ad = &a(0,0,0);
    BoutReal *bd = &b(0,0,0);
    BoutReal *cd = &c(0,0,0);
    SteadyClock start2 = steady_clock::now();
    for(int i=0, iend=(mesh->LocalNx*mesh->LocalNy*mesh->LocalNz)-1; i != iend; i++) {
      *rd = 2.*(*ad) + (*bd)*(*cd);
      rd++;
      ad++;
      bd++;
      cd++;
    }
    Duration elapsed2 = steady_clock::now() - start2;
    
    // Template expressions
    SteadyClock start3 = steady_clock::now();
    result3 = eval3D(add(mul(2,a), mul(b,c)));
    Duration elapsed3 = steady_clock::now() - start3;
    
    // Range iterator
    result4.allocate();
    SteadyClock start4 = steady_clock::now();
    for(auto i : result4)
      result4[i] = 2.*a[i] + b[i] * c[i];
    Duration elapsed4 = steady_clock::now() - start4;
    
    output << "TIMING\n======\n";
    output << "Fields: " << elapsed1.count() << endl;
    output << "C loop: " << elapsed2.count() << endl;
    output << "Templates: " << elapsed3.count() << endl;
    output << "Range For: " << elapsed4.count() << endl;
    
    return 1;
  }
};

BOUTMAIN(Arithmetic);
