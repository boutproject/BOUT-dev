/*
 */

#include <bout/physicsmodel.hxx>

class Test : public PhysicsModel {
private:
  Field2D f;
  
public:
  int init(bool) override {
    SOLVE_FOR(f);
    return 0;
  }
  
  int rhs(BoutReal) override {
    ddt(f) = -f;
    return 0;
  }
};

BOUTMAIN(Test);

