/*
 */

#include <bout/physicsmodel.hxx>

class Test : public PhysicsModel {
private:
  Field2D f;
  
public:
  int init(bool restarting) override {
    SOLVE_FOR(f);
    return 0;
  }
  
  int rhs(BoutReal t) override {
    ddt(f) = -f;
    return 0;
  }
};

BOUTMAIN(Test);

