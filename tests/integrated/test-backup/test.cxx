#include <bout/physicsmodel.hxx>
#include "unused.hxx"

class Test : public PhysicsModel {
protected:

  int init(bool UNUSED(restarting)) {

    SOLVE_FOR(n);
  
    return 0;
  }

  int rhs(BoutReal UNUSED(t)) {
    ddt(n) = 0;

    return 0;
  }

private:
  Field3D n;
};

BOUTMAIN(Test);
