#include <bout/physicsmodel.hxx>
#include "unused.hxx"

class TestDelp2 : public PhysicsModel {
protected:

  int init(bool UNUSED(restarting)) {
    Options *opt = Options::getRoot()->getSection("diffusion");
    OPTION(opt, D, 0.1);
  
    SOLVE_FOR(n);
  
    return 0;
  }

  int rhs(BoutReal UNUSED(t)) {
    mesh->communicate(n);
  
    ddt(n) = D * Delp2(n);
  
    return 0;
  }

private:
  Field3D n;
  BoutReal D;
};

BOUTMAIN(TestDelp2);
