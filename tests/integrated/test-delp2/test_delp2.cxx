#include <bout/physicsmodel.hxx>
#include <bout/scorepwrapper.hxx>

class TestDelp2 : public PhysicsModel {
protected:

  int init(bool restarting) {
    Options *opt = Options::getRoot()->getSection("diffusion");
    OPTION(opt, D, 0.1);
  
    SOLVE_FOR(n);
  
    return 0;
  }

  int rhs(BoutReal t) {
    SCOREP0();
    mesh->communicate(n);
  
    ddt(n) = D * Delp2(n);
  
    return 0;
  }

private:
  Field3D n;
  BoutReal D;
};

BOUTMAIN(TestDelp2);
