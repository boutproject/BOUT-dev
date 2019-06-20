#include <bout/physicsmodel.hxx>
#include "unused.hxx"

class TestDelp2 : public PhysicsModel {
protected:

  int init(bool UNUSED(restarting)) override {
    Options *opt = Options::getRoot()->getSection("diffusion");
    OPTION(opt, D, 0.1);
    OPTION(opt, useFFT, true);

    SOLVE_FOR(n);
  
    return 0;
  }

  int rhs(BoutReal UNUSED(t)) override {
    mesh->communicate(n);

    ddt(n) = D * Delp2(n, CELL_DEFAULT, useFFT);

    return 0;
  }

private:
  Field3D n;
  BoutReal D;
  bool useFFT;
};

BOUTMAIN(TestDelp2);
