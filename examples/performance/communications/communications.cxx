#include <bout/physicsmodel.hxx>
#include <bout/sys/timer.hxx>

class CommsSpeed : public PhysicsModel {
public:
  Field3D f;

  int init(bool) {

    SOLVE_FOR(f);

    return 0;
  }

  int rhs(BoutReal) {

    mesh->communicate(f);

    ddt(f) = 1.;

    return 0;
  }
};

BOUTMAIN(CommsSpeed);
