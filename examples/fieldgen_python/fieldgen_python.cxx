#include <bout/physicsmodel.hxx>


class FieldGenPython : public PhysicsModel {
public:
  int init(bool restarting) {
    SOLVE_FOR(f);
    return 0;
  }
  int rhs(BoutReal time) {
    ddt(f) = 0.;
    return 0;
  }
  
private:
  Field3D f;
};


BOUTMAIN(FieldGenPython);
