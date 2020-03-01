#include <bout/physicsmodel.hxx>
#include "unused.hxx"

class RestartTest : public PhysicsModel {
private:
  Field3D f3d;
  Field2D f2d;
protected:
  int init(bool UNUSED(restarting)) override {
    // Evolve a 3D and a 2D field
    SOLVE_FOR2(f3d, f2d);
    return 0;
  }
  int rhs(BoutReal UNUSED(time)) override {
    // Simple time evolution 
    ddt(f3d) = 0.1*f3d;
    ddt(f2d) = -0.1*f2d;
    return 0;
  }
};

BOUTMAIN(RestartTest);
