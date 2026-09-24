#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/unused.hxx>

class TestSolver : public PhysicsModel {
public:
  Field3D f, g;

  int init(bool UNUSED(restarting)) override {
    solver->add(f, "f");
    solver->add(g, "g");

    f = 1.0;
    g = 0.0;

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    ddt(f) = -0.1 * f;
    ddt(g) = 0.5 * f - g;
    return 0;
  }
};

BOUTMAIN(TestSolver);
