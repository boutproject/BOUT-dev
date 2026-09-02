#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/unused.hxx>

namespace {
constexpr BoutReal constraint_factor = 0.5;
}

class TestSnesConstraints : public PhysicsModel {
public:
  Field3D u;
  Field3D phi;

  int init(bool UNUSED(restarting)) override {
    solver->add(u, "u");
    solver->constraint(phi, ddt(phi), "phi");

    u = 1.0;
    phi = constraint_factor;

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    ddt(u) = -u + phi;
    ddt(phi) = phi - constraint_factor * u;

    return 0;
  }
};

BOUTMAIN(TestSnesConstraints);
