/*
 * MMS test of time integration methods
 * 
 * This code does very little, as all sources are added
 * in the solver code itself.
 */

#include "bout/unused.hxx"
#include <bout/physicsmodel.hxx>

class TimeTest : public PhysicsModel {
public:
  int init([[maybe_unused]] bool restart) override {
    solver->add(f, "f"); // Solve a single 3D field

    setSplitOperator();

    return 0;
  }

  int rhs([[maybe_unused]] BoutReal time) override {
    ddt(f) = f;
    return 0;
  }

  int convective([[maybe_unused]] BoutReal time) {
    ddt(f) = 0.5 * f;
    return 0;
  }
  int diffusive([[maybe_unused]] BoutReal time) {
    ddt(f) = 0.5 * f;
    return 0;
  }

private:
  Field3D f;
};

BOUTMAIN(TimeTest);
