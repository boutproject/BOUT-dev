/*
 */

#include <bout/bout_types.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/unused.hxx>

class SquashRun : public PhysicsModel {
protected:
  // Initialisation
  int init(bool UNUSED(restarting)) override {
    solver->add(f2, "f2");
    solver->add(f3, "f3");
    return 0;
  }

  // Calculate time-derivatives
  int rhs(BoutReal UNUSED(t)) override {
    ddt(f2) = 1;
    ddt(f3) = -1;
    f2.applyBoundary();
    f3.applyBoundary();
    return 0;
  }

private:
  Field2D f2;
  Field3D f3;
};

// Create a default main()
BOUTMAIN(SquashRun);
