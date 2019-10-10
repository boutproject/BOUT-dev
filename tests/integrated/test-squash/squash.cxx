/*
 */

#include <bout/physicsmodel.hxx>

class SquashRun : public PhysicsModel {
protected:
  // Initialisation
  int init(bool UNUSED(restarting)) {
    solver->add(f2, "f2");
    solver->add(f3, "f3");
    return 0;
  }
  
  // Calculate time-derivatives
  int rhs(BoutReal UNUSED(t)) {
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
