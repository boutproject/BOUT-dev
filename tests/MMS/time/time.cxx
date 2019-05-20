/*
 * MMS test of time integration methods
 * 
 * This code does very little, as all sources are added
 * in the solver code itself.
 */


#include <bout/physicsmodel.hxx>
#include "unused.hxx"

class TimeTest : public PhysicsModel {
public:
  int init(MAYBE_UNUSED(bool restart)) {
    solver->add(f, "f"); // Solve a single 3D field
    
    setSplitOperator();
    
    return 0;
  }
  
  int rhs(MAYBE_UNUSED(BoutReal time)) {
    ddt(f) = f;
    return 0;
  }
  
  int convective(MAYBE_UNUSED(BoutReal time)) {
    ddt(f) = 0.5*f;
    return 0;
  }
  int diffusive(MAYBE_UNUSED(BoutReal time)) {
    ddt(f) = 0.5*f;
    return 0;
  }
private:
  Field3D f;      
};

BOUTMAIN(TimeTest);

