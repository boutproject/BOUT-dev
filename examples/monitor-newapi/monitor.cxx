/*
 */

#include <bout/physicsmodel.hxx>

class MonitorExample : public PhysicsModel {
protected:
  // Initialisation
  int init(bool UNUSED(restarting)) {
    solver->add(f, "f");
    return 0;
  }
  
  // Calculate time-derivatives
  int rhs(BoutReal UNUSED(t)) {
    ddt(f) = -f;
    return 0;
  } 
  
  // This called every output timestep
  int outputMonitor(BoutReal simtime, int iter, int NOUT) {
    output.write("\nOutput monitor, time = {:e}, step {:d} of {:d}\n",
                 simtime, iter, NOUT);
    return 0;
  }
  
  // This called every timestep
  int timestepMonitor(BoutReal simtime, BoutReal dt) {
    output.write("\nTimestep monitor, time = {:e}, dt = {:e}\n", simtime, dt);
    return 0;
  }
  
private:
  Field2D f;
};

// Create a default main()
BOUTMAIN(MonitorExample);
