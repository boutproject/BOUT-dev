/*
 */

#include <bout.hxx>
#include <boutmain.hxx>

Field2D f;

// Create a function to be called every output step
int my_output_monitor(Solver *solver, BoutReal simtime, int iter, int NOUT) {
  output.write("\nOutput monitor, time = %e, step %d of %d\n", 
               simtime, iter, NOUT);
  return 0;
}

// Create a function to be called every timestep
int my_timestep_monitor(Solver *solver, BoutReal simtime, BoutReal dt) {
  output.write("\nTimestep monitor, time = %e, dt = %e\n", 
               simtime, dt);
  return 0;
}

int physics_init(bool restarting) {
  solver->addMonitor(my_output_monitor);
  solver->addTimestepMonitor(my_timestep_monitor);
  SOLVE_FOR(f);
  return 0;
}

int physics_run(BoutReal t) {
  ddt(f) = -f;
  return 0;
}
