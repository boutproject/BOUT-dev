/*
 */

#include <bout.hxx>
#include <boutmain.hxx>

Field2D f;

int my_monitor(Solver *solver, BoutReal simtime, int iter, int NOUT) {
  output.write("\nMy monitor, time = %e, dt = %e\n", 
               simtime, solver->getCurrentTimestep());
  return 0;
}

int physics_init(bool restarting) {
  solver->addMonitor(my_monitor);
  SOLVE_FOR(f);
  return 0;
}

int physics_run(BoutReal t) {
  ddt(f) = -f;
  return 0;
}
