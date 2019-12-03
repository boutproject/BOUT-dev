/*
 */

#include <bout.hxx>
#include <boutmain.hxx>

Field2D f;

// Create a monitor to be called every output step
class MyOutputMonitor: public Monitor{
public:
  MyOutputMonitor(BoutReal timestep=-1):Monitor(timestep){};
  int call(Solver *solver, BoutReal simtime, int iter, int NOUT) override;
};

int MyOutputMonitor::call(Solver *UNUSED(solver), BoutReal simtime, int iter, int NOUT) {
  output.write("Output monitor, time = %e, step %d of %d\n",
               simtime, iter, NOUT);
  return 0;
}

// Create a function to be called every timestep
int my_timestep_monitor(Solver *UNUSED(solver), BoutReal simtime, BoutReal dt) {
  output.write("\nTimestep monitor, time = %e, dt = %e\n", 
               simtime, dt);
  return 0;
}

MyOutputMonitor my_output_monitor;
MyOutputMonitor my_output_monitor_fast(.5);

int physics_init(bool UNUSED(restarting)) {
  solver->addMonitor(&my_output_monitor);
  solver->addMonitor(&my_output_monitor_fast);
  solver->addTimestepMonitor(my_timestep_monitor);
  SOLVE_FOR(f);
  return 0;
}

int physics_run(BoutReal UNUSED(t)) {
  ddt(f) = -f;
  return 0;
}
