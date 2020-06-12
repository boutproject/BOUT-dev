/*
 */

#include <bout/physicsmodel.hxx>
#include <bout.hxx>

/// Create a function to be called every timestep
int my_timestep_monitor(Solver* UNUSED(solver), BoutReal simtime, BoutReal dt) {
  output.write("\nTimestep monitor, time = {:e}, dt = {:e}\n", simtime, dt);
  return 0;
}

/// Create a monitor to be called every output step
class MyOutputMonitor : public Monitor {
public:
  explicit MyOutputMonitor(BoutReal timestep = -1) : Monitor(timestep){};
  int call(Solver* UNUSED(solver), BoutReal simtime, int iter, int NOUT) override {
    output.write("\nOutput monitor, time = {:e}, step {:d} of {:d}\n", simtime, iter, NOUT);
    return 0;
  }
};

class MonitorExample : public PhysicsModel {
  Field2D f;

  /// Create instances of the output monitor
  MyOutputMonitor my_output_monitor{};
  /// This output monitor is called twice every output step
  MyOutputMonitor my_output_monitor_fast{.5};

protected:
  int init(bool UNUSED(restarting)) override {
    /// Add the output monitors
    solver->addMonitor(&my_output_monitor);
    solver->addMonitor(&my_output_monitor_fast);

    /// Add the timestep monitor. This also needs to enabled by
    /// setting the input value solver:monitor_timestep to true
    solver->addTimestepMonitor(my_timestep_monitor);

    SOLVE_FOR(f);
    return 0;
  }

  int rhs(BoutReal UNUSED(t)) override {
    ddt(f) = -f;
    return 0;
  }
};

BOUTMAIN(MonitorExample)
