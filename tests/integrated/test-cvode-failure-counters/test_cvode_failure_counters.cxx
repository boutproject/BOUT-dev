// Test of the CVODE option solver:cvode_failure_counters
//
// A stiff, forced, nonlinear diffusion equation in y. With a small Krylov
// subspace and no preconditioner, Krylov solves fail. Recoverable RHS failures
// are thrown in the Newton residual at fixed calls after each output, so their
// number is known.

#include <bout/boutexception.hxx>
#include <bout/derivs.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/unused.hxx>

#include <cmath>

class TestCvodeFailureCounters : public PhysicsModel {
public:
  int init(bool UNUSED(restarting)) override {
    auto& options = Options::root()["model"];
    diffusion = options["diffusion"].doc("Diffusion coefficient").withDefault(100.0);
    inject_rhs_failures =
        options["inject_rhs_failures"]
            .doc("Recoverable RHS failures in the Newton residual per output interval")
            .withDefault(2);

    solver->add(f, "f");
    // The initial profile is also the shape of the periodic forcing
    source = f;
    return 0;
  }

  int rhs(BoutReal time, bool linear) override {
    if (not linear) {
      // CVODE uses the nonlinear RHS only in the Newton residual. Fail on every
      // second call after an output, at most inject_rhs_failures times
      ++calls_since_output;
      if (failures_since_output < inject_rhs_failures and calls_since_output % 2 == 0) {
        ++failures_since_output;
        ++rhs_failures;
        throw BoutRhsFail("Injected recoverable RHS failure");
      }
    }
    mesh->communicate(f);
    ddt(f) = diffusion * D2DY2(f) - f * f * f + source * std::cos(time);
    return 0;
  }

  int outputMonitor(BoutReal UNUSED(simtime), int UNUSED(iter),
                    int UNUSED(nout)) override {
    calls_since_output = 0;
    failures_since_output = 0;
    return 0;
  }

  void outputVars(Options& options) override {
    options["rhs_failures"].assignRepeat(rhs_failures);
  }

private:
  Field3D f;
  Field3D source;
  BoutReal diffusion{0.0};
  int inject_rhs_failures{0};
  int calls_since_output{0};
  int failures_since_output{0};
  /// Total number of injected RHS failures
  int rhs_failures{0};
};

BOUTMAIN(TestCvodeFailureCounters);
