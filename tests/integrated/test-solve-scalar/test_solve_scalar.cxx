#include "bout/constants.hxx"
#include "bout/globals.hxx"
#include "bout/physicsmodel.hxx"
#include "bout/scalar.hxx"

#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>

// A simple phyics model for integrating sin^2(t)
class TestSolver : public PhysicsModel {
public:
  bout::Scalar field{0.0};

  int init(bool UNUSED(restarting)) override {
    solver->add(field, "field");
    return 0;
  }

  int rhs(BoutReal time) override {
    ddt(field) = sin(time) * sin(time);
    return 0;
  }
};

int main() {

  // The expected answer to the integral of \f$\int_0^{\pi/2}\sin^2(t)\f$
  constexpr BoutReal expected = PI / 4.;
  // Absolute tolerance for difference between the actual value and the
  // expected value
  constexpr BoutReal tolerance = 1.e-5;

  // Our own output to stdout, as main library will only be writing to log files
  Output output_test;

  auto& root = Options::root();

  root["mesh"]["MXG"] = 1;
  root["mesh"]["MYG"] = 1;
  root["mesh"]["nx"] = 3;
  root["mesh"]["ny"] = 1;
  root["mesh"]["nz"] = 1;

  root["output"]["enabled"] = false;
  root["restart_files"]["enabled"] = false;
  root["datadir"] = "data";

  // Turn off writing to stdout for the main library
  Output::getInstance()->disable();

  auto mpi_wrapper = std::make_unique<MpiWrapper>();
  bout::globals::mpi = &*mpi_wrapper;

  bout::globals::mesh = Mesh::create();
  bout::globals::mesh->load();

  constexpr BoutReal end = PI / 2.;
  constexpr int NOUT = 100;

  // Global options
  root["nout"] = NOUT;
  root["timestep"] = end / NOUT;

  // Don't error just because we haven't used all the options yet
  root["input"]["error_on_unused_options"] = false;

  BoutReal errors = 0.0;

  try {
    auto solver = Solver::create();

    TestSolver model{};
    solver->setModel(&model);

    BoutMonitor bout_monitor{};
    solver->addMonitor(&bout_monitor, Solver::BACK);

    solver->solve();

    if (std::abs(model.field - expected) > tolerance) {
      errors = model.field;
      output_test.write(" FAILED\n    got: {}, expected: {}\n", errors, expected);
    } else {
      output_test.write(" PASSED\n");
    }
  } catch (BoutException& e) {
    output_test.write(" ERROR\n{}\n", e.what());
    errors = 0.;
  }

  return static_cast<int>(std::abs(errors) > 0.0);
}
