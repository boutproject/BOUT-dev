#include "bout/constants.hxx"
#include "bout/physicsmodel.hxx"
#include "bout/solverfactory.hxx"

#include <cmath>
#include <string>
#include <vector>

// A simple phyics model for integrating sin^2(t)
class TestSolver : public PhysicsModel {
public:
  Field3D field;

  int init(bool UNUSED(restarting)) override {
    solver->add(field, "field");
    return 0;
  }

  int rhs(BoutReal time) override {
    ddt(field) = sin(time) * sin(time);
    return 0;
  }
};

int main(int argc, char **argv) {

  // The expected answer to the integral of \f$\int_0^{\pi/2}\sin^2(t)\f$
  BoutReal expected = PI / 4.;
  // Absolute tolerance for difference between the actual value and the
  // expected value
  BoutReal tolerance = 1.e-5;

  // Our own output to stdout, as main library will only be writing to log files
  Output output_test;

  // Currently hardcode solvers we don't want to test
  // Should be able to check which solvers aren't suitable
  std::vector<std::string> eigen_solvers = {"power", "slepc", "snes"};

  for (auto &eigen_solver : eigen_solvers) {
    if (SolverFactory::getInstance()->remove(eigen_solver)) {
      output_test << "Removed '" << eigen_solver << "' eigen solver\n";
    }
  }

  output_test << "\nTesting the following solvers:\n";
  for (auto &solver : SolverFactory::getInstance()->listAvailable()) {
    output_test << "  " << solver << "\n";
  }
  // Explicit flush to make sure list of available solvers gets printed
  output_test << std::endl;

  // DANGER, hack below! BoutInitialise turns on writing to stdout for rank 0,
  // and then immediately prints a load of stuff to stdout. We want this test to
  // be quiet, so we need to hide stdout from the main library before
  // BoutInitialise. After the call to BoutInitialise, we can turn off the main
  // library writing to stdout in a nicer way.

  // Save cout's buffer here
  std::stringstream buffer;
  auto *sbuf = std::cout.rdbuf();
  // Redirect cout to our buffer
  std::cout.rdbuf(buffer.rdbuf());

  int init_err = BoutInitialise(argc, argv);
  if (init_err < 0) {
    return 0;
  } else if (init_err > 0) {
    return init_err;
  }

  // Now BoutInitialise is done, redirect stdout to its old self
  std::cout.rdbuf(sbuf);
  // Turn off writing to stdout for the main library
  Output::getInstance()->disable();

  // Solver and its actual value if it didn't pass
  std::map<std::string, BoutReal> errors;

  for (auto &name : SolverFactory::getInstance()->listAvailable()) {

    output_test << "Testing " << name << " solver:";
    try {
      // Get specific options section for this solver. Can't just use default
      // "solver" section, as we run into problems when solvers use the same
      // name for an option with inconsistent defaults
      auto options = Options::getRoot()->getSection(name);
      auto solver = std::unique_ptr<Solver>{SolverFactory::getInstance()->createSolver(name, options)};

      auto model = bout::utils::make_unique<TestSolver>();
      solver->setModel(model.get());

      auto bout_monitor = bout::utils::make_unique<BoutMonitor>();
      solver->addMonitor(bout_monitor.get(), Solver::BACK);

      solver->solve();

      if (fabs(model->field(1, 1, 0) - expected) > tolerance) {
        output_test << " FAILED\n";
        errors[name] = model->field(1, 1, 0);
      } else {
        output_test << " PASSED\n";
      }
    } catch (BoutException &e) {
      // Don't let one bad solver stop us trying the rest
      output_test << " ERROR\n";
      output_info << "Error encountered with solver " << name << "\n";
      output_info << e.what() << endl;
      errors[name] = 0.;
    }
  }

  BoutFinalise();

  if (!errors.empty()) {
    output_test << "\n => Some failed tests\n\n";
    for (auto &error : errors) {
      output_test << "    " << error.first << " got: " << error.second
                  << ", expected: " << expected << "\n";
    }
  } else {
    output_test << "\n => All tests passed\n";
  }

  return errors.size();
}
