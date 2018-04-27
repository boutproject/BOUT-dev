#include "bout/constants.hxx"
#include "bout/physicsmodel.hxx"
#include "bout/solverfactory.hxx"

#include <cmath>
#include <string>
#include <vector>

class TestSolver : public PhysicsModel {
public:
  Field3D field;

  int init(bool restarting) {
    solver->add(field, "field");
    return 0;
  }

  int rhs(BoutReal time) {
    ddt(field) = sin(time) * sin(time);
    return 0;
  }
};

int main(int argc, char **argv) {

  Output output_test;

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

  // TODO: turn off std::cout
  std::stringstream buffer;
  // Save cout's buffer here
  std::streambuf *sbuf = std::cout.rdbuf();
  // Redirect cout to our stringstream buffer or any other ostream
  std::cout.rdbuf(buffer.rdbuf());

  int init_err = BoutInitialise(argc, argv);
  if (init_err < 0) {
    return 0;
  } else if (init_err > 0) {
    return init_err;
  }

  // When done redirect cout to its old self
  std::cout.rdbuf(sbuf);
  Output::getInstance()->disable();

  for (auto &name : SolverFactory::getInstance()->listAvailable()) {

    output_test << "\n**************************************************\n"
                << "Now running with: " << name << "\n\n";
    try {
      // Get specific options section for this solver. Can't just use default
      // "solver" section, as we run into problems when solvers use the same
      // name for an option with inconsistent defaults
      auto options = Options::getRoot()->getSection(name);
      Solver *solver = SolverFactory::getInstance()->createSolver(name, options);

      TestSolver *model = new TestSolver();
      solver->setModel(model);

      Monitor *bout_monitor = new BoutMonitor();
      solver->addMonitor(bout_monitor, Solver::BACK);

      std::string field_name = "field_" + name;

      solver->solve();

      BoutReal expected = PI / 4.;
      BoutReal tolerance = 1.e-5;
      output_test << "Solver " << name;
      if (fabs(model->field(1, 1, 0) - expected) > tolerance) {
        output_test << ": FAILED\n"
                    << "    Got: "<< model->field(1, 1, 0) << ", expected: " << expected << "\n";
      } else {
        output_test << ": PASSED\n";
      }

      dump.addOnce(model->field, field_name);
      dump.write("data/%s.nc", name.c_str());

      delete model;
      delete solver;
      delete bout_monitor;

    } catch (BoutException &e) {
      // Don't let one bad solver stop us trying the rest
      output_test << "Error encountered with solver " << name << "\n";
      output_test << e.what() << endl;
    }
  }

  BoutFinalise();
}
