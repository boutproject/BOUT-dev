#include "bout/constants.hxx"
#include "bout/physicsmodel.hxx"
#include "bout/petsclib.hxx"
#include "bout/slepclib.hxx"

#include <cmath>
#include <memory>
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

int main(int argc, char** argv) {

  // The expected answer to the integral of \f$\int_0^{\pi/2}\sin^2(t)\f$
  constexpr BoutReal expected = PI / 4.;
  // Absolute tolerance for difference between the actual value and the
  // expected value
  constexpr BoutReal tolerance = 1.e-5;

  // Our own output to stdout, as main library will only be writing to log files
  Output output_test;

  // Currently hardcode solvers we don't want to test
  // Should be able to check which solvers aren't suitable
  std::vector<std::string> eigen_solvers = {"power", "slepc", "snes"};

  for (auto& eigen_solver : eigen_solvers) {
    if (SolverFactory::getInstance().remove(eigen_solver)) {
      output_test << "Removed '" << eigen_solver << "' eigen solver\n";
    }
  }

  output_test << "\nTesting the following solvers:\n";
  for (auto& solver : SolverFactory::getInstance().listAvailable()) {
    output_test << "  " << solver << "\n";
  }
  // Explicit flush to make sure list of available solvers gets printed
  output_test << std::endl;

  auto& root = Options::root();

  root["mesh"]["MXG"] = 1;
  root["mesh"]["MYG"] = 1;
  root["mesh"]["nx"] = 3;
  root["mesh"]["ny"] = 1;
  root["mesh"]["nz"] = 1;

  root["output"]["enabled"] = false;
  root["restart"]["enabled"] = false;
  root["datadir"] = "data";
  root["dump_format"] = "nc";

  // Set the command-line arguments
  SlepcLib::setArgs(argc, argv);
  PetscLib::setArgs(argc, argv);
  Solver::setArgs(argc, argv);
  BoutComm::setArgs(argc, argv);

  // Turn off writing to stdout for the main library
  Output::getInstance()->disable();

  bout::globals::mpi = new MpiWrapper();

  bout::globals::mesh = Mesh::create();
  bout::globals::mesh->load();

  bout::globals::dump =
      bout::experimental::setupDumpFile(Options::root(), *bout::globals::mesh, ".");

  constexpr BoutReal end = PI / 2.;
  constexpr int NOUT = 100;

  // Global options
  root["NOUT"] = NOUT;
  root["TIMESTEP"] = end / NOUT;

  // Solver-specific options
  root["euler"]["mxstep"] = 100000;
  root["euler"]["nout"] = NOUT;
  root["euler"]["timestep"] = end / (NOUT * 1000);

  root["rk4"]["adaptive"] = true;

  root["rkgeneric"]["adaptive"] = true;

  root["imexbdf2"]["adaptive"] = true;
  root["imexbdf2"]["adaptRtol"] = 1.e-5;

  root["karniadakis"]["nout"] = 100;
  root["karniadakis"]["timestep"] = end / (NOUT * 10000);

  root["petsc"]["nout"] = 10000;
  root["petsc"]["output_step"] = end / 10000;

  root["snes"]["adaptive"] = true;

  root["splitrk"]["timestep"] = end / (NOUT * 500);
  root["splitrk"]["nstages"] = 3;
  root["splitrk"]["mxstep"] = 10000;
  root["splitrk"]["adaptive"] = false;

  // Solver and its actual value if it didn't pass
  std::map<std::string, BoutReal> errors;

  for (auto& name : SolverFactory::getInstance().listAvailable()) {

    output_test << "Testing " << name << " solver:";
    try {
      // Get specific options section for this solver. Can't just use default
      // "solver" section, as we run into problems when solvers use the same
      // name for an option with inconsistent defaults
      auto options = Options::getRoot()->getSection(name);
      auto solver = std::unique_ptr<Solver>{Solver::create(name, options)};

      TestSolver model{};
      solver->setModel(&model);

      BoutMonitor bout_monitor{};
      solver->addMonitor(&bout_monitor, Solver::BACK);

      solver->solve();

      if (std::abs(model.field(1, 1, 0) - expected) > tolerance) {
        output_test << " FAILED\n";
        errors[name] = model.field(1, 1, 0);
      } else {
        output_test << " PASSED\n";
      }
    } catch (BoutException& e) {
      // Don't let one bad solver stop us trying the rest
      output_test << " ERROR\n";
      output_info << "Error encountered with solver " << name << "\n";
      output_info << e.what() << endl;
      errors[name] = 0.;
    }
  }

  BoutFinalise(false);

  if (!errors.empty()) {
    output_test << "\n => Some failed tests\n\n";
    for (auto& error : errors) {
      output_test << "    " << error.first << " got: " << error.second
                  << ", expected: " << expected << "\n";
    }
  } else {
    output_test << "\n => All tests passed\n";
  }

  return errors.size();
}
