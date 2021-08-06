#include "bout/constants.hxx"
#include "bout/petsclib.hxx"
#include "bout/physicsmodel.hxx"
#include "bout/slepclib.hxx"
#include "bout/solverfactory.hxx"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

// A simple phyics model with a stiff decay towards a steady state solution
//
class TestSolver : public PhysicsModel {
public:
  Field3D f, g;

  int init(bool UNUSED(restarting)) override {
    solver->add(f, "f");
    solver->add(g, "g");

    f = 1.0;
    g = 0.0;

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    // This should have a steady-state solution f=0, g=1
    ddt(f) = 998 * f + 1998 * (g - 1.0);
    ddt(g) = -999 * f - 1999 * (g - 1.0);

    return 0;
  }

  bool check_solution(BoutReal atol) {
    // Return true if correct solution
    return (std::abs(f(1, 1, 0)) < atol) and (std::abs(g(1, 1, 0) - 1) < atol);
  }
};

int main(int argc, char** argv) {

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
  root["restart"]["enabled"] = false;

  PetscLib::setArgs(argc, argv);
  Solver::setArgs(argc, argv);
  BoutComm::setArgs(argc, argv);

  bout::globals::mesh = Mesh::create();
  bout::globals::mesh->load();

  bout::globals::dump =
      bout::experimental::setupDumpFile(Options::root(), *bout::globals::mesh, ".");

  // Global options
  root["NOUT"] = 20;
  root["TIMESTEP"] = 1;

  // Get specific options section for this solver. Can't just use default
  // "solver" section, as we run into problems when solvers use the same
  // name for an option with inconsistent defaults
  auto options = Options::getRoot()->getSection("beuler");
  auto solver = std::unique_ptr<Solver>{Solver::create("beuler", options)};

  TestSolver model{};
  solver->setModel(&model);

  BoutMonitor bout_monitor{};
  solver->addMonitor(&bout_monitor, Solver::BACK);

  solver->solve();

  BoutFinalise(false);

  if (model.check_solution(tolerance)) {
    output_test << " PASSED\n";
    return 0;
  }
  output_test << " FAILED\n";
  return 1;
}
