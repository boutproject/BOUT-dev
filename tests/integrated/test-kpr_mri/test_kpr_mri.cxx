#include "bout/petsclib.hxx"
#include "bout/physicsmodel.hxx"
#include "bout/solver.hxx"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

// A simple phyics model with a stiff decay towards a steady state solution
//
class TestSolver : public PhysicsModel {
public:
  Field3D f, g;

  BoutReal e  = 0.5;       /* fast/slow coupling strength */
  BoutReal G  = -100.0;    /* stiffness at slow time scale */
  BoutReal w  = 100.0;     /* time-scale separation factor */

  int init(bool UNUSED(restarting)) override {
    solver->add(f, "f");
    solver->add(g, "g");

    f = sqrt(3.0/2.0);
    g = sqrt(3.0);

    return 0;
  }

  int rhs_se(BoutReal t) override {

    ddt(f) = -0.5*sin(t)/(2.0*f);
    ddt(g) = 0.0;

    return 0;
  }

  int rhs_si(BoutReal t) override {
    /* fill in the slow implicit RHS function:
      [G e]*[(-1+u^2-r(t))/(2*u))]
      [0 0] [(-2+v^2-s(t))/(2*v)]  */
    BoutReal tmp1 = (-1.0 + f(0,0,0) * f(0,0,0) - 0.5*cos(t)) / (2.0 * f(0,0,0));
    BoutReal tmp2 = (-2.0 + g(0,0,0) * g(0,0,0) - cos(w*t)) / (2.0 * g(0,0,0));
    ddt(f) = G * tmp1 + e * tmp2;
    ddt(g) = 0.0;

    return 0;
  }

  int rhs_fe(BoutReal t) override {

    ddt(f) = 0.0;
    ddt(g) = 0.0;

    return 0;
  }

  int rhs_fi(BoutReal t) override {

    BoutReal tmp1 = (-1.0 + f(0,0,0) * f(0,0,0) - 0.5*cos(t)) / (2.0 * f(0,0,0));
    BoutReal tmp2 = (-2.0 + g(0,0,0) * g(0,0,0) - cos(w*t)) / (2.0 * g(0,0,0));
    ddt(f) = 0.0;
    ddt(g) = e * tmp1 - tmp2 - w * sin(w*t) / (2.0 * sqrt(2.0 + cos(w * t)));

    return 0;
  }

  int rhs_s(BoutReal t) override {

    BoutReal tmp1 = (-1.0 + f(0,0,0) * f(0,0,0) - 0.5*cos(t)) / (2.0 * f(0,0,0));
    BoutReal tmp2 = (-2.0 + g(0,0,0) * g(0,0,0) - cos(w*t)) / (2.0 * g(0,0,0));
    ddt(f) = G * tmp1 + e * tmp2 - 0.5*sin(t) / (2.0 * f(0,0,0));
    ddt(g) = 0.0;

    return 0;
  }

  int rhs_f(BoutReal t) override {

    BoutReal tmp1 = (-1.0 + f(0,0,0) * f(0,0,0) - 0.5*cos(t)) / (2.0 * f(0,0,0));
    BoutReal tmp2 = (-2.0 + g(0,0,0) * g(0,0,0) - cos(w*t)) / (2.0 * g(0,0,0));
    ddt(f) = 0.0;
    ddt(g) = e * tmp1 - tmp2 - w * sin(w*t) / (2.0 * sqrt(2.0 + cos(w * t)));

    return 0;
  }

  int rhs(BoutReal t) override {

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "f(0,0,0) = " << f(0,0,0) << std::endl;
    std::cout << "g(0,0,0) = " << g(0,0,0) << std::endl;

    BoutReal tmp1 = f(0,0,0);
    BoutReal tmp2 = g(0,0,0);

    tmp1 = (-1.0 + tmp1 * tmp1 - 0.5*cos(t)) / (2.0 * tmp1);
    tmp2 = (-2.0 + tmp2 * tmp2 - cos(w*t)) / (2.0 * tmp2);

    ddt(f) = G * tmp1 + e * tmp2 - 0.5*sin(t) / (2.0 * f(0,0,0));
    ddt(g) = e * tmp1 - tmp2 - w * sin(w*t) / (2.0 * sqrt(2.0 + cos(w * t)));

    return 0;
  }

  bool check_solution(BoutReal atol) {
    // Return true if correct solution
    return (std::abs(f(1, 1, 0)) < atol) and (std::abs(g(1, 1, 0) - 1) < atol);
  }

  BoutReal compute_error(BoutReal t)
  {
    // return (std::max(abs(sqrt(0.5*cos(t) + 1.0) - f(0,0,0)), abs(sqrt(  cos(w*t) + 2.0) - g(0,0,0))));
    
    return sqrt( pow(sqrt(0.5*cos(t) + 1.0) - f(0,0,0), 2.0) + 
                 pow(sqrt(  cos(w*t) + 2.0) - g(0,0,0), 2.0));
  }

  // Don't need any restarting, or options to control data paths
  int postInit(bool) override { return 0; }
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
  root["restart_files"]["enabled"] = false;

  Solver::setArgs(argc, argv);
  BoutComm::setArgs(argc, argv);

  bout::globals::mpi = new MpiWrapper();

  bout::globals::mesh = Mesh::create();
  bout::globals::mesh->load();

  // Global options
  root["nout"] = 100;
  root["timestep"] = 0.01;

  // Get specific options section for this solver. Can't just use default
  // "solver" section, as we run into problems when solvers use the same
  // name for an option with inconsistent defaults
  auto options = Options::getRoot()->getSection("arkode_mri");
  auto solver = std::unique_ptr<Solver>{Solver::create("arkode_mri", options)};
  
  TestSolver model{};
  solver->setModel(&model);

  BoutMonitor bout_monitor{};
  solver->addMonitor(&bout_monitor, Solver::BACK);

  solver->solve();

  BoutReal error = model.compute_error(1.0);

  std::cout << "error = " << error << std::endl;

  if (model.check_solution(tolerance)) {
    output_test << " PASSED\n";
    return 0;
  }
  output_test << " FAILED\n";
  return 1;
}
