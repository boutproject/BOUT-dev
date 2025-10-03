#include "bout/sundials_backports.hxx"

#if SUNDIALS_VERSION_AT_LEAST(7, 2, 0)

#include "bout/physicsmodel.hxx"
#include "bout/solver.hxx"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

// A simple phyics model with an analytical solution
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

    setSplitOperatorMRI();

    return 0;
  }

  int rhs_se(BoutReal t) override {
  /* fill in the slow explicit RHS function:
     [-0.5*sin(t)/(2*f)]
     [      0          ] */
    ddt(f) = -0.5*sin(t)/(2.0*f(1,1,0));
    ddt(g) = 0.0;

    return 0;
  }

  int rhs_si(BoutReal t) override {
    /* fill in the slow implicit RHS function:
      [G e]*[(-1+f^2-0.5*cos(t))/(2*f)]
      [0 0] [(-2+g^2-cos(w*t))/(2*g)  ]  */
    BoutReal tmp1 = (-1.0 + f(1,1,0) * f(1,1,0) - 0.5*cos(t)) / (2.0 * f(1,1,0));
    BoutReal tmp2 = (-2.0 + g(1,1,0) * g(1,1,0) - cos(w*t)) / (2.0 * g(1,1,0));
    ddt(f) = G * tmp1 + e * tmp2;
    ddt(g) = 0.0;

    return 0;
  }

  int rhs_fe(BoutReal t) override {
  /* fill in the fast explicit RHS function:
     [0  0]*[(-1+f^2-0.5*cos(t))/(2*f)] + [         0                      ]
     [e -1] [(-2+g^2-cos(w*t))/(2*g)  ]   [-w*sin(w*t)/(2*sqrt(2+cos(w*t)))] */
    BoutReal tmp1 = (-1.0 + f(1,1,0) * f(1,1,0) - 0.5*cos(t)) / (2.0 * f(1,1,0));
    BoutReal tmp2 = (-2.0 + g(1,1,0) * g(1,1,0) - cos(w*t)) / (2.0 * g(1,1,0));
    ddt(f) = 0.0;
    ddt(g) = e * tmp1 - tmp2 - w * sin(w*t) / (2.0 * sqrt(2.0 + cos(w * t)));

    return 0;
  }

  int rhs_fi(BoutReal UNUSED(t)) override {

    ddt(f) = 0.0;
    ddt(g) = 0.0;

    return 0;
  }

  BoutReal compute_error(BoutReal t)
  {
    /* Compute the error with the true solution:
     f(t) = sqrt(0.5*cos(t) + 1.0)
     g(t) = sqrt(cos(w*t) + 2.0) */
    return sqrt( pow(sqrt(0.5*cos(t) + 1.0) - f(1,1,0), 2.0) +
                 pow(sqrt(cos(w*t) + 2.0) - g(1,1,0), 2.0));
  }

  // Don't need any restarting, or options to control data paths
  int postInit(bool) override { return 0; }
};

int main(int argc, char** argv) {
  // Relative and Absolute tolerances
  constexpr BoutReal rtol = 1.e-10;
  constexpr BoutReal atol = 1.e-12;

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

  std::string sunsolver = "arkode_mri";
  int nout = 100;
  BoutReal timestep = 0.05;
  BoutReal finaltime = nout*timestep;

  // Global options
  root["nout"] = nout;
  root["timestep"] = timestep;
  root[sunsolver]["rtol"] = rtol;
  root[sunsolver]["atol"] = atol;

  // Get specific options section for this solver. Can't just use default
  // "solver" section, as we run into problems when solvers use the same
  // name for an option with inconsistent defaults
  auto options = Options::getRoot()->getSection(sunsolver);
  auto solver = std::unique_ptr<Solver>{Solver::create(sunsolver, options)};

  TestSolver model{};
  solver->setModel(&model);

  BoutMonitor bout_monitor{};
  solver->addMonitor(&bout_monitor, Solver::BACK);

  solver->solve();

  BoutReal error = model.compute_error(finaltime);

  std::cout << "error = " << error << std::endl;

  return 0;
}
#else
// ARKODE-MRI option for BOUT++
// is not available with this configuration
// do nothing for this test
int main(int argc, char** argv) {
  return 0;
}
#endif