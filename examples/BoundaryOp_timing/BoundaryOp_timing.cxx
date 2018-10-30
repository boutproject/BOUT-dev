#include <bout/physicsmodel.hxx>
#include <bout/sys/timer.hxx>

class BoundaryOp_timing : public PhysicsModel {
  int init(bool restarting);
  int rhs(BoutReal UNUSED(t)) { return 1; }
  Field3D f_dirichlet, f_neumann, f_dirichlet_o3;
};

int BoundaryOp_timing::init(bool UNUSED(restarting)) {
  SOLVE_FOR3(f_dirichlet, f_neumann, f_dirichlet_o3);

  Options* opt = Options::getRoot();
  int ntests;
  OPTION(opt, ntests, 10000);


  // test Dirichlet_O3

  {
    Timer timer("dirichlet_o3");
    for (int i=0; i<ntests; i++) {
      f_dirichlet_o3.applyBoundary();
    }
  }


  // test Dirichlet

  {
    Timer timer("dirichlet");
    for (int i=0; i<ntests; i++) {
      f_dirichlet.applyBoundary();
    }
  }


  // test Neumann

  {
    Timer timer("neumann");
    for (int i=0; i<ntests; i++) {
      f_neumann.applyBoundary();
    }
  }


  output<<"dirichlet:\t"<<ntests<<" iterations took "<<Timer::getTime("dirichlet")<<" s"<<endl;
  output<<"neumann:\t"<<ntests<<" iterations took "<<Timer::getTime("neumann")<<" s"<<endl;
  output<<"dirichlet_o3:\t"<<ntests<<" iterations took "<<Timer::getTime("dirichlet_o3")<<" s"<<endl;

  return 1;
}

BOUTMAIN(BoundaryOp_timing);
