#include <bout/physicsmodel.hxx>
#include <bout/sys/timer.hxx>

class BoundaryOp_timing : public PhysicsModel {
  int init(bool restarting);
  int rhs(BoutReal UNUSED(t)) { return 1; }
  Field3D f_dirichlet, f_neumann, f_dirichlet_o3;
  Field3D f_dirichlet_val, f_neumann_val, f_dirichlet_o3_val;
  Field3D f_dirichlet_expr, f_neumann_expr, f_dirichlet_o3_expr;
};

int BoundaryOp_timing::init(bool UNUSED(restarting)) {
  SOLVE_FOR3(f_dirichlet, f_neumann, f_dirichlet_o3);
  SOLVE_FOR3(f_dirichlet_val, f_neumann_val, f_dirichlet_o3_val);
  SOLVE_FOR3(f_dirichlet_expr, f_neumann_expr, f_dirichlet_o3_expr);

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


  // test Dirichlet_O3

  {
    Timer timer("dirichlet_o3_val");
    for (int i=0; i<ntests; i++) {
      f_dirichlet_o3_val.applyBoundary();
    }
  }


  // test Dirichlet

  {
    Timer timer("dirichlet_val");
    for (int i=0; i<ntests; i++) {
      f_dirichlet_val.applyBoundary();
    }
  }


  // test Neumann

  {
    Timer timer("neumann_val");
    for (int i=0; i<ntests; i++) {
      f_neumann_val.applyBoundary();
    }
  }


  // test Dirichlet_O3

  {
    Timer timer("dirichlet_o3_expr");
    for (int i=0; i<ntests; i++) {
      f_dirichlet_o3_expr.applyBoundary();
    }
  }


  // test Dirichlet

  {
    Timer timer("dirichlet_expr");
    for (int i=0; i<ntests; i++) {
      f_dirichlet_expr.applyBoundary();
    }
  }


  // test Neumann

  {
    Timer timer("neumann_expr");
    for (int i=0; i<ntests; i++) {
      f_neumann_expr.applyBoundary();
    }
  }


  output<<"dirichlet:\t"<<ntests<<" iterations took "<<Timer::getTime("dirichlet")<<" s"<<endl;
  output<<"neumann:\t"<<ntests<<" iterations took "<<Timer::getTime("neumann")<<" s"<<endl;
  output<<"dirichlet_o3:\t"<<ntests<<" iterations took "<<Timer::getTime("dirichlet_o3")<<" s"<<endl;
  output<<endl;

  output<<"dirichlet_val:\t"<<ntests<<" iterations took "<<Timer::getTime("dirichlet_val")<<" s"<<endl;
  output<<"neumann_val:\t"<<ntests<<" iterations took "<<Timer::getTime("neumann_val")<<" s"<<endl;
  output<<"dirichlet_o3_val:\t"<<ntests<<" iterations took "<<Timer::getTime("dirichlet_o3_val")<<" s"<<endl;
  output<<endl;

  output<<"dirichlet_expr:\t"<<ntests<<" iterations took "<<Timer::getTime("dirichlet_expr")<<" s"<<endl;
  output<<"neumann_expr:\t"<<ntests<<" iterations took "<<Timer::getTime("neumann_expr")<<" s"<<endl;
  output<<"dirichlet_o3_expr:\t"<<ntests<<" iterations took "<<Timer::getTime("dirichlet_o3_expr")<<" s"<<endl;

  BoutFinalise();
  exit(0);

  return 1;
}

BOUTMAIN(BoundaryOp_timing);
