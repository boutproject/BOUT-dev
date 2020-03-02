#include <bout/physicsmodel.hxx>
#include <bout/sys/timer.hxx>

class BoundaryOp_timing : public PhysicsModel {
  int init(bool restarting);
  int rhs(BoutReal UNUSED(t)) {
    ddt(f)=0;
    return 0; }
  Field3D f;
};

int BoundaryOp_timing::init(bool UNUSED(restarting)) {
  SOLVE_FOR(f);

  Options* opt = Options::getRoot();
  int ntests;
  OPTION(opt, ntests, 10000);


  {
    Timer timer("dirichlet_o3");
    for (int i=0; i<ntests; i++) {
      f.applyBoundary();
    }
  }

  std::string name;
  opt->getSection("f")->get("bndry_all",name,"error");

  double time = Timer::getTime("dirichlet_o3");
  printf("%-20s | %12.8f s", name.c_str(), time);
  
  return 0;
}

BOUTMAIN(BoundaryOp_timing);
