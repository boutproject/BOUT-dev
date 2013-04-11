#include <bout.hxx>
#include <boutmain.hxx>

Field3D n;
BoutReal D;

int physics_init(bool restarting) {
  
  Options *opt = Options::getRoot()->getSection("diffusion");
  OPTION(opt, D, 0.1);
  
  SOLVE_FOR(n);
  
  return 0;
}

int physics_run(BoutReal t) {
  mesh->communicate(n);
  
  ddt(n) = D * Delp2(n);
  
  return 0;
}
