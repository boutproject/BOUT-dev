/*
 * Check stop tests
 * 
 */

#include <bout.hxx>
#include <bout/boutmain.hxx>

Field3D N;

int physics_init(bool restarting) {
  solver->add(N,"N");
  return 0;
}

int physics_run(BoutReal t) {
  mesh->communicate(N);
  ddt(N) = 0.;
  return 0;
}
