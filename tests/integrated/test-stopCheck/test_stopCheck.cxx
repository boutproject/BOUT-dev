/*
 * Check stop tests
 * 
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include "unused.hxx"

Field3D N;

int physics_init(bool UNUSED(restarting)) {
  solver->add(N,"N");
  return 0;
}

int physics_run(BoutReal UNUSED(t)) {
  mesh->communicate(N);
  ddt(N) = 0.;
  return 0;
}
