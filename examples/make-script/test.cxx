/*
 */

#include <bout.hxx>
#include <boutmain.hxx>

Field2D f;

int physics_init(bool restarting) {
  SOLVE_FOR(f);
  return 0;
}

int physics_run(BoutReal t) {
  ddt(f) = -f;
  return 0;
}
