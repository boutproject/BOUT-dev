#include <boutmain.hxx>

Field3D n;

// functions I need in subprojects
#include "fuu/fuu.hxx"
void bar();

int physics_init(bool restart) {
  n=1;
  SOLVE_FOR(n);
  // test functions
  fuu();
  bar();
  return 0;
}

int physics_run(BoutReal time) {
  ddt(n)=0;
  return 0;
}
