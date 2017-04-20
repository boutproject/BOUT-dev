#include <boutmain.hxx>

Field3D n;

int physics_init(bool restart) {
  SOLVE_FOR(n); 
  output.write("Everything is fine");
  return 0;
}


int physics_run(BoutReal time) {
  ddt(n)=0;
  return 0;
}
