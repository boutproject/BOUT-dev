#include <boutmain.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>

#include <assert.h>

Field3D n;

void compare(Field3D a, Field3D b);

int physics_init(bool restart) {
  n=1;
  int iter;
  OPTION(Options::getRoot(),iter,10);
  for (int i =0; i < iter; ++i){
    mesh->indexDDX(n,CELL_DEFAULT,DIFF_DEFAULT);
    mesh->indexDDY(n,CELL_DEFAULT,DIFF_DEFAULT);
    mesh->indexDDZ(n,CELL_DEFAULT,DIFF_DEFAULT,true);
    mesh->indexD2DX2(n,CELL_DEFAULT,DIFF_DEFAULT);
    mesh->indexD2DY2(n,CELL_DEFAULT,DIFF_DEFAULT);
    mesh->indexD2DZ2(n,CELL_DEFAULT,DIFF_DEFAULT,true);
  }
  SOLVE_FOR(n);
  return 0;
}


int physics_run(BoutReal time) {
  ddt(n)=0;
  return 0;
}
