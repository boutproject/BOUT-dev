#include <boutmain.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>
#include <bout/mesh.hxx>
#include <assert.h>

#include <interpolation.hxx>

Field3D n;

void save_diff(Field3D a, Field3D b,CELL_LOC in=CELL_CENTRE,CELL_LOC out=CELL_CENTRE);


Field3D error;
Field3D diff;
Field3D exact;

int physics_init(bool restart) {
  SOLVE_FOR(n);
  auto ff = FieldFactory::get();
  auto oo = Options::getRoot();
  n=0;
  error=0;
  SAVE_REPEAT(error);
  SAVE_REPEAT(diff);
  SAVE_REPEAT(exact);

  int inloc;
  OPTION(oo,inloc,-1);
  int outloc;
  OPTION(oo,outloc,-1);
  
  auto inf = ff->create3D("all:function",oo,mesh,(CELL_LOC)inloc,0);
  //for (auto out : pos ){
  save_diff(interp_to(inf,(CELL_LOC)outloc),
                ff->create3D("all:function",
                             oo,
                             mesh, (CELL_LOC)outloc, 0));

  return 0;
}


int physics_run(BoutReal time) {
  ddt(n)=0;
  return 0;
}
void save_diff(Field3D a, Field3D b,CELL_LOC in,CELL_LOC out){
  exact=b;
  diff=a;
  error= a-b;
  checkData(error);
  dump.write();
}
