/*
 * Demonstrates how to use staggered grids with boundary conditions
 */

#include <bout.hxx>
#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

class Test_staggered : public PhysicsModel {
protected:
  int init(bool UNUSED(restart)) override;
  int rhs(BoutReal UNUSED(time)) override;
};


using bout::globals::mesh;

Field3D n, v;
CELL_LOC maybe_ylow{CELL_CENTRE};

int Test_staggered::init(bool UNUSED(restart)) {

  if (mesh->StaggerGrids) {
    maybe_ylow = CELL_YLOW;
  }
  
  v.setLocation(maybe_ylow); // Staggered relative to n
  
  SOLVE_FOR(n, v);

  return 0;
}

int Test_staggered::rhs(BoutReal UNUSED(time)) {
  mesh->communicate(n, v);
  
  //ddt(n) = -Div_par_flux(v, n, CELL_CENTRE);
  ddt(n) = -n*Grad_par(v, CELL_CENTRE) - Vpar_Grad_par(v, n, CELL_CENTRE);
  
  ddt(v) = -Grad_par(n, maybe_ylow);
 
  // Have to manually apply the lower Y boundary region, using a width of 3
  for( RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
    for(int y=2;y>=0;y--) 
      for(int z=0;z<mesh->LocalNz;z++) {
        ddt(v)(rlow.ind,y,z) = ddt(v)(rlow.ind,y+1,z);
      }
  
  return 0;
}



BOUTMAIN(Test_staggered)
