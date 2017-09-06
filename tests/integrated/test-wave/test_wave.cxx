/*
 * Parallel wave propagation test, using different methods
 *
 * solving 
 * 
 * f' = d/dy(g)
 * g' = d/dy(f)
 * 
 */

#include <bout.hxx>
#include <boutmain.hxx>

Field3D f, g, f2, g2;

BoutReal v;

int physics_init(bool restarting) {
  
  SOLVE_FOR4(f, g, f2, g2);
  
  v = 1.;
  
  return 0;
}

int physics_run(BoutReal time) {
  
  mesh->communicate(f,g,f2,g2);

  // Standard differencing method
  ddt(f) = v*Grad_par(g);
  ddt(g) = v*Grad_par(f);
  
  // MUSCL scheme
  ddt(f2) = v*Grad_par(g2, f2, 1.); // final argument is fastest velocity
  ddt(g2) = v*Grad_par(f2, g2, 1.);
  
  return 0;
}
