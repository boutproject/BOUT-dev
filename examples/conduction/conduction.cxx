/*
 * 1D temperature conduction problem
 * 
 */

#include <bout.hxx>
#include <boutmain.hxx>

Field3D T; // Evolving temperature equation only

BoutReal chi; // Parallel conduction coefficient

int physics_init(bool restarting) {
  
  // Get the options
  Options *options = Options::getRoot()->getSection("conduction");
  
  OPTION(options, chi, 1.0); // Read from BOUT.inp, setting default to 1.0

  // Tell BOUT++ to solve T
  SOLVE_FOR(T);
  
  return 0;
}

int physics_run(BoutReal t) {
  mesh->communicate(T); // Communicate guard cells
  
  ddt(T) = Div_par_K_Grad_par(chi, T); // Parallel diffusion Div_{||}( chi * Grad_{||}(T) )
  
  return 0;
}
