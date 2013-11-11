/*
 * Global fields for gather/scatter
 * 
 * 
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include <bout/globalfield.hxx>

int physics_init(bool restarting) {
  
  GlobalField2D g2d(mesh);
 
  return 1; // Signal an error, so quits
}

int physics_run(BoutReal t) {
  // Doesn't do anything
  return 1;
}
