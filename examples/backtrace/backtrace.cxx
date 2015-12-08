/*
 * Check for backtrace after throw 
 */

#include <bout.hxx>
#include <boutmain.hxx>


int physics_init(bool restarting) {

  BoutException("Tomatoes are red?\n");
  
  return 1;
}

int physics_run(BoutReal time) {
  
  
  return 1;
}
