/*
 * Check for backtrace after throw 
 */

#include <bout.hxx>
#include <boutmain.hxx>


int physics_init(bool restarting) {
  BoutReal a=1;
  BoutReal b=0;
  BoutReal c = a/b;
  output.write("c is %f\n",c);
  
  BoutException("Tomatoes are red?\n");
  
  return 1;
}

int physics_run(BoutReal time) {
  
  
  return 1;
}
