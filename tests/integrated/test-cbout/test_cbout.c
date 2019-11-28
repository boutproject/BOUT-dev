#include <cbout/bout.h>

int main(int argc, char** argv) {

  /* Start BOUT++ */
  bout_initialise(argc, argv);
  
  Field3D *f = Field3D_new();

  Field3D_delete(f);
  
  /* Shut down BOUT++ */
  bout_finalise();
  return 0;
}
