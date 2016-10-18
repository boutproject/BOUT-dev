/*
 * Test Cyclic Reduction parallel solver
 *
 */

#include <bout.hxx>
#include <boutmain.hxx>

//#include <fielditerator.hxx>
#include <assert.h>

#include <unistd.h>

Field3D n;


int physics_init(bool restarting) {

  /**********************************************
   * * *                                    * * *
   * * *          Field3D test              * * *
   * * *                                    * * *
   **********************************************/
  
  Field3D   d3 = 1.;
  Field2D   d2 = 1.;
  FieldPerp dp = sliceXZ(d3,mesh->xstart);

  #include "test3d.cxx"
  
  /**********************************************
   * * *                                    * * *
   * * *          Field2D test              * * *
   * * *                                    * * *
   **********************************************/


  SOLVE_FOR(n);
  return 0;
}

int physics_run(BoutReal t) {
  // Doesn't do anything
  ddt(n)=0;
  return 0;
}
