#include <bout.hxx>
#include <bout/derivs.hxx>
#include "bout/globals.hxx"


int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(V);


  //ddt(V) = D2DX2(V) + 0.5*DDX(V) + D2DY2(V);
  ddt(V) = DDX(V);

  
  return 0;
}
