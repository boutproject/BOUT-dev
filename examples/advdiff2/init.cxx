#include <bout.hxx>
#include <bout/derivs.hxx>
#include <bout/boutmain.hxx>
#include "bout/globals.hxx"


int physics_init(bool restarting)
{
  // 2D initial profiles
  Field2D V0;


  // Read initial conditions

  mesh->get(V0, "V0");
  mesh->get(mesh->dx,   "dx");
  mesh->get(mesh->dy,   "dy");


  // read options


  // Set evolving variables
  bout_solve(V, "V");

  
  if(!restarting) {
    // Set variables to these values (+ the initial perturbation)
    // NOTE: This must be after the calls to bout_solve
    V += V0;
  }
  
  return 0;
}

