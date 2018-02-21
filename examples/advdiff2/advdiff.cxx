/*******************************************************************
 * Advection-Diffusion Example
 *
 * MVU 19-july-2011
 *******************************************************************/

#include <bout.hxx>
#include <bout/derivs.hxx>
#include <bout/boutmain.hxx>

// Evolving variables 
Field3D V;


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



int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(V);


  //ddt(V) = D2DX2(V) + 0.5*DDX(V) + D2DY2(V);
  ddt(V) = DDX(V);

  
  return 0;
}
