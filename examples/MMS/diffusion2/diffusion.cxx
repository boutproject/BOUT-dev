#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <math.h>
#include <bout/constants.hxx>

Field3D N;

BoutReal mu_N; // Parallel collisional diffusion coefficient
BoutReal Lx, Ly, Lz;

int physics_init(bool restarting) {
  // Get the options
  Options *meshoptions = Options::getRoot()->getSection("mesh");

  meshoptions->get("Lx",Lx,1.0);
  meshoptions->get("Ly",Ly,1.0);

  /*this assumes equidistant grid*/
  int nguard = mesh->xstart;
  mesh->dx = Lx/(mesh->GlobalNx - 2*nguard);
  mesh->dy = Ly/(mesh->GlobalNy - 2*nguard);

  SAVE_ONCE2(Lx,Ly);

  Options *cytooptions = Options::getRoot()->getSection("cyto");
  cytooptions->get("dis", mu_N, 1);

  SAVE_ONCE(mu_N);

  //set mesh
  mesh->g11 = 1.0;
  mesh->g22 = 1.0;
  mesh->g33 = 1.0;
  mesh->g12 = 0.0;
  mesh->g13 = 0.0;
  mesh->g23 = 0.0;

  mesh->g_11 = 1.0;
  mesh->g_22 = 1.0;
  mesh->g_33 = 1.0;
  mesh->g_12 = 0.0;
  mesh->g_13 = 0.0;
  mesh->g_23 = 0.0;
  mesh->geometry();

  // Tell BOUT++ to solve N
  SOLVE_FOR(N);

  return 0;
}

int physics_run(BoutReal t) {
  mesh->communicate(N); // Communicate guard cells

  //update time-dependent boundary conditions
  output.write("APPLYING BOUNDARY\n");
  N.applyBoundary(t);

  ddt(N) = mu_N* D2DX2(N);
  
  output.write("FINISHED RHS\n");
  return 0;
}

