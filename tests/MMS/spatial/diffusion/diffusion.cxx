#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <math.h>
#include <bout/constants.hxx>

Field3D N;

BoutReal Dx, Dy, Dz;
BoutReal Lx, Ly, Lz;

int physics_init(bool UNUSED(restarting)) {
  // Get the options
  Options *meshoptions = Options::getRoot()->getSection("mesh");

  Coordinates *coords = mesh->getCoordinates();

  meshoptions->get("Lx",Lx,1.0);
  meshoptions->get("Ly",Ly,1.0);

  /*this assumes equidistant grid*/
  coords->dx = Lx/(mesh->GlobalNx - 2*mesh->xstart);
  
  coords->dy = Ly/(mesh->GlobalNy - 2*mesh->ystart);
  
  output.write("SIZES: {:d}, {:d}, {:e}\n", mesh->GlobalNy, (mesh->GlobalNy - 2*mesh->ystart), coords->dy(0,0,0));

  SAVE_ONCE2(Lx,Ly);

  Options *cytooptions = Options::getRoot()->getSection("cyto");
  OPTION(cytooptions, Dx, 1.0);
  OPTION(cytooptions, Dy, -1.0);
  OPTION(cytooptions, Dz, -1.0);

  SAVE_ONCE3(Dx, Dy, Dz);

  //set mesh
  coords->g11 = 1.0;
  coords->g22 = 1.0;
  coords->g33 = 1.0;
  coords->g12 = 0.0;
  coords->g13 = 0.0;
  coords->g23 = 0.0;

  coords->g_11 = 1.0;
  coords->g_22 = 1.0;
  coords->g_33 = 1.0;
  coords->g_12 = 0.0;
  coords->g_13 = 0.0;
  coords->g_23 = 0.0;
  coords->geometry();

  // Tell BOUT++ to solve N
  SOLVE_FOR(N);

  return 0;
}

int physics_run(BoutReal UNUSED(t)) {
  mesh->communicate(N); // Communicate guard cells

  ddt(N) = 0.0;
  
  if(Dx > 0.0)
    ddt(N) += Dx * D2DX2(N);
  
  if(Dy > 0.0)
    ddt(N) += Dy * D2DY2(N);
  
  if(Dz > 0.0)
    ddt(N) += Dz * D2DZ2(N);
  
  return 0;
}

