#include <bout.hxx>
#include <bout/physicsmodel.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <cmath>
#include <bout/constants.hxx>
#include <unused.hxx>

class Diffusion : public PhysicsModel {
protected:
  int init(bool UNUSED(restarting)) override;
  int rhs(BoutReal t) override;
};


using bout::globals::mesh;

Field3D N;

BoutReal mu_N; // Parallel collisional diffusion coefficient
BoutReal Lx, Ly, Lz;

Coordinates *coord;

int Diffusion::init(bool UNUSED(restarting)) {
  // Get the options
  Options *meshoptions = Options::getRoot()->getSection("mesh");

  coord = mesh->getCoordinates();
  
  meshoptions->get("Lx",Lx,1.0);
  meshoptions->get("Ly",Ly,1.0);

  /*this assumes equidistant grid*/
  int nguard = mesh->xstart;
  coord->dx = Lx/(mesh->GlobalNx - 2*nguard);
  coord->dy = Ly/(mesh->GlobalNy - 2*nguard);

  SAVE_ONCE2(Lx,Ly);

  Options *cytooptions = Options::getRoot()->getSection("cyto");
  cytooptions->get("dis", mu_N, 1);

  SAVE_ONCE(mu_N);

  //set mesh
  coord->g11 = 1.0;
  coord->g22 = 1.0;
  coord->g33 = 1.0;
  coord->g12 = 0.0;
  coord->g13 = 0.0;
  coord->g23 = 0.0;

  coord->g_11 = 1.0;
  coord->g_22 = 1.0;
  coord->g_33 = 1.0;
  coord->g_12 = 0.0;
  coord->g_13 = 0.0;
  coord->g_23 = 0.0;
  coord->geometry();

  // Tell BOUT++ to solve N
  SOLVE_FOR(N);

  return 0;
}

int Diffusion::rhs(BoutReal t) {
  mesh->communicate(N); // Communicate guard cells

  N.applyBoundary(t);

  ddt(N) = mu_N* D2DX2(N);

  return 0;
}

BOUTMAIN(Diffusion)
