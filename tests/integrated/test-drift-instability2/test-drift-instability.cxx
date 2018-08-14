/* 
 * Test frequency and growth rate of the resistive drift wave instability
 */

#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>

class DW : public PhysicsModel {
protected:
  int init(bool restarting);
  int rhs(BoutReal time);
private:
  Field3D n; /// perturbed density
  Field3D w; /// vorticity
  Field3D phi; /// electrostatic potential
  Field3D n0; /// background density
  Field3D jpar; /// parallel current
  BoutReal mu; /// mass ratio me/mi
  BoutReal eta; /// parallel resistivity
  Laplacian* phi_solver;
  CELL_LOC maybe_ylow;
  BoutReal omega, gamma; /// analytic expressions for mode frequency and growth rate
  BoutReal hyperdiff; /// coefficient for hyperdiffusion, added for numerical stability
};

int DW::init(bool restarting) {
  Options* options = Options::getRoot()->getSection("driftwave");

  OPTION(options, mu, 1./3645.776967344542);
  OPTION(options, eta, .01);
  initial_profile("n0", n0);

  // set the background density gradient
  n0.setBoundary("n0");
  n0.applyBoundary();

  OPTION(options, omega, 0.);
  OPTION(options, gamma, 0.);
  SAVE_ONCE2(omega, gamma);

  OPTION(options, hyperdiff, 1.);
  SAVE_ONCE(hyperdiff);

  output<<"expected omega="<<omega<<endl;
  output<<"expected gamma="<<gamma<<endl;

  if (mesh->StaggerGrids) {
    maybe_ylow = CELL_YLOW;
    jpar.setLocation(CELL_YLOW);
  } else {
    maybe_ylow = CELL_CENTRE;
  }

  SOLVE_FOR2(n, w);
  SAVE_REPEAT2(phi, jpar);
  SAVE_ONCE(n0);

  phi_solver = Laplacian::create();

  return 0;
}

int DW::rhs(BoutReal time) {

  mesh->communicate(n, w);

  phi = phi_solver->solve(w);
  mesh->communicate(phi);

  jpar = (Grad_par(n, maybe_ylow) - Grad_par(phi, maybe_ylow))/mu/eta;
  mesh->communicate(jpar);

  Coordinates* coords = mesh->coordinates();

  ddt(n) = -DDZ(phi)*DDX(n0);

  ddt(w) = Div_par(jpar, CELL_CENTRE) - hyperdiff*(pow(coords->dy, 4)*D4DY4(w) + pow(coords->dz, 4)*D4DZ4(w));

  return 0;
}

BOUTMAIN(DW);

