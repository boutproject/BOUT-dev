/*
  Simple wave test in a sheared slab domain.

  Uses the same field-aligned Clebsch coordinate system as
  most BOUT++ tokamak simulations. See the coordinates manual
  for details.

  Note: Here the only components of the coordinate system which
  are tested are g_22 (for Grad_par), and the twist shift angle.

 */

#include "bout/tokamak_coordinates.hxx"
#include <bout/physicsmodel.hxx>

class WaveTest : public PhysicsModel {
public:
  int init(bool UNUSED(restarting)) {
    bout::set_tokamak_coordinates(*mesh, 1.0, 1.0, true);
    solver->add(f, "f");
    solver->add(g, "g");

    return 0;
  }
  int rhs(BoutReal UNUSED(time)) {
    mesh->communicate(f, g);

    ddt(f) = Grad_par(g);
    ddt(g) = Grad_par(f);

    return 0;
  }

private:
  Field3D f, g;
};

BOUTMAIN(WaveTest);
