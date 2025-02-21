/*
  Simple wave test in a sheared slab domain.
  
  Uses the same field-aligned Clebsch coordinate system as
  most BOUT++ tokamak simulations. See the coordinates manual
  for details.
  
  Note: Here the only components of the coordinate system which 
  are tested are g_22 (for Grad_par), and the twist shift angle.
  
 */

#include <bout/physicsmodel.hxx>
#include <bout/tokamak_coordinates.hxx>

class WaveTest : public PhysicsModel {
public:
  int init(bool UNUSED(restarting)) {

    auto tokamak_options = TokamakOptions(*mesh);

    int ShiftXderivs = 0;
    mesh->get(ShiftXderivs, "false");
    BoutReal shearFactor = 1.0;
    if (ShiftXderivs) {
      // No integrated shear in metric
      shearFactor = 0.0;
    }
    set_tokamak_coordinates_on_mesh(tokamak_options, *mesh, 1.0, 1.0, shearFactor);
    
    auto* coords = mesh->getCoordinates();

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
