/*
  Simple wave test in a sheared slab domain.
  
  Uses the same field-aligned Clebsch coordinate system as
  most BOUT++ tokamak simulations. See the coordinates manual
  for details.
  
  Note: Here the only components of the coordinate system which 
  are tested are g_22 (for Grad_par), and the twist shift angle.
  
 */

#include <bout/physicsmodel.hxx>
#include <bout/tokamak_coordinates_factory.hxx>

class WaveTest : public PhysicsModel {
public:
  int init(bool UNUSED(restarting)) {
    Field2D Rxy, Bpxy, Btxy, hthe, I;
    GRID_LOAD(Rxy);
    GRID_LOAD(Bpxy);
    GRID_LOAD(Btxy);
    GRID_LOAD(hthe);
    const auto& Bxy = mesh->get("Bxy");
    int ShiftXderivs = 0;
    mesh->get(ShiftXderivs, "false");
    if (ShiftXderivs) {
      // No integrated shear in metric
      I = 0.0;
    } else {
      mesh->get(I, "sinty");
    }

    const auto tokamak_coordinates_factory = TokamakCoordinatesFactory(*mesh, Rxy, Bpxy, Btxy, Bxy);
    const auto& coords = tokamak_coordinates_factory.make_tokamak_coordinates(hthe, I);

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
