/*
 * MMS test in tokamak geometry
 *
 * Independently tests the bracket operator, perpendicular diffusion,
 * and parallel diffusion operators.
 * 
 */

#include <bout/derivs.hxx>
#include <bout/field_factory.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/tokamak_coordinates.hxx>

class TokamakMMS : public PhysicsModel {
public:
  int init(bool UNUSED(restarting)) {
    solver->add(laplacepar, "laplacepar");
    solver->add(delp2, "delp2");
    solver->add(advect, "advect");

    // Load the metric tensor
    LoadMetric(1.0, 1.0);

    return 0;
  }
  int rhs(BoutReal time) {
    mesh->communicate(advect, delp2, laplacepar);

    drive = FieldFactory::get()->create3D("drive:solution", Options::getRoot(), mesh,
                                          CELL_CENTRE, time);

    // Test bracket advection operator
    ddt(advect) = -1e-3 * bracket(drive, advect, BRACKET_ARAKAWA)
                  - 10.
                        * (SQ(SQ(mesh->getCoordinates()->dx())) * D4DX4(advect)
                           + SQ(SQ(mesh->getCoordinates()->dz())) * D4DZ4(advect));

    // Test perpendicular diffusion operator
    ddt(delp2) = 1e-5 * Delp2(delp2);

    // Test parallel diffusion operator
    ddt(laplacepar) = Laplace_par(laplacepar);

    return 0;
  }
  void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {

    auto tokamak_options = TokamakOptions(*mesh);
    set_tokamak_coordinates_on_mesh(tokamak_options, *mesh, true, Lnorm, Bnorm);

    Coordinates* coords = mesh->getCoordinates();
  }

private:
  Field3D drive;
  Field3D laplacepar, delp2, advect;
};

BOUTMAIN(TokamakMMS);
