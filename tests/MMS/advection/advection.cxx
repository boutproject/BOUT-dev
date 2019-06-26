/*
 * MMS test of advection operators
 *
 */

#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>
#include <derivs.hxx>

class AdvectMMS : public PhysicsModel {
public:
  int init(bool) {
    solver->add(f, "f");
    
    Options::getRoot()->get("method", method, 0);
    g = FieldFactory::get()->create3D("g:solution", Options::getRoot(), mesh, CELL_CENTRE);

    Coordinates *coords = mesh->getCoordinates();

    dx_sq_sq = SQ(SQ(coords->dx));
    dz_sq_sq = SQ(SQ(coords->dz));
    
    return 0;
  }
  int rhs(BoutReal) {
    mesh->communicate(f);
    ddt(f) = -bracket(g, f, (BRACKET_METHOD) method)
      - 20.*(dx_sq_sq*D4DX4(f) + dz_sq_sq*D4DZ4(f))
      ;

    return 0;
  }
private:
  int method;
  Field3D f,g;
  Field3D dx_sq_sq, dz_sq_sq;
};

BOUTMAIN(AdvectMMS);
