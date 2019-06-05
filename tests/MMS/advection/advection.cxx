/*
 * MMS test of advection operators
 *
 */

#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>
#include <derivs.hxx>

class AdvectMMS : public PhysicsModel {
public:
  int init(bool UNUSED(restarting)) {
    solver->add(f, "f");
    
    Options::getRoot()->get("method", method, 0);
    
    return 0;
  }
  int rhs(BoutReal time) {
    mesh->communicate(f);
    Coordinates *coords = mesh->getCoordinates();
    
    g = FieldFactory::get()->create3D("g:solution", Options::getRoot(), mesh, CELL_CENTRE, time);
    
    ddt(f) = -bracket(g, f, (BRACKET_METHOD) method)
      - 20.*(SQ(SQ(coords->dx))*D4DX4(f) + SQ(SQ(coords->dz))*D4DZ4(f))
      //+ 20.*(SQ(coords->dx)*D2DX2(f) + SQ(coords->dz)*D2DZ2(f))
      ;

    return 0;
  }
private:
  int method;
  Field3D f,g;
};

BOUTMAIN(AdvectMMS);
