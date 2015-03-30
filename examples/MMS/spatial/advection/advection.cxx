/*
 * MMS test of advection operators
 *
 */

#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>
#include <derivs.hxx>

class AdvectMMS : public PhysicsModel {
public:
  int init(bool restarting) {
    solver->add(f, "f");
    
    Options::getRoot()->get("method", method, 0);
    
    return 0;
  }
  int rhs(BoutReal time) {
    mesh->communicate(f);
    
    g = FieldFactory::get()->create3D("g:solution", Options::getRoot(), mesh, CELL_CENTRE, time);
    
    ddt(f) = -bracket(g, f, (BRACKET_METHOD) method)
      - 20.*(SQ(SQ(mesh->dx))*D4DX4(f) + SQ(SQ(mesh->dz))*D4DZ4(f))
      //+ 20.*(SQ(mesh->dx)*D2DX2(f) + SQ(mesh->dz)*D2DZ2(f))
      ;

    return 0;
  }
private:
  int method;
  Field3D f,g;
};

BOUTMAIN(AdvectMMS);
