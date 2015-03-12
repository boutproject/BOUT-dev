
#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

class FieldAlign : public PhysicsModel {
protected:
  int init(bool restarting) {
    solver->add(f, "f");
    return 0;
  }
  
  int rhs(BoutReal t) {
    mesh->communicate(f);

    // df/dt = df/dtheta + df/dphi
    ddt(f) = sqrt(mesh->g22)*DDY(f) + (mesh->g23/sqrt(mesh->g22) + sqrt(mesh->g33))*DDZ(f) - SQ(SQ(mesh->dy))*D4DY4(f) - SQ(SQ(mesh->dz))*D4DY4(f);
    
    return 0;
  }
  
private:
  Field3D f;
};

BOUTMAIN(FieldAlign);
