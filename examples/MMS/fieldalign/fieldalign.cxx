
#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

class FieldAlign : public PhysicsModel {
protected:
  int init(bool restarting) {
    mesh->get(vy, "vy");
    mesh->get(vz, "vz");
    solver->add(f, "f");
    return 0;
  }
  
  int rhs(BoutReal t) {
    mesh->communicate(f);

    // df/dt = df/dtheta + df/dphi

    ddt(f) = 
      vy * (mesh->g22*DDY(f) + mesh->g23*DDZ(f)) +    // Upwinding with second-order central differencing 
      vz * (mesh->g33*DDZ(f) + mesh->g23*DDY(f))      // (unstable without additional dissipation)
      - SQ(SQ(mesh->dy))*D4DY4(f) - SQ(SQ(mesh->dz))*D4DY4(f);   // Numerical dissipation terms
    
    
    return 0;
  }
  
private:
  Field3D f;
  Field2D vy, vz;
};

BOUTMAIN(FieldAlign);
