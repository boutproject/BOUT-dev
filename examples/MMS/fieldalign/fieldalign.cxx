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
    Coordinates *coords = mesh->coordinates();

    // df/dt = df/dtheta + df/dphi

    ddt(f) =
      vy * (coords->g22*DDY(f) + coords->g23*DDZ(f)) +    // Upwinding with second-order central differencing
      vz * (coords->g33*DDZ(f) + coords->g23*DDY(f))      // (unstable without additional dissipation)
      /*- SQ(SQ(coords->dy))*D4DY4(f)*/ - SQ(SQ(coords->dz))*D4DZ4(f);   // Numerical dissipation terms

    return 0;
  }

private:
  Field3D f;
  Field2D vy, vz;
};

BOUTMAIN(FieldAlign);
