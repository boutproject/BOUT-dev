#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

class FieldAlign : public PhysicsModel {
protected:
  int init(bool restarting) {
    mesh->get(vx, "vx");
    mesh->get(vy, "vy");
    mesh->get(vz, "vz");
    mesh->get(G, "G");
    solver->add(f, "f");
    return 0;
  }

  int rhs(BoutReal t) {
    Coordinates *metric = mesh->getCoordinates();
    mesh->communicate(f);
    f.applyBoundary(t);

    // df/dt = df/dtheta + df/dphi

    ddt(f) =
        vx / G * (metric->g11*DDX(f) + metric->g12*DDY(f) + metric->g13*DDZ(f)) +
        vy / G * (metric->g12*DDX(f) + metric->g22*DDY(f) + metric->g23*DDZ(f)) +    // Upwinding with second-order central differencing
        vz / G * (metric->g13*DDX(f) + metric->g23*DDY(f) + metric->g33*DDZ(f));      // (unstable without additional dissipation)
    - SQ(SQ(metric->dx))*D4DX4(f) /*- SQ(SQ(metric->dy))*D4DY4(f)*/ - SQ(SQ(metric->dz))*D4DZ4(f);   // Numerical dissipation terms

    return 0;
  }

private:
  Field3D f;
  Field2D vx, vy, vz, G;
};

BOUTMAIN(FieldAlign);
