#include <bout/derivs.hxx>
#include <bout/physicsmodel.hxx>

class FieldAlign : public PhysicsModel {
protected:
  int init(bool UNUSED(restarting)) {
    mesh->get(vx, "vx");
    mesh->get(vy, "vy");
    mesh->get(vz, "vz");
    mesh->get(G, "G");
    solver->add(f, "f");
    return 0;
  }

  int rhs(BoutReal t) {
    Coordinates* metric = mesh->getCoordinates();
    mesh->communicate(f);
    f.applyBoundary(t);

    // df/dt = df/dtheta + df/dphi
    
    Coordinates::MetricTensor metric_tensor = metric->getContravariantMetricTensor();
    ddt(f) =
        vx / G * (metric_tensor.g11 * DDX(f) + metric_tensor.g12 * DDY(f) + metric_tensor.g13 * DDZ(f))
        + vy / G * (metric_tensor.g12 * DDX(f) + metric_tensor.g22 * DDY(f) + metric_tensor.g23 * DDZ(f))
        + // Upwinding with second-order central differencing
        vz / G
            * (metric_tensor.g13 * DDX(f) + metric_tensor.g23 * DDY(f)
               + metric_tensor.g33 * DDZ(f));  // (unstable without additional dissipation)
    -SQ(SQ(metric->dx)) * D4DX4(f)       /*- SQ(SQ(metric->dy))*D4DY4(f)*/
        - SQ(SQ(metric->dz)) * D4DZ4(f); // Numerical dissipation terms

    return 0;
  }

private:
  Field3D f;
  Field2D vx, vy, vz, G;
};

BOUTMAIN(FieldAlign);
