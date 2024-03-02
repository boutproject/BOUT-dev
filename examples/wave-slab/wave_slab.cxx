/*
  Simple wave test in a sheared slab domain.
  
  Uses the same field-aligned Clebsch coordinate system as
  most BOUT++ tokamak simulations. See the coordinates manual
  for details.
  
  Note: Here the only components of the coordinate system which 
  are tested are g_22 (for Grad_par), and the twist shift angle.
  
 */

#include <bout/physicsmodel.hxx>

class WaveTest : public PhysicsModel {
public:
  int init(bool UNUSED(restarting)) {
    auto* coords = mesh->getCoordinates();
    Field2D Rxy, Bpxy, Btxy, hthe, I;
    GRID_LOAD(Rxy);
    GRID_LOAD(Bpxy);
    GRID_LOAD(Btxy);
    GRID_LOAD(hthe);
    coords->setBxy(mesh->get("Bxy"));
    int ShiftXderivs = 0;
    mesh->get(ShiftXderivs, "false");
    if (ShiftXderivs) {
      // No integrated shear in metric
      I = 0.0;
    } else {
      mesh->get(I, "sinty");
    }

    MetricTensor::FieldMetric g11, g22, g33, g12, g13, g23;
    g11 = pow(Rxy * Bpxy, 2.0);
    g22 = 1.0 / pow(hthe, 2.0);
    g33 = pow(I, 2.0) * coords->g11() + pow(coords->Bxy(), 2.0) / coords->g11();
    g12 = 0.0;
    g13 = -I * coords->g11();
    g23 = -Btxy / (hthe * Bpxy * Rxy);
    coords->setContravariantMetricTensor(MetricTensor(g11, g22, g33, g12, g13, g23));

    coords->setJ(hthe / Bpxy);

    MetricTensor::FieldMetric g_11, g_22, g_33, g_12, g_13, g_23;
    g_11 = 1.0 / coords->g11() + (pow(I * Rxy, 2.0));
    g_22 = pow(coords->Bxy() * hthe / Bpxy, 2.0);
    g_33 = Rxy * Rxy;
    g_12 = Btxy * hthe * I * Rxy / Bpxy;
    g_13 = I * Rxy * Rxy;
    g_23 = Btxy * hthe * Rxy / Bpxy;
    coords->setCovariantMetricTensor(MetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));

    coords->calculateGeometry();

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
