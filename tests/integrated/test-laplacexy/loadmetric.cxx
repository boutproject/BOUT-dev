#include "bout/field2d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"
#include "bout/utils.hxx"

#include "loadmetric.hxx"

void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
  // Load metric coefficients from the mesh
  Field2D Rxy, Bpxy, Btxy, hthe, sinty;

  auto mesh = bout::globals::mesh;
  auto coords = mesh->getCoordinates();

  GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics

  // Checking for dpsi and qinty used in BOUT grids
  Field2D dx;
  if (!mesh->get(dx, "dpsi")) {
    output << "\tUsing dpsi as the x grid spacing\n";
    coords->setDx(dx); // Only use dpsi if found
  } else {
    // dx will have been read already from the grid
    output << "\tUsing dx as the x grid spacing\n";
  }
  Field2D qinty;

  Rxy /= Lnorm;
  hthe /= Lnorm;
  sinty *= SQ(Lnorm) * Bnorm;
  coords->setDx(coords->dx() / (SQ(Lnorm) * Bnorm));

  Bpxy /= Bnorm;
  Btxy /= Bnorm;
  coords->setBxy(coords->Bxy() / Bnorm);

  // Calculate metric components
  std::string ptstr;
  Options::getRoot()->get("mesh:paralleltransform", ptstr, "identity");
  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);
  if (ptstr == "shifted") {
    sinty = 0.0; // I disappears from metric
  }

  BoutReal sbp = 1.0; // Sign of Bp
  if (min(Bpxy, true) < 0.0) {
    sbp = -1.0;
  }

  const auto g11 = pow(Rxy * Bpxy, 2);
  const auto g22 = 1.0 / pow(hthe, 2);
  const auto g33 = pow(sinty, 2) * coords->g11() + pow(coords->Bxy(), 2) / coords->g11();
  const auto g12 = 0.0;
  const auto g13 = -sinty * coords->g11();
  const auto g23 = -sbp * Btxy / (hthe * Bpxy * Rxy);
  coords->setContravariantMetricTensor(
      ContravariantMetricTensor(g11, g22, g33, g12, g13, g23));

  coords->setJ(hthe / Bpxy);

  const auto g_11 = 1.0 / coords->g11() + pow(sinty * Rxy, 2);
  const auto g_22 = pow(coords->Bxy() * hthe / Bpxy, 2);
  const auto g_33 = Rxy * Rxy;
  const auto g_12 = sbp * Btxy * hthe * sinty * Rxy / Bpxy;
  const auto g_13 = sinty * Rxy * Rxy;
  const auto g_23 = sbp * Btxy * hthe * Rxy / Bpxy;
  coords->setCovariantMetricTensor(
      CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));
}
