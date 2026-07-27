#include "bout/g_values.hxx"
#include "bout/coordinates.hxx"
#include "bout/derivs.hxx"
#include "bout/mesh.hxx"
#include "bout/metric_tensor.hxx"

GValues::GValues(const Coordinates& coordinates) {

  const auto& contravariantMetricTensor = coordinates.getContravariantMetricTensor();
  const auto& J = coordinates.J();

  const auto& g11 = contravariantMetricTensor.g11();
  const auto& g22 = contravariantMetricTensor.g22();
  const auto& g33 = contravariantMetricTensor.g33();
  const auto& g12 = contravariantMetricTensor.g12();
  const auto& g13 = contravariantMetricTensor.g13();
  const auto& g23 = contravariantMetricTensor.g23();

  auto* mesh = J.getMesh();

  G1_m = (DDX(J * g11) + DDY(J.asField3DParallel() * g12) + DDZ(J * g13)) / J;
  G2_m = (DDX(J * g12) + DDY(J.asField3DParallel() * g22) + DDZ(J * g23)) / J;
  G3_m = (DDX(J * g13) + DDY(J.asField3DParallel() * g23) + DDZ(J * g33)) / J;

  mesh->communicate_no_slices(G1_m, G2_m, G3_m);
}
