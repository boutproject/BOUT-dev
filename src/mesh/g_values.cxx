#include "bout/g_values.hxx"
#include "bout/coordinates.hxx"
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

  bout::FieldMetric Jg12 = J * g12;
  mesh->communicate(Jg12);
  G1_m =
      (coordinates.DDX(J * g11) + coordinates.DDY(Jg12) + coordinates.DDZ(J * g13)) / J;
  bout::FieldMetric Jg22 = J * g22;
  mesh->communicate(Jg22);
  G2_m =
      (coordinates.DDX(J * g12) + coordinates.DDY(Jg22) + coordinates.DDZ(J * g23)) / J;
  bout::FieldMetric Jg23 = J * g23;
  mesh->communicate(Jg23);
  G3_m =
      (coordinates.DDX(J * g13) + coordinates.DDY(Jg23) + coordinates.DDZ(J * g33)) / J;

  mesh->communicate(G1_m, G2_m, G3_m);
}
