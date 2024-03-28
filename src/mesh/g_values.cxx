
#include "bout/g_values.hxx"
#include "bout/coordinates.hxx"

GValues::GValues(FieldMetric G1, FieldMetric G2, FieldMetric G3)
    : G1_(std::move(G1)), G2_(std::move(G2)), G3_(std::move(G3)){};

GValues::GValues(const Coordinates& coordinates) {

  const auto& contravariantMetricTensor = coordinates.getContravariantMetricTensor();
  const auto& J = coordinates.J();

  const auto& g11 = contravariantMetricTensor.g11();
  const auto& g22 = contravariantMetricTensor.g22();
  const auto& g33 = contravariantMetricTensor.g33();
  const auto& g12 = contravariantMetricTensor.g12();
  const auto& g13 = contravariantMetricTensor.g13();
  const auto& g23 = contravariantMetricTensor.g23();

  auto tmp = J * g12;
  coordinates.communicate(tmp);
  setG1((coordinates.DDX(J * g11) + coordinates.DDY(tmp) + coordinates.DDZ(J * g13)) / J);
  tmp = J * g22;
  coordinates.communicate(tmp);
  setG2((coordinates.DDX(J * g12) + coordinates.DDY(tmp) + coordinates.DDZ(J * g23)) / J);
  tmp = J * g23;
  coordinates.communicate(tmp);
  setG3((coordinates.DDX(J * g13) + coordinates.DDY(tmp) + coordinates.DDZ(J * g33)) / J);
}
