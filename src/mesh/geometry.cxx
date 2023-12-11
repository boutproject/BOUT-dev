/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

#include "bout/geometry.hxx"
#include <utility>

Geometry::Geometry(FieldMetric J, FieldMetric Bxy, const FieldMetric& g11,
                   const FieldMetric& g22, const FieldMetric& g33, const FieldMetric& g12,
                   const FieldMetric& g13, const FieldMetric& g23,
                   const FieldMetric& g_11, const FieldMetric& g_22,
                   const FieldMetric& g_33, const FieldMetric& g_12,
                   const FieldMetric& g_13, const FieldMetric& g_23,
                   DifferentialOperators* differential_operators)
    : christoffel_symbols(differential_operators),
      contravariantMetricTensor(g11, g22, g33, g12, g13, g23),
      covariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23), this_J(std::move(J)),
      this_Bxy(std::move(Bxy)), differential_operators(differential_operators){ASSERT0(
                                    differential_operators != nullptr)}

      Geometry::Geometry(Mesh * mesh, DifferentialOperators * differential_operators)
    //bool force_interpolate_from_centre
    : christoffel_symbols(mesh, differential_operators),
      contravariantMetricTensor(1., 1., 1., 0, 0, 0, mesh),
      // Identity metric tensor
      covariantMetricTensor(1., 1., 1., 0, 0, 0, mesh), this_J(1., mesh),
      this_Bxy(1., mesh), differential_operators(differential_operators) {
  ASSERT0(differential_operators != nullptr)
}

void Geometry::CalculateChristoffelSymbols(const FieldMetric& dx, const FieldMetric& dy) {
  christoffel_symbols.CalculateChristoffelSymbols(dx, dy, contravariantMetricTensor,
                                                  covariantMetricTensor);
}

MetricTensor::FieldMetric Geometry::recalculateJacobian() {

  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)
  auto g = contravariantMetricTensor.g11() * contravariantMetricTensor.g22()
               * contravariantMetricTensor.g33()
           + 2.0 * contravariantMetricTensor.g12() * contravariantMetricTensor.g13()
                 * contravariantMetricTensor.g23()
           - contravariantMetricTensor.g11() * contravariantMetricTensor.g23()
                 * contravariantMetricTensor.g23()
           - contravariantMetricTensor.g22() * contravariantMetricTensor.g13()
                 * contravariantMetricTensor.g13()
           - contravariantMetricTensor.g33() * contravariantMetricTensor.g12()
                 * contravariantMetricTensor.g12();

  // Check that g is positive
  bout::checkPositive(g, "The determinant of g^ij", "RGN_NOBNDRY");

  return 1. / sqrt(g);
}

MetricTensor::FieldMetric Geometry::recalculateBxy() {

  return sqrt(covariantMetricTensor.g22()) / this_J;
}

void Geometry::checkCovariant(int ystart) { covariantMetricTensor.check(ystart); }

void Geometry::checkContravariant(int ystart) { contravariantMetricTensor.check(ystart); }

void Geometry::setContravariantMetricTensor(const MetricTensor& metric_tensor,
                                            const std::string& region) {
  contravariantMetricTensor.setMetricTensor(metric_tensor);
  covariantMetricTensor.setMetricTensor(contravariantMetricTensor.inverse(region));
}

void Geometry::setCovariantMetricTensor(const MetricTensor& metric_tensor,
                                        const std::string& region) {
  covariantMetricTensor.setMetricTensor(metric_tensor);
  contravariantMetricTensor.setMetricTensor(covariantMetricTensor.inverse(region));
}

const MetricTensor::FieldMetric& Geometry::g_11() const {
  return covariantMetricTensor.g11();
}
const MetricTensor::FieldMetric& Geometry::g_22() const {
  return covariantMetricTensor.g22();
}
const MetricTensor::FieldMetric& Geometry::g_33() const {
  return covariantMetricTensor.g33();
}
const MetricTensor::FieldMetric& Geometry::g_12() const {
  return covariantMetricTensor.g12();
}
const MetricTensor::FieldMetric& Geometry::g_13() const {
  return covariantMetricTensor.g13();
}
const MetricTensor::FieldMetric& Geometry::g_23() const {
  return covariantMetricTensor.g23();
}

const MetricTensor::FieldMetric& Geometry::g11() const {
  return contravariantMetricTensor.g11();
}
const MetricTensor::FieldMetric& Geometry::g22() const {
  return contravariantMetricTensor.g22();
}
const MetricTensor::FieldMetric& Geometry::g33() const {
  return contravariantMetricTensor.g33();
}
const MetricTensor::FieldMetric& Geometry::g12() const {
  return contravariantMetricTensor.g12();
}
const MetricTensor::FieldMetric& Geometry::g13() const {
  return contravariantMetricTensor.g13();
}
const MetricTensor::FieldMetric& Geometry::g23() const {
  return contravariantMetricTensor.g23();
}

const FieldMetric& Geometry::G1_11() const { return christoffel_symbols.G1_11(); }
const FieldMetric& Geometry::G1_22() const { return christoffel_symbols.G1_22(); }
const FieldMetric& Geometry::G1_33() const { return christoffel_symbols.G1_33(); }
const FieldMetric& Geometry::G1_12() const { return christoffel_symbols.G1_12(); }
const FieldMetric& Geometry::G1_13() const { return christoffel_symbols.G1_13(); }
const FieldMetric& Geometry::G1_23() const { return christoffel_symbols.G1_23(); }

const FieldMetric& Geometry::G2_11() const { return christoffel_symbols.G2_11(); }
const FieldMetric& Geometry::G2_22() const { return christoffel_symbols.G2_22(); }
const FieldMetric& Geometry::G2_33() const { return christoffel_symbols.G2_33(); }
const FieldMetric& Geometry::G2_12() const { return christoffel_symbols.G2_12(); }
const FieldMetric& Geometry::G2_13() const { return christoffel_symbols.G2_13(); }
const FieldMetric& Geometry::G2_23() const { return christoffel_symbols.G2_23(); }

const FieldMetric& Geometry::G3_11() const { return christoffel_symbols.G3_11(); }
const FieldMetric& Geometry::G3_22() const { return christoffel_symbols.G3_22(); }
const FieldMetric& Geometry::G3_33() const { return christoffel_symbols.G3_33(); }
const FieldMetric& Geometry::G3_12() const { return christoffel_symbols.G3_12(); }
const FieldMetric& Geometry::G3_13() const { return christoffel_symbols.G3_13(); }
const FieldMetric& Geometry::G3_23() const { return christoffel_symbols.G3_23(); }

const FieldMetric& Geometry::G1() const { return christoffel_symbols.G1(); }
const FieldMetric& Geometry::G2() const { return christoffel_symbols.G2(); }
const FieldMetric& Geometry::G3() const { return christoffel_symbols.G3(); }

const FieldMetric& Geometry::J() const { return this_J; }

const FieldMetric& Geometry::Bxy() const { return this_Bxy; }

void Geometry::setG1(FieldMetric G1) { christoffel_symbols.setG1(G1); }
void Geometry::setG2(FieldMetric G2) { christoffel_symbols.setG2(G2); }
void Geometry::setG3(FieldMetric G3) { christoffel_symbols.setG3(G3); }

void Geometry::setJ(FieldMetric J) {
  //TODO: Calculate J and check value is close
  this_J = std::move(J);
}

void Geometry::setJ(BoutReal value, int x, int y) {
  //TODO: Calculate Bxy and check value is close
  this_J(x, y) = value;
}

void Geometry::setBxy(FieldMetric Bxy) {
  //TODO: Calculate Bxy and check value is close
  this_Bxy = std::move(Bxy);
}

const MetricTensor& Geometry::getContravariantMetricTensor() const {
  return contravariantMetricTensor;
}

const MetricTensor& Geometry::getCovariantMetricTensor() const {
  return covariantMetricTensor;
}

void Geometry::applyToContravariantMetricTensor(
    const std::function<const FieldMetric(const FieldMetric)>& function) {
  contravariantMetricTensor.map(function);
}

void Geometry::applyToCovariantMetricTensor(
    const std::function<const FieldMetric(const FieldMetric)>& function) {
  covariantMetricTensor.map(function);
}

void Geometry::applyToChristoffelSymbols(
    const std::function<const FieldMetric(const FieldMetric)>& function) {
  christoffel_symbols.map(function);
}
