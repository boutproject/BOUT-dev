/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

#include "bout/geometry.hxx"

Geometry::Geometry(const FieldMetric& J, const FieldMetric& Bxy, const FieldMetric& g11,
                   const FieldMetric& g22, const FieldMetric& g33, const FieldMetric& g12,
                   const FieldMetric& g13, const FieldMetric& g23,
                   const FieldMetric& g_11, const FieldMetric& g_22,
                   const FieldMetric& g_33, const FieldMetric& g_12,
                   const FieldMetric& g_13, const FieldMetric& g_23,
                   DifferentialOperators* differential_operators, Mesh* mesh)
    : christoffel_symbols(differential_operators),
      contravariantMetricTensor(g11, g22, g33, g12, g13, g23),
      covariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23), this_J(J), this_Bxy(Bxy),
      differential_operators(differential_operators) {
  ASSERT0(differential_operators != nullptr);
}

Geometry::Geometry(Mesh* mesh, DifferentialOperators* differential_operators)
    //bool force_interpolate_from_centre
    : christoffel_symbols(mesh, differential_operators),
      contravariantMetricTensor(1., 1., 1., 0, 0, 0, mesh),
      // Identity metric tensor
      covariantMetricTensor(1., 1., 1., 0, 0, 0, mesh), this_J(1., mesh),
      this_Bxy(1., mesh), differential_operators(differential_operators) {
  ASSERT0(differential_operators != nullptr);
}

void Geometry::CalculateChristoffelSymbols(const FieldMetric& dx, const FieldMetric& dy) {
  christoffel_symbols.CalculateChristoffelSymbols(dx, dy, contravariantMetricTensor,
                                                  covariantMetricTensor);
}

MetricTensor::FieldMetric Geometry::recalculateJacobian() {

  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)
  auto g = contravariantMetricTensor.Getg11() * contravariantMetricTensor.Getg22()
               * contravariantMetricTensor.Getg33()
           + 2.0 * contravariantMetricTensor.Getg12() * contravariantMetricTensor.Getg13()
                 * contravariantMetricTensor.Getg23()
           - contravariantMetricTensor.Getg11() * contravariantMetricTensor.Getg23()
                 * contravariantMetricTensor.Getg23()
           - contravariantMetricTensor.Getg22() * contravariantMetricTensor.Getg13()
                 * contravariantMetricTensor.Getg13()
           - contravariantMetricTensor.Getg33() * contravariantMetricTensor.Getg12()
                 * contravariantMetricTensor.Getg12();

  // Check that g is positive
  bout::checkPositive(g, "The determinant of g^ij", "RGN_NOBNDRY");

  return 1. / sqrt(g);
}

MetricTensor::FieldMetric Geometry::recalculateBxy() {

  return sqrt(covariantMetricTensor.Getg22()) / this_J;
}

void Geometry::checkCovariant(int ystart) { covariantMetricTensor.check(ystart); }

void Geometry::checkContravariant(int ystart) { contravariantMetricTensor.check(ystart); }

void Geometry::setContravariantMetricTensor(MetricTensor metric_tensor,
                                            const std::string& region) {
  contravariantMetricTensor.setMetricTensor(metric_tensor);
  covariantMetricTensor.setMetricTensor(contravariantMetricTensor.inverse(region));
}

void Geometry::setCovariantMetricTensor(MetricTensor metric_tensor,
                                        const std::string& region) {
  covariantMetricTensor.setMetricTensor(metric_tensor);
  contravariantMetricTensor.setMetricTensor(covariantMetricTensor.inverse(region));
}

const MetricTensor::FieldMetric& Geometry::g_11() const {
  return covariantMetricTensor.Getg11();
}
const MetricTensor::FieldMetric& Geometry::g_22() const {
  return covariantMetricTensor.Getg22();
}
const MetricTensor::FieldMetric& Geometry::g_33() const {
  return covariantMetricTensor.Getg33();
}
const MetricTensor::FieldMetric& Geometry::g_12() const {
  return covariantMetricTensor.Getg12();
}
const MetricTensor::FieldMetric& Geometry::g_13() const {
  return covariantMetricTensor.Getg13();
}
const MetricTensor::FieldMetric& Geometry::g_23() const {
  return covariantMetricTensor.Getg23();
}

const MetricTensor::FieldMetric& Geometry::g11() const {
  return contravariantMetricTensor.Getg11();
}
const MetricTensor::FieldMetric& Geometry::g22() const {
  return contravariantMetricTensor.Getg22();
}
const MetricTensor::FieldMetric& Geometry::g33() const {
  return contravariantMetricTensor.Getg33();
}
const MetricTensor::FieldMetric& Geometry::g12() const {
  return contravariantMetricTensor.Getg12();
}
const MetricTensor::FieldMetric& Geometry::g13() const {
  return contravariantMetricTensor.Getg13();
}
const MetricTensor::FieldMetric& Geometry::g23() const {
  return contravariantMetricTensor.Getg23();
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
  this_J = J;
}

void Geometry::setJ(BoutReal value, int x, int y) {
  //TODO: Calculate Bxy and check value is close
  this_J(x, y) = value;
}

void Geometry::setBxy(FieldMetric Bxy) {
  //TODO: Calculate Bxy and check value is close
  this_Bxy = Bxy;
}

const MetricTensor& Geometry::getContravariantMetricTensor() const {
  return contravariantMetricTensor;
}

const MetricTensor& Geometry::getCovariantMetricTensor() const {
  return covariantMetricTensor;
}

void Geometry::applyToContravariantMetricTensor(
    std::function<const FieldMetric(const FieldMetric)> function) {
  contravariantMetricTensor.map(std::move(function));
}

void Geometry::applyToCovariantMetricTensor(
    std::function<const FieldMetric(const FieldMetric)> function) {
  covariantMetricTensor.map(std::move(function));
}

void Geometry::applyToChristoffelSymbols(
    std::function<const FieldMetric(const FieldMetric)> function) {
  christoffel_symbols.map(std::move(function));
}
