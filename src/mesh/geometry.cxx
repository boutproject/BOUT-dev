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
                   DifferentialOperators* differential_operators)
    : contravariantMetricTensor(g11, g22, g33, g12, g13, g23),
      covariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23), this_J(J), this_Bxy(Bxy),
      differential_operators(differential_operators) {
  ASSERT0(differential_operators != nullptr);
}

Geometry::Geometry(Mesh* mesh, DifferentialOperators* differential_operators)
    //bool force_interpolate_from_centre
    : G1_11_(mesh), G1_22_(mesh), G1_33_(mesh), G1_12_(mesh), G1_13_(mesh), G1_23_(mesh),
      G2_11_(mesh), G2_22_(mesh), G2_33_(mesh), G2_12_(mesh), G2_13_(mesh), G2_23_(mesh),
      G3_11_(mesh), G3_22_(mesh), G3_33_(mesh), G3_12_(mesh), G3_13_(mesh), G3_23_(mesh),
      G1_(mesh), G2_(mesh), G3_(mesh),
      contravariantMetricTensor(1., 1., 1., 0, 0, 0, mesh),
      covariantMetricTensor(1., 1., 1., 0, 0, 0, mesh),
      // Identity metric tensor
      this_J(1., mesh), this_Bxy(1., mesh),
      differential_operators(differential_operators) {
  ASSERT0(differential_operators != nullptr);
}

void Geometry::CalculateChristoffelSymbols(FieldMetric& dx, FieldMetric& dy) {
  // Calculate Christoffel symbol terms (18 independent values)
  // Note: This calculation is completely general: metric
  // tensor can be 2D or 3D. For 2D, all DDZ terms are zero

  G1_11_ =
      0.5 * contravariantMetricTensor.Getg11()
          * differential_operators->DDX(covariantMetricTensor.Getg11(), dx)
      + contravariantMetricTensor.Getg12()
            * (differential_operators->DDX(covariantMetricTensor.Getg12(), dx)
               - 0.5 * differential_operators->DDY(covariantMetricTensor.Getg11(), dy))
      + contravariantMetricTensor.Getg13()
            * (differential_operators->DDX(covariantMetricTensor.Getg13(), dx)
               - 0.5 * differential_operators->DDZ(covariantMetricTensor.Getg11()));
  G1_22_ = contravariantMetricTensor.Getg11()
               * (differential_operators->DDY(covariantMetricTensor.Getg12(), dy)
                  - 0.5 * differential_operators->DDX(covariantMetricTensor.Getg22(), dx))
           + 0.5 * contravariantMetricTensor.Getg12()
                 * differential_operators->DDY(covariantMetricTensor.Getg22(), dy)
           + contravariantMetricTensor.Getg13()
                 * (differential_operators->DDY(covariantMetricTensor.Getg23(), dy)
                    - 0.5 * differential_operators->DDZ(covariantMetricTensor.Getg22()));
  G1_33_ =
      contravariantMetricTensor.Getg11()
          * (differential_operators->DDZ(covariantMetricTensor.Getg13())
             - 0.5 * differential_operators->DDX(covariantMetricTensor.Getg33(), dx))
      + contravariantMetricTensor.Getg12()
            * (differential_operators->DDZ(covariantMetricTensor.Getg23())
               - 0.5 * differential_operators->DDY(covariantMetricTensor.Getg33(), dy))
      + 0.5 * contravariantMetricTensor.Getg13()
            * differential_operators->DDZ(covariantMetricTensor.Getg33());
  G1_12_ = 0.5 * contravariantMetricTensor.Getg11()
               * differential_operators->DDY(covariantMetricTensor.Getg11(), dy)
           + 0.5 * contravariantMetricTensor.Getg12()
                 * differential_operators->DDX(covariantMetricTensor.Getg22(), dx)
           + 0.5 * contravariantMetricTensor.Getg13()
                 * (differential_operators->DDY(covariantMetricTensor.Getg13(), dy)
                    + differential_operators->DDX(covariantMetricTensor.Getg23(), dx)
                    - differential_operators->DDZ(covariantMetricTensor.Getg12()));
  G1_13_ = 0.5 * contravariantMetricTensor.Getg11()
               * differential_operators->DDZ(covariantMetricTensor.Getg11())
           + 0.5 * contravariantMetricTensor.Getg12()
                 * (differential_operators->DDZ(covariantMetricTensor.Getg12())
                    + differential_operators->DDX(covariantMetricTensor.Getg23(), dx)
                    - differential_operators->DDY(covariantMetricTensor.Getg13(), dy))
           + 0.5 * contravariantMetricTensor.Getg13()
                 * differential_operators->DDX(covariantMetricTensor.Getg33(), dx);
  G1_23_ =
      0.5 * contravariantMetricTensor.Getg11()
          * (differential_operators->DDZ(covariantMetricTensor.Getg12())
             + differential_operators->DDY(covariantMetricTensor.Getg13(), dy)
             - differential_operators->DDX(covariantMetricTensor.Getg23(), dx))
      + 0.5 * contravariantMetricTensor.Getg12()
            * (differential_operators->DDZ(covariantMetricTensor.Getg22())
               + differential_operators->DDY(covariantMetricTensor.Getg23(), dy)
               - differential_operators->DDY(covariantMetricTensor.Getg23(), dy))
      // + 0.5 *g13*(differential_operators->DDZ(g_32) + differential_operators->DDY(g_33) - differential_operators->DDZ(g_23));
      // which equals
      + 0.5 * contravariantMetricTensor.Getg13()
            * differential_operators->DDY(covariantMetricTensor.Getg33(), dy);

  G2_11_ =
      0.5 * contravariantMetricTensor.Getg12()
          * differential_operators->DDX(covariantMetricTensor.Getg11(), dx)
      + contravariantMetricTensor.Getg22()
            * (differential_operators->DDX(covariantMetricTensor.Getg12(), dx)
               - 0.5 * differential_operators->DDY(covariantMetricTensor.Getg11(), dy))
      + contravariantMetricTensor.Getg23()
            * (differential_operators->DDX(covariantMetricTensor.Getg13(), dx)
               - 0.5 * differential_operators->DDZ(covariantMetricTensor.Getg11()));
  G2_22_ = contravariantMetricTensor.Getg12()
               * (differential_operators->DDY(covariantMetricTensor.Getg12(), dy)
                  - 0.5 * differential_operators->DDX(covariantMetricTensor.Getg22(), dx))
           + 0.5 * contravariantMetricTensor.Getg22()
                 * differential_operators->DDY(covariantMetricTensor.Getg22(), dy)
           + contravariantMetricTensor.Getg23()
                 * (differential_operators->DDY(contravariantMetricTensor.Getg23(), dy)
                    - 0.5 * differential_operators->DDZ(covariantMetricTensor.Getg22()));
  G2_33_ =
      contravariantMetricTensor.Getg12()
          * (differential_operators->DDZ(covariantMetricTensor.Getg13())
             - 0.5 * differential_operators->DDX(covariantMetricTensor.Getg33(), dx))
      + contravariantMetricTensor.Getg22()
            * (differential_operators->DDZ(covariantMetricTensor.Getg23())
               - 0.5 * differential_operators->DDY(covariantMetricTensor.Getg33(), dy))
      + 0.5 * contravariantMetricTensor.Getg23()
            * differential_operators->DDZ(covariantMetricTensor.Getg33());
  G2_12_ = 0.5 * contravariantMetricTensor.Getg12()
               * differential_operators->DDY(covariantMetricTensor.Getg11(), dy)
           + 0.5 * contravariantMetricTensor.Getg22()
                 * differential_operators->DDX(covariantMetricTensor.Getg22(), dx)
           + 0.5 * contravariantMetricTensor.Getg23()
                 * (differential_operators->DDY(covariantMetricTensor.Getg13(), dy)
                    + differential_operators->DDX(covariantMetricTensor.Getg23(), dx)
                    - differential_operators->DDZ(covariantMetricTensor.Getg12()));
  G2_13_ =
      // 0.5 *g21*(differential_operators->DDZ(covariantMetricTensor.Getg11()) + differential_operators->DDX(covariantMetricTensor.Getg13()) - differential_operators->DDX(covariantMetricTensor.Getg13()))
      // which equals
      0.5 * contravariantMetricTensor.Getg12()
          * (differential_operators->DDZ(covariantMetricTensor.Getg11())
             + differential_operators->DDX(covariantMetricTensor.Getg13(), dx)
             - differential_operators->DDX(covariantMetricTensor.Getg13(), dx))
      // + 0.5 *g22*(differential_operators->DDZ(covariantMetricTensor.Getg21()) + differential_operators->DDX(covariantMetricTensor.Getg23()) - differential_operators->DDY(covariantMetricTensor.Getg13()))
      // which equals
      + 0.5 * contravariantMetricTensor.Getg22()
            * (differential_operators->DDZ(covariantMetricTensor.Getg12())
               + differential_operators->DDX(covariantMetricTensor.Getg23(), dx)
               - differential_operators->DDY(covariantMetricTensor.Getg13(), dy))
      // + 0.5 *g23*(differential_operators->DDZ(covariantMetricTensor.Getg31()) + differential_operators->DDX(covariantMetricTensor.Getg33()) - differential_operators->DDZ(g_13));
      // which equals
      + 0.5 * contravariantMetricTensor.Getg23()
            * differential_operators->DDX(covariantMetricTensor.Getg33(), dx);
  G2_23_ = 0.5 * contravariantMetricTensor.Getg12()
               * (differential_operators->DDZ(covariantMetricTensor.Getg12())
                  + differential_operators->DDY(covariantMetricTensor.Getg13(), dy)
                  - differential_operators->DDX(covariantMetricTensor.Getg23(), dx))
           + 0.5 * contravariantMetricTensor.Getg22()
                 * differential_operators->DDZ(covariantMetricTensor.Getg22())
           + 0.5 * contravariantMetricTensor.Getg23()
                 * differential_operators->DDY(covariantMetricTensor.Getg33(), dy);

  G3_11_ =
      0.5 * contravariantMetricTensor.Getg13()
          * differential_operators->DDX(covariantMetricTensor.Getg11(), dx)
      + contravariantMetricTensor.Getg23()
            * (differential_operators->DDX(covariantMetricTensor.Getg12(), dx)
               - 0.5 * differential_operators->DDY(covariantMetricTensor.Getg11(), dy))
      + contravariantMetricTensor.Getg33()
            * (differential_operators->DDX(covariantMetricTensor.Getg13(), dx)
               - 0.5 * differential_operators->DDZ(covariantMetricTensor.Getg11()));
  G3_22_ = contravariantMetricTensor.Getg13()
               * (differential_operators->DDY(covariantMetricTensor.Getg12(), dy)
                  - 0.5 * differential_operators->DDX(covariantMetricTensor.Getg22(), dx))
           + 0.5 * contravariantMetricTensor.Getg23()
                 * differential_operators->DDY(covariantMetricTensor.Getg22(), dy)
           + contravariantMetricTensor.Getg33()
                 * (differential_operators->DDY(covariantMetricTensor.Getg23(), dy)
                    - 0.5 * differential_operators->DDZ(covariantMetricTensor.Getg22()));
  G3_33_ =
      contravariantMetricTensor.Getg13()
          * (differential_operators->DDZ(covariantMetricTensor.Getg13())
             - 0.5 * differential_operators->DDX(covariantMetricTensor.Getg33(), dx))
      + contravariantMetricTensor.Getg23()
            * (differential_operators->DDZ(covariantMetricTensor.Getg23())
               - 0.5 * differential_operators->DDY(covariantMetricTensor.Getg33(), dy))
      + 0.5 * contravariantMetricTensor.Getg33()
            * differential_operators->DDZ(covariantMetricTensor.Getg33());
  G3_12_ =
      // 0.5 *g31*(differential_operators->DDY(covariantMetricTensor.Getg11()) + differential_operators->DDX(covariantMetricTensor.Getg12()) - differential_operators->DDX(covariantMetricTensor.Getg12()))
      // which equals to
      0.5 * contravariantMetricTensor.Getg13()
          * differential_operators->DDY(covariantMetricTensor.Getg11(), dy)
      // + 0.5 *g32*(differential_operators->DDY(covariantMetricTensor.Getg21()) + differential_operators->DDX(covariantMetricTensor.Getg22()) - differential_operators->DDY(covariantMetricTensor.Getg12()))
      // which equals to
      + 0.5 * contravariantMetricTensor.Getg23()
            * differential_operators->DDX(covariantMetricTensor.Getg22(), dx)
      //+ 0.5 *g33*(differential_operators->DDY(covariantMetricTensor.Getg31()) + differential_operators->DDX(covariantMetricTensor.Getg32()) - differential_operators->DDZ(covariantMetricTensor.Getg12()));
      // which equals to
      + 0.5 * contravariantMetricTensor.Getg33()
            * (differential_operators->DDY(covariantMetricTensor.Getg13(), dy))
      + differential_operators->DDX(covariantMetricTensor.Getg23(), dx)
      - differential_operators->DDZ(covariantMetricTensor.Getg12());
  G3_13_ = 0.5 * contravariantMetricTensor.Getg13()
               * differential_operators->DDZ(covariantMetricTensor.Getg11())
           + 0.5 * contravariantMetricTensor.Getg23()
                 * (differential_operators->DDZ(covariantMetricTensor.Getg12())
                    + differential_operators->DDX(covariantMetricTensor.Getg23(), dx)
                    - differential_operators->DDY(covariantMetricTensor.Getg13(), dy))
           + 0.5 * contravariantMetricTensor.Getg33()
                 * differential_operators->DDX(covariantMetricTensor.Getg33(), dx);
  G3_23_ = 0.5 * contravariantMetricTensor.Getg13()
               * (differential_operators->DDZ(covariantMetricTensor.Getg12())
                  + differential_operators->DDY(covariantMetricTensor.Getg13(), dy))
           - differential_operators->DDX(covariantMetricTensor.Getg23(), dx)
           + 0.5 * contravariantMetricTensor.Getg23()
                 * differential_operators->DDZ(covariantMetricTensor.Getg22())
           + 0.5 * contravariantMetricTensor.Getg33()
                 * differential_operators->DDY(covariantMetricTensor.Getg33(), dy);
}

void Geometry::calcCovariant(const std::string& region) {
  TRACE("Geometry::calcCovariant");
  covariantMetricTensor.setMetricTensor(contravariantMetricTensor.inverse(region));
}

void Geometry::calcContravariant(const std::string& region) {
  TRACE("Geometry::calcContravariant");
  contravariantMetricTensor.setMetricTensor(covariantMetricTensor.inverse(region));
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

const FieldMetric& Geometry::G1_11() const { return G1_11_; }
const FieldMetric& Geometry::G1_22() const { return G1_22_; }
const FieldMetric& Geometry::G1_33() const { return G1_33_; }
const FieldMetric& Geometry::G1_12() const { return G1_12_; }
const FieldMetric& Geometry::G1_13() const { return G1_13_; }
const FieldMetric& Geometry::G1_23() const { return G1_23_; }

const FieldMetric& Geometry::G2_11() const { return G2_11_; }
const FieldMetric& Geometry::G2_22() const { return G2_22_; }
const FieldMetric& Geometry::G2_33() const { return G2_33_; }
const FieldMetric& Geometry::G2_12() const { return G2_12_; }
const FieldMetric& Geometry::G2_13() const { return G2_13_; }
const FieldMetric& Geometry::G2_23() const { return G2_23_; }

const FieldMetric& Geometry::G3_11() const { return G3_11_; }
const FieldMetric& Geometry::G3_22() const { return G3_22_; }
const FieldMetric& Geometry::G3_33() const { return G3_33_; }
const FieldMetric& Geometry::G3_12() const { return G3_12_; }
const FieldMetric& Geometry::G3_13() const { return G3_13_; }
const FieldMetric& Geometry::G3_23() const { return G3_23_; }

const FieldMetric& Geometry::G1() const { return G1_; }
const FieldMetric& Geometry::G2() const { return G2_; }
const FieldMetric& Geometry::G3() const { return G3_; }

const FieldMetric& Geometry::J() const { return this_J; }

const FieldMetric& Geometry::Bxy() const { return this_Bxy; }

void Geometry::setG1_11(FieldMetric G1_11) { G1_11_ = G1_11; }
void Geometry::setG1_22(FieldMetric G1_22) { G1_11_ = G1_22; }
void Geometry::setG1_33(FieldMetric G1_33) { G1_11_ = G1_33; }
void Geometry::setG1_12(FieldMetric G1_12) { G1_11_ = G1_12; }
void Geometry::setG1_13(FieldMetric G1_13) { G1_11_ = G1_13; }
void Geometry::setG1_23(FieldMetric G1_23) { G1_11_ = G1_23; }

void Geometry::setG2_11(FieldMetric G2_11) { G2_11_ = G2_11; }
void Geometry::setG2_22(FieldMetric G2_22) { G2_11_ = G2_22; }
void Geometry::setG2_33(FieldMetric G2_33) { G2_11_ = G2_33; }
void Geometry::setG2_12(FieldMetric G2_12) { G2_11_ = G2_12; }
void Geometry::setG2_13(FieldMetric G2_13) { G2_11_ = G2_13; }
void Geometry::setG2_23(FieldMetric G2_23) { G2_11_ = G2_23; }

void Geometry::setG3_11(FieldMetric G3_11) { G3_11_ = G3_11; }
void Geometry::setG3_22(FieldMetric G3_22) { G3_11_ = G3_22; }
void Geometry::setG3_33(FieldMetric G3_33) { G3_11_ = G3_33; }
void Geometry::setG3_12(FieldMetric G3_12) { G3_11_ = G3_12; }
void Geometry::setG3_13(FieldMetric G3_13) { G3_11_ = G3_13; }
void Geometry::setG3_23(FieldMetric G3_23) { G3_11_ = G3_23; }

void Geometry::setG1(FieldMetric G1) { G1_ = G1; }
void Geometry::setG2(FieldMetric G2) { G2_ = G2; }
void Geometry::setG3(FieldMetric G3) { G3_ = G3; }

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
  contravariantMetricTensor.map(function);
}

void Geometry::applyToCovariantMetricTensor(
    std::function<const FieldMetric(const FieldMetric)> function) {
  covariantMetricTensor.map(std::move(function));
}

void Geometry::applyToChristoffelSymbols(
    const std::function<const FieldMetric(const FieldMetric)> function) {

  const auto components_in =
      std::vector<FieldMetric>{G1_11_, G1_22_, G1_33_, G1_12_, G1_13_, G1_23_, G2_11_,
                               G2_22_, G2_33_, G2_12_, G2_13_, G2_23_, G3_11_, G3_22_,
                               G3_33_, G3_12_, G3_13_, G3_23_, G1_,    G2_,    G3_};

  FieldMetric components_out[21] = {G1_11_, G1_22_, G1_33_, G1_12_, G1_13_, G1_23_,
                                    G2_11_, G2_22_, G2_33_, G2_12_, G2_13_, G2_23_,
                                    G3_11_, G3_22_, G3_33_, G3_12_, G3_13_, G3_23_,
                                    G1_,    G2_,    G3_};

  std::transform(components_in.begin(), components_in.end(), components_out, function);
}
