
#include "bout/christoffel_symbols.hxx"
#include "bout/coordinates.hxx"
#include <utility>

ChristoffelSymbols::ChristoffelSymbols(
    FieldMetric G1_11, FieldMetric G1_22, FieldMetric G1_33, FieldMetric G1_12,
    FieldMetric G1_13, FieldMetric G1_23, FieldMetric G2_11, FieldMetric G2_22,
    FieldMetric G2_33, FieldMetric G2_12, FieldMetric G2_13, FieldMetric G2_23,
    FieldMetric G3_11, FieldMetric G3_22, FieldMetric G3_33, FieldMetric G3_12,
    FieldMetric G3_13, FieldMetric G3_23, DifferentialOperators* differential_operators)
    : G1_11_(std::move(G1_11)), G1_22_(std::move(G1_22)), G1_33_(std::move(G1_33)),
      G1_12_(std::move(G1_12)), G1_13_(std::move(G1_13)), G1_23_(std::move(G1_23)),
      G2_11_(std::move(G2_11)), G2_22_(std::move(G2_22)), G2_33_(std::move(G2_33)),
      G2_12_(std::move(G2_12)), G2_13_(std::move(G2_13)), G2_23_(std::move(G2_23)),
      G3_11_(std::move(G3_11)), G3_22_(std::move(G3_22)), G3_33_(std::move(G3_33)),
      G3_12_(std::move(G3_12)), G3_13_(std::move(G3_13)), G3_23_(std::move(G3_23)),
      differential_operators(differential_operators){
          ASSERT0(differential_operators != nullptr)};

ChristoffelSymbols::ChristoffelSymbols(const Coordinates& coordinates,
                                       DifferentialOperators* differential_operators)
    : differential_operators(differential_operators) {
  ASSERT0(differential_operators != nullptr);
  // Calculate Christoffel symbol terms (18 independent values)
  // Note: This calculation is completely general: metric
  // tensor can be 2D or 3D. For 2D, all DDZ terms are zero

  const auto& dx = coordinates.dx();
  const auto& dy = coordinates.dy();
  const auto& contravariantMetricTensor = coordinates.getContravariantMetricTensor();
  const auto& covariantMetricTensor = coordinates.getCovariantMetricTensor();

  const auto& g11 = contravariantMetricTensor.g11();
  const auto& g22 = contravariantMetricTensor.g22();
  const auto& g33 = contravariantMetricTensor.g33();
  const auto& g12 = contravariantMetricTensor.g12();
  const auto& g13 = contravariantMetricTensor.g13();
  const auto& g23 = contravariantMetricTensor.g23();

  const auto& g_11 = covariantMetricTensor.g11();
  const auto& g_22 = covariantMetricTensor.g22();
  const auto& g_33 = covariantMetricTensor.g33();
  const auto& g_12 = covariantMetricTensor.g12();
  const auto& g_13 = covariantMetricTensor.g13();
  const auto& g_23 = covariantMetricTensor.g23();

  G1_11_ = 0.5 * g11 * differential_operators->DDX(g_11, dx)
           + g12
                 * (differential_operators->DDX(g_12, dx)
                    - 0.5 * differential_operators->DDY(g_11, dy))
           + g13
                 * (differential_operators->DDX(g_13, dx)
                    - 0.5 * differential_operators->DDZ(g_11));
  G1_22_ = g11
               * (differential_operators->DDY(g_12, dy)
                  - 0.5 * differential_operators->DDX(g_22, dx))
           + 0.5 * g12 * differential_operators->DDY(g_22, dy)
           + g13
                 * (differential_operators->DDY(g_23, dy)
                    - 0.5 * differential_operators->DDZ(g_22));
  G1_33_ = g11
               * (differential_operators->DDZ(g_13)
                  - 0.5 * differential_operators->DDX(g_33, dx))
           + g12
                 * (differential_operators->DDZ(g_23)
                    - 0.5 * differential_operators->DDY(g_33, dy))
           + 0.5 * g13 * differential_operators->DDZ(g_33);
  G1_12_ = 0.5 * g11 * differential_operators->DDY(g_11, dy)
           + 0.5 * g12 * differential_operators->DDX(g_22, dx)
           + 0.5 * g13
                 * (differential_operators->DDY(g_13, dy)
                    + differential_operators->DDX(g_23, dx)
                    - differential_operators->DDZ(g_12));
  G1_13_ =
      0.5 * g11 * differential_operators->DDZ(g_11)
      + 0.5 * g12
            * (differential_operators->DDZ(g_12) + differential_operators->DDX(g_23, dx)
               - differential_operators->DDY(g_13, dy))
      + 0.5 * g13 * differential_operators->DDX(g_33, dx);
  G1_23_ =
      0.5 * g11
          * (differential_operators->DDZ(g_12) + differential_operators->DDY(g_13, dy)
             - differential_operators->DDX(g_23, dx))
      + 0.5 * g12
            * (differential_operators->DDZ(g_22) + differential_operators->DDY(g_23, dy)
               - differential_operators->DDY(g_23, dy))
      // + 0.5 *g13*(differential_operators->DDZ(g_32) + differential_operators->DDY(g_33) - differential_operators->DDZ(g_23));
      // which equals
      + 0.5 * g13 * differential_operators->DDY(g_33, dy);

  G2_11_ = 0.5 * g12 * differential_operators->DDX(g_11, dx)
           + g22
                 * (differential_operators->DDX(g_12, dx)
                    - 0.5 * differential_operators->DDY(g_11, dy))
           + g23
                 * (differential_operators->DDX(g_13, dx)
                    - 0.5 * differential_operators->DDZ(g_11));
  G2_22_ = g12
               * (differential_operators->DDY(g_12, dy)
                  - 0.5 * differential_operators->DDX(g_22, dx))
           + 0.5 * g22 * differential_operators->DDY(g_22, dy)
           + g23
                 * (differential_operators->DDY(g23, dy)
                    - 0.5 * differential_operators->DDZ(g_22));
  G2_33_ = g12
               * (differential_operators->DDZ(g_13)
                  - 0.5 * differential_operators->DDX(g_33, dx))
           + g22
                 * (differential_operators->DDZ(g_23)
                    - 0.5 * differential_operators->DDY(g_33, dy))
           + 0.5 * g23 * differential_operators->DDZ(g_33);
  G2_12_ = 0.5 * g12 * differential_operators->DDY(g_11, dy)
           + 0.5 * g22 * differential_operators->DDX(g_22, dx)
           + 0.5 * g23
                 * (differential_operators->DDY(g_13, dy)
                    + differential_operators->DDX(g_23, dx)
                    - differential_operators->DDZ(g_12));
  G2_13_ =
      // 0.5 *g21*(differential_operators->DDZ(g_11) + differential_operators->DDX(covariantMetricTensor.Getg13()) - differential_operators->DDX(g_13))
      // which equals
      0.5 * g12
          * (differential_operators->DDZ(g_11) + differential_operators->DDX(g_13, dx)
             - differential_operators->DDX(g_13, dx))
      // + 0.5 *g22*(differential_operators->DDZ(covariantMetricTensor.Getg21()) + differential_operators->DDX(g_23) - differential_operators->DDY(g_13))
      // which equals
      + 0.5 * g22
            * (differential_operators->DDZ(g_12) + differential_operators->DDX(g_23, dx)
               - differential_operators->DDY(g_13, dy))
      // + 0.5 *g23*(differential_operators->DDZ(covariantMetricTensor.Getg31()) + differential_operators->DDX(g_33) - differential_operators->DDZ(g_13));
      // which equals
      + 0.5 * g23 * differential_operators->DDX(g_33, dx);
  G2_23_ =
      0.5 * g12
          * (differential_operators->DDZ(g_12) + differential_operators->DDY(g_13, dy)
             - differential_operators->DDX(g_23, dx))
      + 0.5 * g22 * differential_operators->DDZ(g_22)
      + 0.5 * g23 * differential_operators->DDY(g_33, dy);

  G3_11_ = 0.5 * g13 * differential_operators->DDX(g_11, dx)
           + g23
                 * (differential_operators->DDX(g_12, dx)
                    - 0.5 * differential_operators->DDY(g_11, dy))
           + g33
                 * (differential_operators->DDX(g_13, dx)
                    - 0.5 * differential_operators->DDZ(g_11));
  G3_22_ = g13
               * (differential_operators->DDY(g_12, dy)
                  - 0.5 * differential_operators->DDX(g_22, dx))
           + 0.5 * g23 * differential_operators->DDY(g_22, dy)
           + g33
                 * (differential_operators->DDY(g_23, dy)
                    - 0.5 * differential_operators->DDZ(g_22));
  G3_33_ = g13
               * (differential_operators->DDZ(g_13)
                  - 0.5 * differential_operators->DDX(g_33, dx))
           + g23
                 * (differential_operators->DDZ(g_23)
                    - 0.5 * differential_operators->DDY(g_33, dy))
           + 0.5 * g33 * differential_operators->DDZ(g_33);
  G3_12_ =
      // 0.5 *g31*(differential_operators->DDY(g_11) + differential_operators->DDX(covariantMetricTensor.Getg12()) - differential_operators->DDX(g_12))
      // which equals to
      0.5 * g13 * differential_operators->DDY(g_11, dy)
      // + 0.5 *g32*(differential_operators->DDY(covariantMetricTensor.Getg21()) + differential_operators->DDX(g_22) - differential_operators->DDY(g_12))
      // which equals to
      + 0.5 * g23 * differential_operators->DDX(g_22, dx)
      //+ 0.5 *g33*(differential_operators->DDY(covariantMetricTensor.Getg31()) + differential_operators->DDX(covariantMetricTensor.Getg32()) - differential_operators->DDZ(g_12));
      // which equals to
      + 0.5 * g33 * (differential_operators->DDY(g_13, dy))
      + differential_operators->DDX(g_23, dx) - differential_operators->DDZ(g_12);
  G3_13_ =
      0.5 * g13 * differential_operators->DDZ(g_11)
      + 0.5 * g23
            * (differential_operators->DDZ(g_12) + differential_operators->DDX(g_23, dx)
               - differential_operators->DDY(g_13, dy))
      + 0.5 * g33 * differential_operators->DDX(g_33, dx);
  G3_23_ =
      0.5 * g13
          * (differential_operators->DDZ(g_12) + differential_operators->DDY(g_13, dy))
      - differential_operators->DDX(g_23, dx)
      + 0.5 * g23 * differential_operators->DDZ(g_22)
      + 0.5 * g33 * differential_operators->DDY(g_33, dy);
}

ChristoffelSymbols::ChristoffelSymbols(DifferentialOperators* differential_operators)
    : differential_operators(differential_operators) {
  ASSERT0(differential_operators != nullptr);
}

const FieldMetric& ChristoffelSymbols::G1_11() const { return G1_11_; }
const FieldMetric& ChristoffelSymbols::G1_22() const { return G1_22_; }
const FieldMetric& ChristoffelSymbols::G1_33() const { return G1_33_; }
const FieldMetric& ChristoffelSymbols::G1_12() const { return G1_12_; }
const FieldMetric& ChristoffelSymbols::G1_13() const { return G1_13_; }
const FieldMetric& ChristoffelSymbols::G1_23() const { return G1_23_; }

const FieldMetric& ChristoffelSymbols::G2_11() const { return G2_11_; }
const FieldMetric& ChristoffelSymbols::G2_22() const { return G2_22_; }
const FieldMetric& ChristoffelSymbols::G2_33() const { return G2_33_; }
const FieldMetric& ChristoffelSymbols::G2_12() const { return G2_12_; }
const FieldMetric& ChristoffelSymbols::G2_13() const { return G2_13_; }
const FieldMetric& ChristoffelSymbols::G2_23() const { return G2_23_; }

const FieldMetric& ChristoffelSymbols::G3_11() const { return G3_11_; }
const FieldMetric& ChristoffelSymbols::G3_22() const { return G3_22_; }
const FieldMetric& ChristoffelSymbols::G3_33() const { return G3_33_; }
const FieldMetric& ChristoffelSymbols::G3_12() const { return G3_12_; }
const FieldMetric& ChristoffelSymbols::G3_13() const { return G3_13_; }
const FieldMetric& ChristoffelSymbols::G3_23() const { return G3_23_; }

void ChristoffelSymbols::setChristoffelSymbols(
    const FieldMetric& G1_11, const FieldMetric& G1_22, const FieldMetric& G1_33,
    const FieldMetric& G1_12, const FieldMetric& G1_13, const FieldMetric& G1_23,
    const FieldMetric& G2_11, const FieldMetric& G2_22, const FieldMetric& G2_33,
    const FieldMetric& G2_12, const FieldMetric& G2_13, const FieldMetric& G2_23,
    const FieldMetric& G3_11, const FieldMetric& G3_22, const FieldMetric& G3_33,
    const FieldMetric& G3_12, const FieldMetric& G3_13, const FieldMetric& G3_23) {
  G1_11_ = G1_11;
  G1_22_ = G1_22;
  G1_33_ = G1_33;
  G1_12_ = G1_12;
  G1_13_ = G1_13;
  G1_23_ = G1_23;
  G2_11_ = G2_11;
  G2_22_ = G2_22;
  G2_33_ = G2_33;
  G2_12_ = G2_12;
  G2_13_ = G2_13;
  G2_23_ = G2_23;
  G3_11_ = G3_11;
  G3_22_ = G3_22;
  G3_33_ = G3_33;
  G3_12_ = G3_12;
  G3_13_ = G3_13;
  G3_23_ = G3_23;
}

void ChristoffelSymbols::map(
    const std::function<const FieldMetric(const FieldMetric)>& function) {

  const ChristoffelSymbols updated_christoffel_symbols = applyToComponents(function);

  setChristoffelSymbols(
      updated_christoffel_symbols.G1_11_, updated_christoffel_symbols.G1_22_,
      updated_christoffel_symbols.G1_33_, updated_christoffel_symbols.G1_12_,
      updated_christoffel_symbols.G1_13_, updated_christoffel_symbols.G1_23_,
      updated_christoffel_symbols.G2_11_, updated_christoffel_symbols.G2_22_,
      updated_christoffel_symbols.G2_33_, updated_christoffel_symbols.G2_12_,
      updated_christoffel_symbols.G2_13_, updated_christoffel_symbols.G2_23_,
      updated_christoffel_symbols.G3_11_, updated_christoffel_symbols.G3_22_,
      updated_christoffel_symbols.G3_33_, updated_christoffel_symbols.G3_12_,
      updated_christoffel_symbols.G3_13_, updated_christoffel_symbols.G3_23_);
}

ChristoffelSymbols ChristoffelSymbols::applyToComponents(
    const std::function<const FieldMetric(const FieldMetric)>& function) const {

  const auto components_in = std::vector<FieldMetric>{
      G1_11_, G1_22_, G1_33_, G1_12_, G1_13_, G1_23_, G2_11_, G2_22_, G2_33_,
      G2_12_, G2_13_, G2_23_, G3_11_, G3_22_, G3_33_, G3_12_, G3_13_, G3_23_};

  FieldMetric components_out[18];

  std::transform(components_in.begin(), components_in.end(), components_out, function);
  auto [G1_11, G1_22, G1_33, G1_12, G1_13, G1_23, G2_11, G2_22, G2_33, G2_12, G2_13,
        G2_23, G3_11, G3_22, G3_33, G3_12, G3_13, G3_23] = components_out;

  return ChristoffelSymbols(G1_11, G1_22, G1_33, G1_12, G1_13, G1_23, G2_11, G2_22, G2_33,
                            G2_12, G2_13, G2_23, G3_11, G3_22, G3_33, G3_12, G3_13, G3_23,
                            differential_operators);
}
