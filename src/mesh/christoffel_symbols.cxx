
#include <utility>
#include "bout/christoffel_symbols.hxx"

ChristoffelSymbols::ChristoffelSymbols(
    FieldMetric G1_11, FieldMetric G1_22, FieldMetric G1_33, FieldMetric G1_12,
    FieldMetric G1_13, FieldMetric G1_23, FieldMetric G2_11, FieldMetric G2_22,
    FieldMetric G2_33, FieldMetric G2_12, FieldMetric G2_13, FieldMetric G2_23,
    FieldMetric G3_11, FieldMetric G3_22, FieldMetric G3_33, FieldMetric G3_12,
    FieldMetric G3_13, FieldMetric G3_23, FieldMetric G1, FieldMetric G2, FieldMetric G3,
    DifferentialOperators* differential_operators)
    : G1_11_(std::move(G1_11)), G1_22_(std::move(G1_22)), G1_33_(std::move(G1_33)),
      G1_12_(std::move(G1_12)), G1_13_(std::move(G1_13)), G1_23_(std::move(G1_23)),
      G2_11_(std::move(G2_11)), G2_22_(std::move(G2_22)), G2_33_(std::move(G2_33)),
      G2_12_(std::move(G2_12)), G2_13_(std::move(G2_13)), G2_23_(std::move(G2_23)),
      G3_11_(std::move(G3_11)), G3_22_(std::move(G3_22)), G3_33_(std::move(G3_33)),
      G3_12_(std::move(G3_12)), G3_13_(std::move(G3_13)), G3_23_(std::move(G3_23)),
      G1_(std::move(G1)), G2_(std::move(G2)), G3_(std::move(G3)),
      differential_operators(differential_operators){
          ASSERT0(differential_operators != nullptr)}

      /// This constructor (taking no mesh) is needed because mesh may be a FakeMesh,
      /// which has no 'source' set, leading to a SEGFAULT if the constructor leads
      /// to initialisation of FieldMetric objects.
      ChristoffelSymbols::ChristoffelSymbols(DifferentialOperators
                                             * differential_operators)
    : differential_operators(
          differential_operators){ASSERT0(differential_operators != nullptr)}

      ChristoffelSymbols::ChristoffelSymbols(Mesh * mesh, DifferentialOperators
                                                              * differential_operators)
    : G1_11_(mesh), G1_22_(mesh), G1_33_(mesh), G1_12_(mesh), G1_13_(mesh), G1_23_(mesh),
      G2_11_(mesh), G2_22_(mesh), G2_33_(mesh), G2_12_(mesh), G2_13_(mesh), G2_23_(mesh),
      G3_11_(mesh), G3_22_(mesh), G3_33_(mesh), G3_12_(mesh), G3_13_(mesh), G3_23_(mesh),
      G1_(mesh), G2_(mesh), G3_(mesh), differential_operators(differential_operators) {
  ASSERT0(differential_operators != nullptr)
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

const FieldMetric& ChristoffelSymbols::G1() const { return G1_; }
const FieldMetric& ChristoffelSymbols::G2() const { return G2_; }
const FieldMetric& ChristoffelSymbols::G3() const { return G3_; }

void ChristoffelSymbols::setG1(FieldMetric& G1) { G1_ = G1; }
void ChristoffelSymbols::setG2(FieldMetric& G2) { G2_ = G2; }
void ChristoffelSymbols::setG3(FieldMetric& G3) { G3_ = G3; }

void ChristoffelSymbols::setChristoffelSymbols(
    const FieldMetric& G1_11, const FieldMetric& G1_22, const FieldMetric& G1_33,
    const FieldMetric& G1_12, const FieldMetric& G1_13, const FieldMetric& G1_23,
    const FieldMetric& G2_11, const FieldMetric& G2_22, const FieldMetric& G2_33,
    const FieldMetric& G2_12, const FieldMetric& G2_13, const FieldMetric& G2_23,
    const FieldMetric& G3_11, const FieldMetric& G3_22, const FieldMetric& G3_33,
    const FieldMetric& G3_12, const FieldMetric& G3_13, const FieldMetric& G3_23,
    const FieldMetric& G1, const FieldMetric& G2, const FieldMetric& G3) {
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
  G1_ = G1;
  G2_ = G2;
  G3_ = G3;
}

void ChristoffelSymbols::CalculateChristoffelSymbols(
    const FieldMetric& dx, const FieldMetric& dy,
    const MetricTensor& contravariantMetricTensor,
    const MetricTensor& covariantMetricTensor) {
  // Calculate Christoffel symbol terms (18 independent values)
  // Note: This calculation is completely general: metric
  // tensor can be 2D or 3D. For 2D, all DDZ terms are zero

  G1_11_ = 0.5 * contravariantMetricTensor.g11()
               * differential_operators->DDX(covariantMetricTensor.g11(), dx)
           + contravariantMetricTensor.g12()
                 * (differential_operators->DDX(covariantMetricTensor.g12(), dx)
                    - 0.5 * differential_operators->DDY(covariantMetricTensor.g11(), dy))
           + contravariantMetricTensor.g13()
                 * (differential_operators->DDX(covariantMetricTensor.g13(), dx)
                    - 0.5 * differential_operators->DDZ(covariantMetricTensor.g11()));
  G1_22_ = contravariantMetricTensor.g11()
               * (differential_operators->DDY(covariantMetricTensor.g12(), dy)
                  - 0.5 * differential_operators->DDX(covariantMetricTensor.g22(), dx))
           + 0.5 * contravariantMetricTensor.g12()
                 * differential_operators->DDY(covariantMetricTensor.g22(), dy)
           + contravariantMetricTensor.g13()
                 * (differential_operators->DDY(covariantMetricTensor.g23(), dy)
                    - 0.5 * differential_operators->DDZ(covariantMetricTensor.g22()));
  G1_33_ = contravariantMetricTensor.g11()
               * (differential_operators->DDZ(covariantMetricTensor.g13())
                  - 0.5 * differential_operators->DDX(covariantMetricTensor.g33(), dx))
           + contravariantMetricTensor.g12()
                 * (differential_operators->DDZ(covariantMetricTensor.g23())
                    - 0.5 * differential_operators->DDY(covariantMetricTensor.g33(), dy))
           + 0.5 * contravariantMetricTensor.g13()
                 * differential_operators->DDZ(covariantMetricTensor.g33());
  G1_12_ = 0.5 * contravariantMetricTensor.g11()
               * differential_operators->DDY(covariantMetricTensor.g11(), dy)
           + 0.5 * contravariantMetricTensor.g12()
                 * differential_operators->DDX(covariantMetricTensor.g22(), dx)
           + 0.5 * contravariantMetricTensor.g13()
                 * (differential_operators->DDY(covariantMetricTensor.g13(), dy)
                    + differential_operators->DDX(covariantMetricTensor.g23(), dx)
                    - differential_operators->DDZ(covariantMetricTensor.g12()));
  G1_13_ = 0.5 * contravariantMetricTensor.g11()
               * differential_operators->DDZ(covariantMetricTensor.g11())
           + 0.5 * contravariantMetricTensor.g12()
                 * (differential_operators->DDZ(covariantMetricTensor.g12())
                    + differential_operators->DDX(covariantMetricTensor.g23(), dx)
                    - differential_operators->DDY(covariantMetricTensor.g13(), dy))
           + 0.5 * contravariantMetricTensor.g13()
                 * differential_operators->DDX(covariantMetricTensor.g33(), dx);
  G1_23_ =
      0.5 * contravariantMetricTensor.g11()
          * (differential_operators->DDZ(covariantMetricTensor.g12())
             + differential_operators->DDY(covariantMetricTensor.g13(), dy)
             - differential_operators->DDX(covariantMetricTensor.g23(), dx))
      + 0.5 * contravariantMetricTensor.g12()
            * (differential_operators->DDZ(covariantMetricTensor.g22())
               + differential_operators->DDY(covariantMetricTensor.g23(), dy)
               - differential_operators->DDY(covariantMetricTensor.g23(), dy))
      // + 0.5 *g13*(differential_operators->DDZ(g_32) + differential_operators->DDY(g_33) - differential_operators->DDZ(g_23));
      // which equals
      + 0.5 * contravariantMetricTensor.g13()
            * differential_operators->DDY(covariantMetricTensor.g33(), dy);

  G2_11_ = 0.5 * contravariantMetricTensor.g12()
               * differential_operators->DDX(covariantMetricTensor.g11(), dx)
           + contravariantMetricTensor.g22()
                 * (differential_operators->DDX(covariantMetricTensor.g12(), dx)
                    - 0.5 * differential_operators->DDY(covariantMetricTensor.g11(), dy))
           + contravariantMetricTensor.g23()
                 * (differential_operators->DDX(covariantMetricTensor.g13(), dx)
                    - 0.5 * differential_operators->DDZ(covariantMetricTensor.g11()));
  G2_22_ = contravariantMetricTensor.g12()
               * (differential_operators->DDY(covariantMetricTensor.g12(), dy)
                  - 0.5 * differential_operators->DDX(covariantMetricTensor.g22(), dx))
           + 0.5 * contravariantMetricTensor.g22()
                 * differential_operators->DDY(covariantMetricTensor.g22(), dy)
           + contravariantMetricTensor.g23()
                 * (differential_operators->DDY(contravariantMetricTensor.g23(), dy)
                    - 0.5 * differential_operators->DDZ(covariantMetricTensor.g22()));
  G2_33_ = contravariantMetricTensor.g12()
               * (differential_operators->DDZ(covariantMetricTensor.g13())
                  - 0.5 * differential_operators->DDX(covariantMetricTensor.g33(), dx))
           + contravariantMetricTensor.g22()
                 * (differential_operators->DDZ(covariantMetricTensor.g23())
                    - 0.5 * differential_operators->DDY(covariantMetricTensor.g33(), dy))
           + 0.5 * contravariantMetricTensor.g23()
                 * differential_operators->DDZ(covariantMetricTensor.g33());
  G2_12_ = 0.5 * contravariantMetricTensor.g12()
               * differential_operators->DDY(covariantMetricTensor.g11(), dy)
           + 0.5 * contravariantMetricTensor.g22()
                 * differential_operators->DDX(covariantMetricTensor.g22(), dx)
           + 0.5 * contravariantMetricTensor.g23()
                 * (differential_operators->DDY(covariantMetricTensor.g13(), dy)
                    + differential_operators->DDX(covariantMetricTensor.g23(), dx)
                    - differential_operators->DDZ(covariantMetricTensor.g12()));
  G2_13_ =
      // 0.5 *g21*(differential_operators->DDZ(covariantMetricTensor.g11()) + differential_operators->DDX(covariantMetricTensor.Getg13()) - differential_operators->DDX(covariantMetricTensor.g13()))
      // which equals
      0.5 * contravariantMetricTensor.g12()
          * (differential_operators->DDZ(covariantMetricTensor.g11())
             + differential_operators->DDX(covariantMetricTensor.g13(), dx)
             - differential_operators->DDX(covariantMetricTensor.g13(), dx))
      // + 0.5 *g22*(differential_operators->DDZ(covariantMetricTensor.Getg21()) + differential_operators->DDX(covariantMetricTensor.g23()) - differential_operators->DDY(covariantMetricTensor.g13()))
      // which equals
      + 0.5 * contravariantMetricTensor.g22()
            * (differential_operators->DDZ(covariantMetricTensor.g12())
               + differential_operators->DDX(covariantMetricTensor.g23(), dx)
               - differential_operators->DDY(covariantMetricTensor.g13(), dy))
      // + 0.5 *g23*(differential_operators->DDZ(covariantMetricTensor.Getg31()) + differential_operators->DDX(covariantMetricTensor.g33()) - differential_operators->DDZ(g_13));
      // which equals
      + 0.5 * contravariantMetricTensor.g23()
            * differential_operators->DDX(covariantMetricTensor.g33(), dx);
  G2_23_ = 0.5 * contravariantMetricTensor.g12()
               * (differential_operators->DDZ(covariantMetricTensor.g12())
                  + differential_operators->DDY(covariantMetricTensor.g13(), dy)
                  - differential_operators->DDX(covariantMetricTensor.g23(), dx))
           + 0.5 * contravariantMetricTensor.g22()
                 * differential_operators->DDZ(covariantMetricTensor.g22())
           + 0.5 * contravariantMetricTensor.g23()
                 * differential_operators->DDY(covariantMetricTensor.g33(), dy);

  G3_11_ = 0.5 * contravariantMetricTensor.g13()
               * differential_operators->DDX(covariantMetricTensor.g11(), dx)
           + contravariantMetricTensor.g23()
                 * (differential_operators->DDX(covariantMetricTensor.g12(), dx)
                    - 0.5 * differential_operators->DDY(covariantMetricTensor.g11(), dy))
           + contravariantMetricTensor.g33()
                 * (differential_operators->DDX(covariantMetricTensor.g13(), dx)
                    - 0.5 * differential_operators->DDZ(covariantMetricTensor.g11()));
  G3_22_ = contravariantMetricTensor.g13()
               * (differential_operators->DDY(covariantMetricTensor.g12(), dy)
                  - 0.5 * differential_operators->DDX(covariantMetricTensor.g22(), dx))
           + 0.5 * contravariantMetricTensor.g23()
                 * differential_operators->DDY(covariantMetricTensor.g22(), dy)
           + contravariantMetricTensor.g33()
                 * (differential_operators->DDY(covariantMetricTensor.g23(), dy)
                    - 0.5 * differential_operators->DDZ(covariantMetricTensor.g22()));
  G3_33_ = contravariantMetricTensor.g13()
               * (differential_operators->DDZ(covariantMetricTensor.g13())
                  - 0.5 * differential_operators->DDX(covariantMetricTensor.g33(), dx))
           + contravariantMetricTensor.g23()
                 * (differential_operators->DDZ(covariantMetricTensor.g23())
                    - 0.5 * differential_operators->DDY(covariantMetricTensor.g33(), dy))
           + 0.5 * contravariantMetricTensor.g33()
                 * differential_operators->DDZ(covariantMetricTensor.g33());
  G3_12_ =
      // 0.5 *g31*(differential_operators->DDY(covariantMetricTensor.g11()) + differential_operators->DDX(covariantMetricTensor.Getg12()) - differential_operators->DDX(covariantMetricTensor.g12()))
      // which equals to
      0.5 * contravariantMetricTensor.g13()
          * differential_operators->DDY(covariantMetricTensor.g11(), dy)
      // + 0.5 *g32*(differential_operators->DDY(covariantMetricTensor.Getg21()) + differential_operators->DDX(covariantMetricTensor.g22()) - differential_operators->DDY(covariantMetricTensor.g12()))
      // which equals to
      + 0.5 * contravariantMetricTensor.g23()
            * differential_operators->DDX(covariantMetricTensor.g22(), dx)
      //+ 0.5 *g33*(differential_operators->DDY(covariantMetricTensor.Getg31()) + differential_operators->DDX(covariantMetricTensor.Getg32()) - differential_operators->DDZ(covariantMetricTensor.g12()));
      // which equals to
      + 0.5 * contravariantMetricTensor.g33()
            * (differential_operators->DDY(covariantMetricTensor.g13(), dy))
      + differential_operators->DDX(covariantMetricTensor.g23(), dx)
      - differential_operators->DDZ(covariantMetricTensor.g12());
  G3_13_ = 0.5 * contravariantMetricTensor.g13()
               * differential_operators->DDZ(covariantMetricTensor.g11())
           + 0.5 * contravariantMetricTensor.g23()
                 * (differential_operators->DDZ(covariantMetricTensor.g12())
                    + differential_operators->DDX(covariantMetricTensor.g23(), dx)
                    - differential_operators->DDY(covariantMetricTensor.g13(), dy))
           + 0.5 * contravariantMetricTensor.g33()
                 * differential_operators->DDX(covariantMetricTensor.g33(), dx);
  G3_23_ = 0.5 * contravariantMetricTensor.g13()
               * (differential_operators->DDZ(covariantMetricTensor.g12())
                  + differential_operators->DDY(covariantMetricTensor.g13(), dy))
           - differential_operators->DDX(covariantMetricTensor.g23(), dx)
           + 0.5 * contravariantMetricTensor.g23()
                 * differential_operators->DDZ(covariantMetricTensor.g22())
           + 0.5 * contravariantMetricTensor.g33()
                 * differential_operators->DDY(covariantMetricTensor.g33(), dy);
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
      updated_christoffel_symbols.G3_13_, updated_christoffel_symbols.G3_23_,
      updated_christoffel_symbols.G1_, updated_christoffel_symbols.G2_,
      updated_christoffel_symbols.G3_);
}

ChristoffelSymbols ChristoffelSymbols::applyToComponents(
    const std::function<const FieldMetric(const FieldMetric)>& function) const {

  const auto components_in =
      std::vector<FieldMetric>{G1_11_, G1_22_, G1_33_, G1_12_, G1_13_, G1_23_, G2_11_,
                               G2_22_, G2_33_, G2_12_, G2_13_, G2_23_, G3_11_, G3_22_,
                               G3_33_, G3_12_, G3_13_, G3_23_, G1_,    G2_,    G3_};

  FieldMetric components_out[21];

  std::transform(components_in.begin(), components_in.end(), components_out, function);
  auto [G1_11, G1_22, G1_33, G1_12, G1_13, G1_23, G2_11, G2_22, G2_33, G2_12, G2_13,
        G2_23, G3_11, G3_22, G3_33, G3_12, G3_13, G3_23, G1, G2, G3] = components_out;

  return ChristoffelSymbols(G1_11, G1_22, G1_33, G1_12, G1_13, G1_23, G2_11, G2_22, G2_33,
                            G2_12, G2_13, G2_23, G3_11, G3_22, G3_33, G3_12, G3_13, G3_23,
                            G1, G2, G3, differential_operators);
}
