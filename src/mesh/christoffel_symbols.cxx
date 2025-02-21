#include "bout/christoffel_symbols.hxx"
#include "bout/coordinates.hxx"
#include "bout/mesh.hxx"

ChristoffelSymbols::ChristoffelSymbols(Coordinates& coordinates) {
  // Calculate Christoffel symbol terms (18 independent values)
  // Note: This calculation is completely general: metric
  // tensor can be 2D or 3D. For 2D, all DDZ terms are zero

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

  G1_11_m = 0.5 * g11 * coordinates.DDX(g_11)
            + g12 * (coordinates.DDX(g_12) - 0.5 * coordinates.DDY(g_11))
            + g13 * (coordinates.DDX(g_13) - 0.5 * coordinates.DDZ(g_11));
  G1_22_m = g11 * (coordinates.DDY(g_12) - 0.5 * coordinates.DDX(g_22))
            + 0.5 * g12 * coordinates.DDY(g_22)
            + g13 * (coordinates.DDY(g_23) - 0.5 * coordinates.DDZ(g_22));
  G1_33_m = g11 * (coordinates.DDZ(g_13) - 0.5 * coordinates.DDX(g_33))
            + g12 * (coordinates.DDZ(g_23) - 0.5 * coordinates.DDY(g_33))
            + 0.5 * g13 * coordinates.DDZ(g_33);
  G1_12_m =
      0.5 * g11 * coordinates.DDY(g_11) + 0.5 * g12 * coordinates.DDX(g_22)
      + 0.5 * g13
            * (coordinates.DDY(g_13) + coordinates.DDX(g_23) - coordinates.DDZ(g_12));
  G1_13_m =
      0.5 * g11 * coordinates.DDZ(g_11)
      + 0.5 * g12
            * (coordinates.DDZ(g_12) + coordinates.DDX(g_23) - coordinates.DDY(g_13))
      + 0.5 * g13 * coordinates.DDX(g_33);
  G1_23_m =
      0.5 * g11 * (coordinates.DDZ(g_12) + coordinates.DDY(g_13) - coordinates.DDX(g_23))
      + 0.5 * g12
            * (coordinates.DDZ(g_22) + coordinates.DDY(g_23) - coordinates.DDY(g_23))
      // + 0.5 *g13*(coordinates.DDZ(g_32) + coordinates.DDY(g_33) - coordinates.DDZ(g_23));
      // which equals
      + 0.5 * g13 * coordinates.DDY(g_33);

  G2_11_m = 0.5 * g12 * coordinates.DDX(g_11)
            + g22 * (coordinates.DDX(g_12) - 0.5 * coordinates.DDY(g_11))
            + g23 * (coordinates.DDX(g_13) - 0.5 * coordinates.DDZ(g_11));
  G2_22_m = g12 * (coordinates.DDY(g_12) - 0.5 * coordinates.DDX(g_22))
            + 0.5 * g22 * coordinates.DDY(g_22)
            + g23 * (coordinates.DDY(g23) - 0.5 * coordinates.DDZ(g_22));
  G2_33_m = g12 * (coordinates.DDZ(g_13) - 0.5 * coordinates.DDX(g_33))
            + g22 * (coordinates.DDZ(g_23) - 0.5 * coordinates.DDY(g_33))
            + 0.5 * g23 * coordinates.DDZ(g_33);
  G2_12_m =
      0.5 * g12 * coordinates.DDY(g_11) + 0.5 * g22 * coordinates.DDX(g_22)
      + 0.5 * g23
            * (coordinates.DDY(g_13) + coordinates.DDX(g_23) - coordinates.DDZ(g_12));
  G2_13_m =
      // 0.5 *g21*(coordinates.DDZ(g_11) + coordinates.DDX(covariantMetricTensor.Getg13()) - coordinates.DDX(g_13))
      // which equals
      0.5 * g12 * (coordinates.DDZ(g_11) + coordinates.DDX(g_13) - coordinates.DDX(g_13))
      // + 0.5 *g22*(coordinates.DDZ(covariantMetricTensor.Getg21()) + coordinates.DDX(g_23) - coordinates.DDY(g_13))
      // which equals
      + 0.5 * g22
            * (coordinates.DDZ(g_12) + coordinates.DDX(g_23) - coordinates.DDY(g_13))
      // + 0.5 *g23*(coordinates.DDZ(covariantMetricTensor.Getg31()) + coordinates.DDX(g_33) - coordinates.DDZ(g_13));
      // which equals
      + 0.5 * g23 * coordinates.DDX(g_33);
  G2_23_m =
      0.5 * g12 * (coordinates.DDZ(g_12) + coordinates.DDY(g_13) - coordinates.DDX(g_23))
      + 0.5 * g22 * coordinates.DDZ(g_22) + 0.5 * g23 * coordinates.DDY(g_33);

  G3_11_m = 0.5 * g13 * coordinates.DDX(g_11)
            + g23 * (coordinates.DDX(g_12) - 0.5 * coordinates.DDY(g_11))
            + g33 * (coordinates.DDX(g_13) - 0.5 * coordinates.DDZ(g_11));
  G3_22_m = g13 * (coordinates.DDY(g_12) - 0.5 * coordinates.DDX(g_22))
            + 0.5 * g23 * coordinates.DDY(g_22)
            + g33 * (coordinates.DDY(g_23) - 0.5 * coordinates.DDZ(g_22));
  G3_33_m = g13 * (coordinates.DDZ(g_13) - 0.5 * coordinates.DDX(g_33))
            + g23 * (coordinates.DDZ(g_23) - 0.5 * coordinates.DDY(g_33))
            + 0.5 * g33 * coordinates.DDZ(g_33);
  G3_12_m =
      // 0.5 *g31*(coordinates.DDY(g_11) + coordinates.DDX(covariantMetricTensor.Getg12()) - coordinates.DDX(g_12))
      // which equals to
      0.5 * g13 * coordinates.DDY(g_11)
      // + 0.5 *g32*(coordinates.DDY(covariantMetricTensor.Getg21()) + coordinates.DDX(g_22) - coordinates.DDY(g_12))
      // which equals to
      + 0.5 * g23 * coordinates.DDX(g_22)
      //+ 0.5 *g33*(coordinates.DDY(covariantMetricTensor.Getg31()) + coordinates.DDX(covariantMetricTensor.Getg32()) - coordinates.DDZ(g_12));
      // which equals to
      + 0.5 * g33 * (coordinates.DDY(g_13)) + coordinates.DDX(g_23) - coordinates.DDZ(g_12);
  G3_13_m = 0.5 * g13 * coordinates.DDZ(g_11) + 0.5 * g23 * (coordinates.DDZ(g_12) + coordinates.DDX(g_23) - coordinates.DDY(g_13))
            + 0.5 * g33 * coordinates.DDX(g_33);
  G3_23_m = 0.5 * g13 * (coordinates.DDZ(g_12) + coordinates.DDY(g_13)) - coordinates.DDX(g_23) + 0.5 * g23 * coordinates.DDZ(g_22)
            + 0.5 * g33 * coordinates.DDY(g_33);

  output_progress.write("\tCommunicating connection terms\n");

  G1_11_m.getMesh()->communicate(G1_11_m, G1_22_m, G1_33_m, G1_12_m, G1_13_m, G1_23_m,
                                 G2_11_m, G2_22_m, G2_33_m, G2_12_m, G2_13_m, G2_23_m,
                                 G3_11_m, G3_22_m, G3_33_m, G3_12_m, G3_13_m, G3_23_m);
}
