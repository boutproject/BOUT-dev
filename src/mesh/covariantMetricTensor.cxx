
#include "bout/covariantMetricTensor.hxx"
#include "bout/contravariantMetricTensor.hxx"
#include "bout/coordinates.hxx"
#include "bout/output.hxx"

CovariantMetricTensor::CovariantMetricTensor(
    const FieldMetric& g_11, const FieldMetric& g_22, const FieldMetric& g_33,
    const FieldMetric& g_12, const FieldMetric& g_13, const FieldMetric& g_23)
    : MetricTensor::MetricTensor(g_11, g_22, g_33, g_12, g_13, g_23) {}

CovariantMetricTensor::CovariantMetricTensor(const BoutReal g_11, const BoutReal g_22,
                                             const BoutReal g_33, const BoutReal g_12,
                                             const BoutReal g_13, const BoutReal g_23,
                                             Mesh* mesh)
    : MetricTensor::MetricTensor(g_11, g_22, g_33, g_12, g_13, g_23, mesh) {}

void CovariantMetricTensor::calcCovariant(
    ContravariantMetricTensor contravariantMetricTensor, const CELL_LOC location,
    const std::string& region) {
  TRACE("CovariantMetricTensor::calcCovariant");

  // Perform inversion of g^{ij} to get g{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  BOUT_FOR_SERIAL(i, contravariantMetricTensor.Getg11().getRegion(region)) {
    a(0, 0) = contravariantMetricTensor.Getg11()[i];
    a(1, 1) = contravariantMetricTensor.Getg22()[i];
    a(2, 2) = contravariantMetricTensor.Getg33()[i];

    a(0, 1) = a(1, 0) = contravariantMetricTensor.Getg12()[i];
    a(1, 2) = a(2, 1) = contravariantMetricTensor.Getg23()[i];
    a(0, 2) = a(2, 0) = contravariantMetricTensor.Getg13()[i];

    if (invert3x3(a)) {
      const auto error_message = "\tERROR: metric tensor is singular at ({:d}, {:d})\n";
      output_error.write(error_message, i.x(), i.y());
      throw BoutException(error_message);
    }
  }

  g11 = a(0, 0);
  g22 = a(1, 1);
  g33 = a(2, 2);
  g12 = a(0, 1);
  g13 = a(0, 2);
  g23 = a(1, 2);
  //  covariant_components = CovariantComponents{(a(0, 0), a(1, 1), a(2, 2), a(0, 1), a(0, 2), a(1, 2))};

  setLocation(location);

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((g11 * contravariantMetricTensor.Getg11()
                            + g12 * contravariantMetricTensor.Getg12()
                            + g13 * contravariantMetricTensor.Getg13())
                           - 1)),
                   max(abs((g12 * contravariantMetricTensor.Getg12()
                            + g22 * contravariantMetricTensor.Getg22()
                            + g23 * contravariantMetricTensor.Getg23())
                           - 1)),
                   max(abs((g13 * contravariantMetricTensor.Getg13()
                            + g23 * contravariantMetricTensor.Getg23()
                            + g33 * contravariantMetricTensor.Getg33())
                           - 1)));

  output_info.write("\tLocal maximum error in diagonal inversion is {:e}\n", maxerr);

  maxerr = BOUTMAX(max(abs(g11 * contravariantMetricTensor.Getg12()
                           + g12 * contravariantMetricTensor.Getg22()
                           + g13 * contravariantMetricTensor.Getg23())),
                   max(abs(g11 * contravariantMetricTensor.Getg13()
                           + g12 * contravariantMetricTensor.Getg23()
                           + g13 * contravariantMetricTensor.Getg33())),
                   max(abs(g12 * contravariantMetricTensor.Getg13()
                           + g22 * contravariantMetricTensor.Getg23()
                           + g23 * contravariantMetricTensor.Getg33())));

  output_info.write("\tLocal maximum error in off-diagonal inversion is {:e}\n", maxerr);
}
