
#include "bout/contravariantMetricTensor.hxx"
#include "bout/covariantMetricTensor.hxx"
#include "bout/output.hxx"

ContravariantMetricTensor::ContravariantMetricTensor(
    const FieldMetric& g11, const FieldMetric& g22, const FieldMetric& g33,
    const FieldMetric& g12, const FieldMetric& g13, const FieldMetric& g23)
    : MetricTensor::MetricTensor(g11, g22, g33, g12, g13, g23) {}

ContravariantMetricTensor::ContravariantMetricTensor(
    const BoutReal g11, const BoutReal g22, const BoutReal g33, const BoutReal g12,
    const BoutReal g13, const BoutReal g23, Mesh* mesh)
    : MetricTensor::MetricTensor(g11, g22, g33, g12, g13, g23, mesh) {}

void ContravariantMetricTensor::calcContravariant(
    CovariantMetricTensor covariantMetricTensor, CELL_LOC location,
    const std::string& region) {

  TRACE("ContravariantMetricTensor::calcContravariant");

  // Perform inversion of g{ij} to get g^{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  BOUT_FOR_SERIAL(i, covariantMetricTensor.Getg11().getRegion(region)) {
    a(0, 0) = covariantMetricTensor.Getg11()[i];
    a(1, 1) = covariantMetricTensor.Getg22()[i];
    a(2, 2) = covariantMetricTensor.Getg33()[i];

    a(0, 1) = a(1, 0) = covariantMetricTensor.Getg12()[i];
    a(1, 2) = a(2, 1) = covariantMetricTensor.Getg23()[i];
    a(0, 2) = a(2, 0) = covariantMetricTensor.Getg13()[i];

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
  //  contravariant_components = ContravariantComponents{a(0, 0), a(1, 1), a(2, 2), a(0, 1), a(0, 2), a(1, 2)};

  setLocation(location);

  BoutReal maxerr;
  maxerr = BOUTMAX(
      max(abs((covariantMetricTensor.Getg11() * g11 + covariantMetricTensor.Getg12() * g12
               + covariantMetricTensor.Getg13() * g13)
              - 1)),
      max(abs((covariantMetricTensor.Getg12() * g12 + covariantMetricTensor.Getg22() * g22
               + covariantMetricTensor.Getg23() * g23)
              - 1)),
      max(abs((covariantMetricTensor.Getg13() * g13 + covariantMetricTensor.Getg23() * g23
               + covariantMetricTensor.Getg33() * g33)
              - 1)));

  output_info.write("\tMaximum error in diagonal inversion is {:e}\n", maxerr);

  maxerr = BOUTMAX(
      max(abs(covariantMetricTensor.Getg11() * g12 + covariantMetricTensor.Getg12() * g22
              + covariantMetricTensor.Getg13() * g23)),
      max(abs(covariantMetricTensor.Getg11() * g13 + covariantMetricTensor.Getg12() * g23
              + covariantMetricTensor.Getg13() * g33)),
      max(abs(covariantMetricTensor.Getg12() * g13 + covariantMetricTensor.Getg22() * g23
              + covariantMetricTensor.Getg23() * g33)));

  output_info.write("\tMaximum error in off-diagonal inversion is {:e}\n", maxerr);
}
