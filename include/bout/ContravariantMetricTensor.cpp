
#include "ContravariantMetricTensor.h"
#include "CovariantMetricTensor.h"

ContravariantMetricTensor::ContravariantMetricTensor(
    const FieldMetric g11, const FieldMetric g22, const FieldMetric g33,
    const FieldMetric g12, const FieldMetric g13, const FieldMetric g23) {
  contravariant_components = {g11, g22, g33, g12, g13, g23};
  Allocate(); // Make sure metric elements are allocated //  ; TODO: Required?
}

ContravariantMetricTensor::ContravariantComponents
ContravariantMetricTensor::getContravariantMetricTensor() const {
  return ContravariantComponents{
      contravariant_components.g11, contravariant_components.g22,
      contravariant_components.g33, contravariant_components.g12,
      contravariant_components.g13, contravariant_components.g23};
}

void ContravariantMetricTensor::setContravariantMetricTensor(
    CELL_LOC location, const ContravariantMetricTensor& metric_tensor) {

  const auto new_components = metric_tensor.getContravariantMetricTensor();
  contravariant_components.g11 = new_components.g11;
  contravariant_components.g22 = new_components.g22;
  contravariant_components.g33 = new_components.g33;
  contravariant_components.g12 = new_components.g12;
  contravariant_components.g13 = new_components.g13;
  contravariant_components.g23 = new_components.g23;
  calcCovariant(location);
}

int ContravariantMetricTensor::calcCovariant(const CELL_LOC location,
                                             std::string& region) {
  TRACE("ContravariantMetricTensor::calcCovariant");

  // Perform inversion of g^{ij} to get g_{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  BOUT_FOR_SERIAL(i, contravariant_components.g11.getRegion(region)) {
    a(0, 0) = contravariant_components.g11[i];
    a(1, 1) = contravariant_components.g22[i];
    a(2, 2) = contravariant_components.g33[i];

    a(0, 1) = a(1, 0) = contravariant_components.g12[i];
    a(1, 2) = a(2, 1) = contravariant_components.g23[i];
    a(0, 2) = a(2, 0) = contravariant_components.g13[i];

    if (invert3x3(a)) {
      output_error.write("\tERROR: metric tensor is singular at ({:d}, {:d})\n", i.x(),
                         i.y());
      return 1;
    }
  }

  CovariantMetricTensor const covariantMetricTensor =
      CovariantMetricTensor(a(0, 0), a(1, 1), a(2, 2), a(0, 1), a(0, 2), a(1, 2));

  auto const covariant_components = covariantMetricTensor.getCovariantMetricTensor();

  covariant_components.g_11.setLocation(location);
  covariant_components.g_22.setLocation(location);
  covariant_components.g_33.setLocation(location);
  covariant_components.g_12.setLocation(location);
  covariant_components.g_13.setLocation(location);
  covariant_components.g_23.setLocation(location);

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((covariant_components.g_11 * contravariant_components.g11
                            + covariant_components.g_12 * contravariant_components.g12
                            + covariant_components.g_13 * contravariant_components.g13)
                           - 1)),
                   max(abs((covariant_components.g_12 * contravariant_components.g12
                            + covariant_components.g_22 * contravariant_components.g22
                            + covariant_components.g_23 * contravariant_components.g23)
                           - 1)),
                   max(abs((covariant_components.g_13 * contravariant_components.g13
                            + covariant_components.g_23 * contravariant_components.g23
                            + covariant_components.g_33 * contravariant_components.g33)
                           - 1)));

  output_info.write("\tLocal maximum error in diagonal inversion is {:e}\n", maxerr);

  maxerr = BOUTMAX(max(abs(covariant_components.g_11 * contravariant_components.g12
                           + covariant_components.g_12 * contravariant_components.g22
                           + covariant_components.g_13 * contravariant_components.g23)),
                   max(abs(covariant_components.g_11 * contravariant_components.g13
                           + covariant_components.g_12 * contravariant_components.g23
                           + covariant_components.g_13 * contravariant_components.g33)),
                   max(abs(covariant_components.g_12 * contravariant_components.g13
                           + covariant_components.g_22 * contravariant_components.g23
                           + covariant_components.g_23 * contravariant_components.g33)));

  output_info.write("\tLocal maximum error in off-diagonal inversion is {:e}\n", maxerr);

  return 0;
}

void ContravariantMetricTensor::checkContravariant(int ystart) {
  // Diagonal metric components should be finite
  bout::checkFinite(contravariant_components.g11, "g11", "RGN_NOCORNERS");
  bout::checkFinite(contravariant_components.g22, "g22", "RGN_NOCORNERS");
  bout::checkFinite(contravariant_components.g33, "g33", "RGN_NOCORNERS");
  if (contravariant_components.g11.hasParallelSlices()
      && &contravariant_components.g11.ynext(1) != &contravariant_components.g11) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(contravariant_components.g11.ynext(sign * dy), "g11.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(contravariant_components.g22.ynext(sign * dy), "g22.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(contravariant_components.g33.ynext(sign * dy), "g33.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
  // Diagonal metric components should be positive
  bout::checkPositive(contravariant_components.g11, "g11", "RGN_NOCORNERS");
  bout::checkPositive(contravariant_components.g22, "g22", "RGN_NOCORNERS");
  bout::checkPositive(contravariant_components.g33, "g33", "RGN_NOCORNERS");
  if (contravariant_components.g11.hasParallelSlices()
      && &contravariant_components.g11.ynext(1) != &contravariant_components.g11) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkPositive(contravariant_components.g11.ynext(sign * dy), "g11.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(contravariant_components.g22.ynext(sign * dy), "g22.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(contravariant_components.g33.ynext(sign * dy), "g33.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }

  // Off-diagonal metric components should be finite
  bout::checkFinite(contravariant_components.g12, "g12", "RGN_NOCORNERS");
  bout::checkFinite(contravariant_components.g13, "g13", "RGN_NOCORNERS");
  bout::checkFinite(contravariant_components.g23, "g23", "RGN_NOCORNERS");
  if (contravariant_components.g23.hasParallelSlices()
      && &contravariant_components.g23.ynext(1) != &contravariant_components.g23) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(contravariant_components.g12.ynext(sign * dy), "g12.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(contravariant_components.g13.ynext(sign * dy), "g13.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(contravariant_components.g23.ynext(sign * dy), "g23.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
}

void ContravariantMetricTensor::Allocate() { //  ; TODO: Required?
  contravariant_components.g11.allocate();
  contravariant_components.g22.allocate();
  contravariant_components.g33.allocate();
  contravariant_components.g12.allocate();
  contravariant_components.g13.allocate();
  contravariant_components.g23.allocate();
}
