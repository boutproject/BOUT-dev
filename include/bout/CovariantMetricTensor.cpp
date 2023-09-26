
#include "CovariantMetricTensor.h"

CovariantMetricTensor::CovariantMetricTensor(
    const FieldMetric g_11, const FieldMetric g_22, const FieldMetric g_33,
    const FieldMetric g_12, const FieldMetric g_13, const FieldMetric g_23) {
  covariant_components = {g_11, g_22, g_33, g_12, g_13, g_23};
}

CovariantMetricTensor::CovariantComponents
CovariantMetricTensor::getCovariantMetricTensor() const {
  return CovariantComponents{covariant_components.g_11, covariant_components.g_22,
                             covariant_components.g_33, covariant_components.g_12,
                             covariant_components.g_13, covariant_components.g_23};
}

void CovariantMetricTensor::setCovariantMetricTensor(
    const CovariantMetricTensor& metric_tensor) {

  const auto new_components = metric_tensor.getCovariantMetricTensor();
  covariant_components.g_11 = new_components.g_11;
  covariant_components.g_22 = new_components.g_22;
  covariant_components.g_33 = new_components.g_33;
  covariant_components.g_12 = new_components.g_12;
  covariant_components.g_13 = new_components.g_13;
  covariant_components.g_23 = new_components.g_23;
  calcContravariant();
}

int CovariantMetricTensor::calcContravariant(const std::string& region) {
  TRACE("CovariantMetricTensor::calcContravariant");

  // Make sure metric elements are allocated
  g11.allocate();
  g22.allocate();
  g33.allocate();
  g12.allocate();
  g13.allocate();
  g23.allocate();

  // Perform inversion of g_{ij} to get g^{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  BOUT_FOR_SERIAL(i, g_11.getRegion(region)) {
    a(0, 0) = g_11[i];
    a(1, 1) = g_22[i];
    a(2, 2) = g_33[i];

    a(0, 1) = a(1, 0) = g_12[i];
    a(1, 2) = a(2, 1) = g_23[i];
    a(0, 2) = a(2, 0) = g_13[i];

    if (invert3x3(a)) {
      output_error.write("\tERROR: metric tensor is singular at ({:d}, {:d})\n", i.x(),
                         i.y());
      return 1;
    }

    g11[i] = a(0, 0);
    g22[i] = a(1, 1);
    g33[i] = a(2, 2);

    g12[i] = a(0, 1);
    g13[i] = a(0, 2);
    g23[i] = a(1, 2);
  }

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((g_11 * g11 + g_12 * g12 + g_13 * g13) - 1)),
                   max(abs((g_12 * g12 + g_22 * g22 + g_23 * g23) - 1)),
                   max(abs((g_13 * g13 + g_23 * g23 + g_33 * g33) - 1)));

  output_info.write("\tMaximum error in diagonal inversion is {:e}\n", maxerr);

  maxerr = BOUTMAX(max(abs(g_11 * g12 + g_12 * g22 + g_13 * g23)),
                   max(abs(g_11 * g13 + g_12 * g23 + g_13 * g33)),
                   max(abs(g_12 * g13 + g_22 * g23 + g_23 * g33)));

  output_info.write("\tMaximum error in off-diagonal inversion is {:e}\n", maxerr);
  return 0;
}

void CovariantMetricTensor::checkCovariant() {
  // Diagonal metric components should be finite
  bout::checkFinite(g_11, "g_11", "RGN_NOCORNERS");
  bout::checkFinite(g_22, "g_22", "RGN_NOCORNERS");
  bout::checkFinite(g_33, "g_33", "RGN_NOCORNERS");
  if (g_11.hasParallelSlices() && &g_11.ynext(1) != &g_11) {
    for (int dy = 1; dy <= localmesh->ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(g_11.ynext(sign * dy), "g_11.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g_22.ynext(sign * dy), "g_22.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g_33.ynext(sign * dy), "g_33.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
  // Diagonal metric components should be positive
  bout::checkPositive(g_11, "g_11", "RGN_NOCORNERS");
  bout::checkPositive(g_22, "g_22", "RGN_NOCORNERS");
  bout::checkPositive(g_33, "g_33", "RGN_NOCORNERS");
  if (g_11.hasParallelSlices() && &g_11.ynext(1) != &g_11) {
    for (int dy = 1; dy <= localmesh->ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkPositive(g_11.ynext(sign * dy), "g_11.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(g_22.ynext(sign * dy), "g_22.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(g_33.ynext(sign * dy), "g_33.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }

  // Off-diagonal metric components should be finite
  bout::checkFinite(g_12, "g_12", "RGN_NOCORNERS");
  bout::checkFinite(g_13, "g_13", "RGN_NOCORNERS");
  bout::checkFinite(g_23, "g_23", "RGN_NOCORNERS");
  if (g_23.hasParallelSlices() && &g_23.ynext(1) != &g_23) {
    for (int dy = 1; dy <= localmesh->ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(g_12.ynext(sign * dy), "g_12.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g_13.ynext(sign * dy), "g_13.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g_23.ynext(sign * dy), "g_23.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
}
