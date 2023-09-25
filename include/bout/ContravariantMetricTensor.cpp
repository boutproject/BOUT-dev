
#include "ContravariantMetricTensor.h"


ContravariantMetricTensor ContravariantMetricTensor::getContravariantMetricTensor() const {
  ContravariantMetricTensor g_contravariant = { g11, g22, g33, g12, g13, g23 };
  return g_contravariant;
}

void ContravariantMetricTensor::setContravariantMetricTensor(
    const ContravariantMetricTensor& metric_tensor) {
  g11 = metric_tensor.g11;
  g22 = metric_tensor.g22;
  g33 = metric_tensor.g33;
  g12 = metric_tensor.g12;
  g13 = metric_tensor.g13;
  g23 = metric_tensor.g23;
  calcCovariant();
}

int ContravariantMetricTensor::calcCovariant(const std::string& region) {
  TRACE("Coordinates::calcCovariant");

  // Make sure metric elements are allocated
  g_11.allocate();
  g_22.allocate();
  g_33.allocate();
  g_12.allocate();
  g_13.allocate();
  g_23.allocate();

  g_11.setLocation(location);
  g_22.setLocation(location);
  g_33.setLocation(location);
  g_12.setLocation(location);
  g_13.setLocation(location);
  g_23.setLocation(location);

  // Perform inversion of g^{ij} to get g_{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  BOUT_FOR_SERIAL(i, g11.getRegion(region)) {
    a(0, 0) = g11[i];
    a(1, 1) = g22[i];
    a(2, 2) = g33[i];

    a(0, 1) = a(1, 0) = g12[i];
    a(1, 2) = a(2, 1) = g23[i];
    a(0, 2) = a(2, 0) = g13[i];

    if (invert3x3(a)) {
      output_error.write("\tERROR: metric tensor is singular at ({:d}, {:d})\n", i.x(),
                         i.y());
      return 1;
    }

    g_11[i] = a(0, 0);
    g_22[i] = a(1, 1);
    g_33[i] = a(2, 2);

    g_12[i] = a(0, 1);
    g_13[i] = a(0, 2);
    g_23[i] = a(1, 2);
  }

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((g_11 * g11 + g_12 * g12 + g_13 * g13) - 1)),
                   max(abs((g_12 * g12 + g_22 * g22 + g_23 * g23) - 1)),
                   max(abs((g_13 * g13 + g_23 * g23 + g_33 * g33) - 1)));

  output_info.write("\tLocal maximum error in diagonal inversion is {:e}\n", maxerr);

  maxerr = BOUTMAX(max(abs(g_11 * g12 + g_12 * g22 + g_13 * g23)),
                   max(abs(g_11 * g13 + g_12 * g23 + g_13 * g33)),
                   max(abs(g_12 * g13 + g_22 * g23 + g_23 * g33)));

  output_info.write("\tLocal maximum error in off-diagonal inversion is {:e}\n", maxerr);

  return 0;
}

void ContravariantMetricTensor::checkContravariant() {
  // Diagonal metric components should be finite
  bout::checkFinite(g11, "g11", "RGN_NOCORNERS");
  bout::checkFinite(g22, "g22", "RGN_NOCORNERS");
  bout::checkFinite(g33, "g33", "RGN_NOCORNERS");
  if (g11.hasParallelSlices() && &g11.ynext(1) != &g11) {
    for (int dy = 1; dy <= localmesh->ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(g11.ynext(sign * dy), "g11.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g22.ynext(sign * dy), "g22.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g33.ynext(sign * dy), "g33.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
  // Diagonal metric components should be positive
  bout::checkPositive(g11, "g11", "RGN_NOCORNERS");
  bout::checkPositive(g22, "g22", "RGN_NOCORNERS");
  bout::checkPositive(g33, "g33", "RGN_NOCORNERS");
  if (g11.hasParallelSlices() && &g11.ynext(1) != &g11) {
    for (int dy = 1; dy <= localmesh->ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkPositive(g11.ynext(sign * dy), "g11.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(g22.ynext(sign * dy), "g22.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(g33.ynext(sign * dy), "g33.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }

  // Off-diagonal metric components should be finite
  bout::checkFinite(g12, "g12", "RGN_NOCORNERS");
  bout::checkFinite(g13, "g13", "RGN_NOCORNERS");
  bout::checkFinite(g23, "g23", "RGN_NOCORNERS");
  if (g23.hasParallelSlices() && &g23.ynext(1) != &g23) {
    for (int dy = 1; dy <= localmesh->ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(g12.ynext(sign * dy), "g12.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g13.ynext(sign * dy), "g13.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g23.ynext(sign * dy), "g23.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
}
