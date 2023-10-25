
#include "bout/CovariantMetricTensor.hxx"
#include "bout/ContravariantMetricTensor.hxx"
#include "bout/coordinates.hxx"
#include "bout/output.hxx"
  
CovariantMetricTensor::CovariantMetricTensor(
    const FieldMetric& g_11, const FieldMetric& g_22, const FieldMetric& g_33,
    const FieldMetric& g_12, const FieldMetric& g_13, const FieldMetric& g_23)
    : g_11(std::move(g_11)), g_22(std::move(g_22)), g_33(std::move(g_33)),
      g_12(std::move(g_12)), g_13(std::move(g_13)), g_23(std::move(g_23)) {

  Allocate(); // Make sure metric elements are allocated //  ; TODO: Required?
}

CovariantMetricTensor::CovariantMetricTensor(const BoutReal g_11, const BoutReal g_22,
                                             const BoutReal g_33, const BoutReal g_12,
                                             const BoutReal g_13, const BoutReal g_23,
                                             Mesh* mesh)
    : g_11(g_11, mesh), g_22(g_22, mesh), g_33(g_33, mesh), g_12(g_12, mesh),
      g_13(g_13, mesh), g_23(g_23, mesh) {

  Allocate(); // Make sure metric elements are allocated //  ; TODO: Required?
}

const CovariantMetricTensor::FieldMetric& CovariantMetricTensor::Getg_11() const {
  return g_11;
}
const CovariantMetricTensor::FieldMetric& CovariantMetricTensor::Getg_22() const {
  return g_22;
}
const CovariantMetricTensor::FieldMetric& CovariantMetricTensor::Getg_33() const {
  return g_33;
}
const CovariantMetricTensor::FieldMetric& CovariantMetricTensor::Getg_12() const {
  return g_12;
}
const CovariantMetricTensor::FieldMetric& CovariantMetricTensor::Getg_13() const {
  return g_13;
}
const CovariantMetricTensor::FieldMetric& CovariantMetricTensor::Getg_23() const {
  return g_23;
}

void CovariantMetricTensor::setCovariantMetricTensor(
    const CovariantMetricTensor& metric_tensor) {

  g_11 = metric_tensor.Getg_11();
  g_22 = metric_tensor.Getg_22();
  g_33 = metric_tensor.Getg_33();
  g_12 = metric_tensor.Getg_12();
  g_13 = metric_tensor.Getg_13();
  g_23 = metric_tensor.Getg_23();
}

void CovariantMetricTensor::checkCovariant(int ystart) {
  // Diagonal metric components should be finite
  bout::checkFinite(g_11, "g_11", "RGN_NOCORNERS");
  bout::checkFinite(g_22, "g_22", "RGN_NOCORNERS");
  bout::checkFinite(g_33, "g_33", "RGN_NOCORNERS");
  if (g_11.hasParallelSlices() && &g_11.ynext(1) != &g_11) {
    for (int dy = 1; dy <= ystart; ++dy) {
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
    for (int dy = 1; dy <= ystart; ++dy) {
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
    for (int dy = 1; dy <= ystart; ++dy) {
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

void CovariantMetricTensor::Allocate() { //  ; TODO: Required?
  g_11.allocate();
  g_22.allocate();
  g_33.allocate();
  g_12.allocate();
  g_13.allocate();
  g_23.allocate();
}

void CovariantMetricTensor::setLocation(const CELL_LOC location) {
  g_11.setLocation(location);
  g_22.setLocation(location);
  g_33.setLocation(location);
  g_12.setLocation(location);
  g_13.setLocation(location);
  g_23.setLocation(location);
}

void CovariantMetricTensor::calcCovariant(
    ContravariantMetricTensor contravariantMetricTensor, const CELL_LOC location,
    const std::string& region) {
  TRACE("CovariantMetricTensor::calcCovariant");

  // Perform inversion of g^{ij} to get g_{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  const auto contravariant_components =
      contravariantMetricTensor.getContravariantMetricTensor();

  auto a = Matrix<BoutReal>(3, 3);

  BOUT_FOR_SERIAL(i, contravariant_components.g11.getRegion(region)) {
    a(0, 0) = contravariant_components.g11[i];
    a(1, 1) = contravariant_components.g22[i];
    a(2, 2) = contravariant_components.g33[i];

    a(0, 1) = a(1, 0) = contravariant_components.g12[i];
    a(1, 2) = a(2, 1) = contravariant_components.g23[i];
    a(0, 2) = a(2, 0) = contravariant_components.g13[i];

    if (invert3x3(a)) {
      const auto error_message = "\tERROR: metric tensor is singular at ({:d}, {:d})\n";
      output_error.write(error_message, i.x(), i.y());
      throw BoutException(error_message);
    }
  }

  g_11 = a(0, 0);
  g_22 = a(1, 1);
  g_33 = a(2, 2);
  g_12 = a(0, 1);
  g_13 = a(0, 2);
  g_23 = a(1, 2);
  //  covariant_components =
  //      CovariantComponents{(a(0, 0), a(1, 1), a(2, 2), a(0, 1), a(0, 2), a(1, 2))};

  setLocation(location);

  BoutReal maxerr;
  maxerr = BOUTMAX(
      max(abs((g_11 * contravariant_components.g11 + g_12 * contravariant_components.g12
               + g_13 * contravariant_components.g13)
              - 1)),
      max(abs((g_12 * contravariant_components.g12 + g_22 * contravariant_components.g22
               + g_23 * contravariant_components.g23)
              - 1)),
      max(abs((g_13 * contravariant_components.g13 + g_23 * contravariant_components.g23
               + g_33 * contravariant_components.g33)
              - 1)));

  output_info.write("\tLocal maximum error in diagonal inversion is {:e}\n", maxerr);

  maxerr = BOUTMAX(
      max(abs(g_11 * contravariant_components.g12 + g_12 * contravariant_components.g22
              + g_13 * contravariant_components.g23)),
      max(abs(g_11 * contravariant_components.g13 + g_12 * contravariant_components.g23
              + g_13 * contravariant_components.g33)),
      max(abs(g_12 * contravariant_components.g13 + g_22 * contravariant_components.g23
              + g_23 * contravariant_components.g33)));

  output_info.write("\tLocal maximum error in off-diagonal inversion is {:e}\n", maxerr);
}
