
#include "bout/CovariantMetricTensor.hxx"
#include "bout/ContravariantMetricTensor.hxx"
#include "bout/coordinates.hxx"
#include "bout/output.hxx"

CovariantMetricTensor::CovariantMetricTensor(
    const FieldMetric g_11, const FieldMetric g_22, const FieldMetric g_33,
    const FieldMetric g_12, const FieldMetric g_13, const FieldMetric g_23)
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
