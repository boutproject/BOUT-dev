
#include "bout/metricTensor.hxx"

MetricTensor::MetricTensor(const FieldMetric& g11, const FieldMetric& g22,
                           const FieldMetric& g33, const FieldMetric& g12,
                           const FieldMetric& g13, const FieldMetric& g23)
    : g11(std::move(g11)), g22(std::move(g22)), g33(std::move(g33)), g12(std::move(g12)),
      g13(std::move(g13)), g23(std::move(g23)) {

  Allocate(); // Make sure metric elements are allocated //  ; TODO: Required?
}

MetricTensor::MetricTensor(const BoutReal g11, const BoutReal g22, const BoutReal g33,
                           const BoutReal g12, const BoutReal g13, const BoutReal g23,
                           Mesh* mesh)
    : g11(g11, mesh), g22(g22, mesh), g33(g33, mesh), g12(g12, mesh), g13(g13, mesh),
      g23(g23, mesh) {

  Allocate(); // Make sure metric elements are allocated //  ; TODO: Required?
}

const MetricTensor::FieldMetric& MetricTensor::Getg11() const { return g11; }
const MetricTensor::FieldMetric& MetricTensor::Getg22() const { return g22; }
const MetricTensor::FieldMetric& MetricTensor::Getg33() const { return g33; }
const MetricTensor::FieldMetric& MetricTensor::Getg12() const { return g12; }
const MetricTensor::FieldMetric& MetricTensor::Getg13() const { return g13; }
const MetricTensor::FieldMetric& MetricTensor::Getg23() const { return g23; }

void MetricTensor::setMetricTensor(const MetricTensor& metric_tensor) {

  g11 = metric_tensor.Getg11();
  g22 = metric_tensor.Getg22();
  g33 = metric_tensor.Getg33();
  g12 = metric_tensor.Getg12();
  g13 = metric_tensor.Getg13();
  g23 = metric_tensor.Getg23();
}

void MetricTensor::check(int ystart) {
  // Diagonal metric components should be finite
  bout::checkFinite(g11, "g11", "RGN_NOCORNERS");
  bout::checkFinite(g22, "g22", "RGN_NOCORNERS");
  bout::checkFinite(g33, "g33", "RGN_NOCORNERS");
  if (g11.hasParallelSlices() && &g11.ynext(1) != &g11) {
    for (int dy = 1; dy <= ystart; ++dy) {
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
    for (int dy = 1; dy <= ystart; ++dy) {
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
    for (int dy = 1; dy <= ystart; ++dy) {
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

void MetricTensor::Allocate() { //  ; TODO: Required?
  g11.allocate();
  g22.allocate();
  g33.allocate();
  g12.allocate();
  g13.allocate();
  g23.allocate();
}

void MetricTensor::setLocation(const CELL_LOC location) {
  g11.setLocation(location);
  g22.setLocation(location);
  g33.setLocation(location);
  g12.setLocation(location);
  g13.setLocation(location);
  g23.setLocation(location);
}
