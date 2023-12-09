
#include "bout/metricTensor.hxx"
#include "bout/output.hxx"

MetricTensor::MetricTensor(const FieldMetric& g11, const FieldMetric& g22,
                           const FieldMetric& g33, const FieldMetric& g12,
                           const FieldMetric& g13, const FieldMetric& g23)
    : g11_(g11), g22_(g22), g33_(g33), g12_(g12), g13_(g13), g23_(g23) {
  Allocate(); // Make sure metric elements are allocated //  ; TODO: Required?
}

MetricTensor::MetricTensor(const BoutReal g11, const BoutReal g22, const BoutReal g33,
                           const BoutReal g12, const BoutReal g13, const BoutReal g23,
                           Mesh* mesh)
    : g11_(g11, mesh), g22_(g22, mesh), g33_(g33, mesh), g12_(g12, mesh), g13_(g13, mesh),
      g23_(g23, mesh) {

  Allocate(); // Make sure metric elements are allocated //  ; TODO: Required?
}

const MetricTensor::FieldMetric& MetricTensor::g11() const { return g11_; }
const MetricTensor::FieldMetric& MetricTensor::g22() const { return g22_; }
const MetricTensor::FieldMetric& MetricTensor::g33() const { return g33_; }
const MetricTensor::FieldMetric& MetricTensor::g12() const { return g12_; }
const MetricTensor::FieldMetric& MetricTensor::g13() const { return g13_; }
const MetricTensor::FieldMetric& MetricTensor::g23() const { return g23_; }

void MetricTensor::setMetricTensor(const MetricTensor& metric_tensor) {

  g11_ = metric_tensor.g11();
  g22_ = metric_tensor.g22();
  g33_ = metric_tensor.g33();
  g12_ = metric_tensor.g12();
  g13_ = metric_tensor.g13();
  g23_ = metric_tensor.g23();
}

void MetricTensor::check(int ystart) {
  // Diagonal metric components should be finite
  bout::checkFinite(g11_, "g11", "RGN_NOCORNERS");
  bout::checkFinite(g22_, "g22", "RGN_NOCORNERS");
  bout::checkFinite(g33_, "g33", "RGN_NOCORNERS");
  if (g11_.hasParallelSlices() && &g11_.ynext(1) != &g11_) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(g11_.ynext(sign * dy), "g11.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g22_.ynext(sign * dy), "g22.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g33_.ynext(sign * dy), "g33.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
  // Diagonal metric components should be positive
  bout::checkPositive(g11_, "g11", "RGN_NOCORNERS");
  bout::checkPositive(g22_, "g22", "RGN_NOCORNERS");
  bout::checkPositive(g33_, "g33", "RGN_NOCORNERS");
  if (g11_.hasParallelSlices() && &g11_.ynext(1) != &g11_) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkPositive(g11_.ynext(sign * dy), "g11.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(g22_.ynext(sign * dy), "g22.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(g33_.ynext(sign * dy), "g33.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }

  // Off-diagonal metric components should be finite
  bout::checkFinite(g12_, "g12", "RGN_NOCORNERS");
  bout::checkFinite(g13_, "g13", "RGN_NOCORNERS");
  bout::checkFinite(g23_, "g23", "RGN_NOCORNERS");
  if (g23_.hasParallelSlices() && &g23_.ynext(1) != &g23_) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(g12_.ynext(sign * dy), "g12.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g13_.ynext(sign * dy), "g13.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(g23_.ynext(sign * dy), "g23.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
}

void MetricTensor::Allocate() { //  ; TODO: Required?
  g11_.allocate();
  g22_.allocate();
  g33_.allocate();
  g12_.allocate();
  g13_.allocate();
  g23_.allocate();
}

void MetricTensor::setLocation(const CELL_LOC location) {
  g11_.setLocation(location);
  g22_.setLocation(location);
  g33_.setLocation(location);
  g12_.setLocation(location);
  g13_.setLocation(location);
  g23_.setLocation(location);
}

MetricTensor MetricTensor::inverse(const std::string& region) {

  TRACE("MetricTensor::CalculateOppositeRepresentation");

  // Perform inversion of g{ij} to get g^{ij}, or vice versa
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  BOUT_FOR_SERIAL(i, g11_.getRegion(region)) {
    a(0, 0) = g11_[i];
    a(1, 1) = g22_[i];
    a(2, 2) = g33_[i];

    a(0, 1) = a(1, 0) = g12_[i];
    a(1, 2) = a(2, 1) = g23_[i];
    a(0, 2) = a(2, 0) = g13_[i];

    if (invert3x3(a)) {
      const auto error_message = "\tERROR: metric tensor is singular at ({:d}, {:d})\n";
      output_error.write(error_message, i.x(), i.y());
      throw BoutException(error_message);
    }
  }

  BoutReal g_11, g_22, g_33, g_12, g_13, g_23;
  g_11 = a(0, 0);
  g_22 = a(1, 1);
  g_33 = a(2, 2);
  g_12 = a(0, 1);
  g_13 = a(0, 2);
  g_23 = a(1, 2);

  //  BoutReal maxerr;
  //  maxerr = BOUTMAX(
  //      max(abs((g_11 * g_11 + g_12 * g_12
  //               + g_13 * g_13)
  //              - 1)),
  //      max(abs((g_12 * g_12 + g_22 * g_22
  //               + g_23 * g_23)
  //              - 1)),
  //      max(abs((g_13 * g_13 + g_23 * g_23
  //               + g_33 * g_33)
  //              - 1)));
  //
  //  output_info.write("\tMaximum error in diagonal inversion is {:e}\n", maxerr);
  //
  //  maxerr = BOUTMAX(
  //      max(abs(g_11 * g_12 + g_12 * g_22
  //              + g_13 * g_23)),
  //      max(abs(g_11 * g_13 + g_12 * g_23
  //              + g_13 * g_33)),
  //      max(abs(g_12 * g_13 + g_22 * g_23
  //              + g_23 * g_33)));
  //
  //  output_info.write("\tMaximum error in off-diagonal inversion is {:e}\n", maxerr);
  const auto mesh = g11_.getMesh(); // All the components have the same mesh
  auto other_representation = MetricTensor(g_11, g_22, g_33, g_12, g_13, g_23, mesh);
  const auto location = g11_.getLocation();
  other_representation.setLocation(location);
  return other_representation;
}

void MetricTensor::map(
    const std::function<const FieldMetric(const FieldMetric)>& function) {

  const MetricTensor updated_metric_tensor = applyToComponents(function);

  setMetricTensor(MetricTensor(updated_metric_tensor.g11_, updated_metric_tensor.g22_,
                               updated_metric_tensor.g33_, updated_metric_tensor.g12_,
                               updated_metric_tensor.g13_, updated_metric_tensor.g23_));
}

MetricTensor MetricTensor::applyToComponents(
    const std::function<const FieldMetric(const FieldMetric)>& function) const {

  const auto components_in = std::vector<FieldMetric>{g11_, g22_, g33_, g12_, g13_, g23_};

  FieldMetric components_out[6];

  std::transform(components_in.begin(), components_in.end(), components_out, function);
  auto [g_11, g_22, g_33, g_12, g_13, g_23] = components_out;

  return MetricTensor(g_11, g_22, g_33, g_12, g_13, g_23);
}
