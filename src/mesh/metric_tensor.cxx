#include "bout/metric_tensor.hxx"
#include "fmt/core.h"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field2d.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"

#include <cstdlib>
#include <utility>

MetricTensor::MetricTensor(FieldMetric g11, FieldMetric g22, FieldMetric g33,
                           FieldMetric g12, FieldMetric g13, FieldMetric g23)
    : g11_m(std::move(g11)), g22_m(std::move(g22)), g33_m(std::move(g33)),
      g12_m(std::move(g12)), g13_m(std::move(g13)), g23_m(std::move(g23)) {}

MetricTensor::MetricTensor(const BoutReal g11, const BoutReal g22, const BoutReal g33,
                           const BoutReal g12, const BoutReal g13, const BoutReal g23,
                           Mesh* mesh)
    : g11_m(g11, mesh), g22_m(g22, mesh), g33_m(g33, mesh), g12_m(g12, mesh),
      g13_m(g13, mesh), g23_m(g23, mesh) {}

void MetricTensor::check(int ystart) {
  const bool non_identity_parallel_transform =
      g11_m.hasParallelSlices() && &g11_m.ynext(1) != &g11_m;

  // Diagonal metric components should be finite
  bout::checkFinite(g11_m, "g11", "RGN_NOCORNERS");
  bout::checkFinite(g22_m, "g22", "RGN_NOCORNERS");
  bout::checkFinite(g33_m, "g33", "RGN_NOCORNERS");
  if (non_identity_parallel_transform) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        const auto region = fmt::format("RGN_YPAR_{:+d}", sign * dy);
        bout::checkFinite(g11_m.ynext(sign * dy), "g11.ynext", region);
        bout::checkFinite(g22_m.ynext(sign * dy), "g22.ynext", region);
        bout::checkFinite(g33_m.ynext(sign * dy), "g33.ynext", region);
      }
    }
  }

  // Diagonal metric components should be positive
  bout::checkPositive(g11_m, "g11", "RGN_NOCORNERS");
  bout::checkPositive(g22_m, "g22", "RGN_NOCORNERS");
  bout::checkPositive(g33_m, "g33", "RGN_NOCORNERS");
  if (non_identity_parallel_transform) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        const auto region = fmt::format("RGN_YPAR_{:+d}", sign * dy);
        bout::checkPositive(g11_m.ynext(sign * dy), "g11.ynext", region);
        bout::checkPositive(g22_m.ynext(sign * dy), "g22.ynext", region);
        bout::checkPositive(g33_m.ynext(sign * dy), "g33.ynext", region);
      }
    }
  }

  // Off-diagonal metric components should be finite
  bout::checkFinite(g12_m, "g12", "RGN_NOCORNERS");
  bout::checkFinite(g13_m, "g13", "RGN_NOCORNERS");
  bout::checkFinite(g23_m, "g23", "RGN_NOCORNERS");
  if (non_identity_parallel_transform) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        const auto region = fmt::format("RGN_YPAR_{:+d}", sign * dy);
        bout::checkFinite(g12_m.ynext(sign * dy), "g12.ynext", region);
        bout::checkFinite(g13_m.ynext(sign * dy), "g13.ynext", region);
        bout::checkFinite(g23_m.ynext(sign * dy), "g23.ynext", region);
      }
    }
  }
}

MetricTensor MetricTensor::inverse(const std::string& region) {

  TRACE("MetricTensor::inverse");

  // Perform inversion of g{ij} to get g^{ij}, or vice versa
  auto matrix = Matrix<BoutReal>(3, 3);

  FieldMetric g_11 = emptyFrom(g11_m);
  FieldMetric g_22 = emptyFrom(g22_m);
  FieldMetric g_33 = emptyFrom(g33_m);
  FieldMetric g_12 = emptyFrom(g12_m);
  FieldMetric g_13 = emptyFrom(g13_m);
  FieldMetric g_23 = emptyFrom(g23_m);

  BOUT_FOR_SERIAL(i, g11_m.getRegion(region)) {
    matrix(0, 0) = g11_m[i];
    matrix(1, 1) = g22_m[i];
    matrix(2, 2) = g33_m[i];

    matrix(0, 1) = matrix(1, 0) = g12_m[i];
    matrix(1, 2) = matrix(2, 1) = g23_m[i];
    matrix(0, 2) = matrix(2, 0) = g13_m[i];

    if (invert3x3(matrix) != 0) {
      const auto error_message = fmt::format(
          "\tERROR: metric tensor is singular at ({:d}, {:d})\n", i.x(), i.y());
      output_error.write(error_message);
      throw BoutException(error_message);
    }

    g_11[i] = matrix(0, 0);
    g_22[i] = matrix(1, 1);
    g_33[i] = matrix(2, 2);
    g_12[i] = matrix(0, 1);
    g_13[i] = matrix(0, 2);
    g_23[i] = matrix(1, 2);
  }

  const BoutReal diagonal_maxerr =
      BOUTMAX(max(abs((g_11 * g_11 + g_12 * g_12 + g_13 * g_13) - 1)),
              max(abs((g_12 * g_12 + g_22 * g_22 + g_23 * g_23) - 1)),
              max(abs((g_13 * g_13 + g_23 * g_23 + g_33 * g_33) - 1)));

  output_info.write("\tMaximum error in diagonal inversion is {:e}\n", diagonal_maxerr);

  const BoutReal off_diagonal_maxerr =
      BOUTMAX(max(abs(g_11 * g_12 + g_12 * g_22 + g_13 * g_23)),
              max(abs(g_11 * g_13 + g_12 * g_23 + g_13 * g_33)),
              max(abs(g_12 * g_13 + g_22 * g_23 + g_23 * g_33)));

  output_info.write("\tMaximum error in off-diagonal inversion is {:e}\n",
                    off_diagonal_maxerr);

  g_11.getMesh()->communicate(g_11, g_22, g_33, g_12, g_13, g_23);
  return MetricTensor(g_11, g_22, g_33, g_12, g_13, g_23);
}
