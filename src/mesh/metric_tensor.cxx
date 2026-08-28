#include "bout/metric_tensor.hxx"
#include "invert3x3.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field2d.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"

#include <fmt/format.h>

#include <cstdlib>
#include <string>
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
  // Check off-diagonal separately, might not have them even if we have parallel
  // slices for the diagonal components
  if (g23_m.hasParallelSlices() && &g23_m.ynext(1) != &g23_m) {
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

namespace {
template <class InverseMetric>
auto inverse_impl(const MetricTensor& metric, const std::string& region)
    -> InverseMetric {
  // Perform inversion of g{ij} to get g^{ij}, or vice versa
  auto matrix = Matrix<BoutReal>(3, 3);

  bout::FieldMetric g_11 = emptyFrom(metric.g11());
  bout::FieldMetric g_22 = emptyFrom(metric.g22());
  bout::FieldMetric g_33 = emptyFrom(metric.g33());
  bout::FieldMetric g_12 = emptyFrom(metric.g12());
  bout::FieldMetric g_13 = emptyFrom(metric.g13());
  bout::FieldMetric g_23 = emptyFrom(metric.g23());

  BOUT_FOR_SERIAL(i, metric.g11().getRegion(region)) {
    matrix(0, 0) = metric.g11()[i];
    matrix(1, 1) = metric.g22()[i];
    matrix(2, 2) = metric.g33()[i];

    matrix(0, 1) = matrix(1, 0) = metric.g12()[i];
    matrix(1, 2) = matrix(2, 1) = metric.g23()[i];
    matrix(0, 2) = matrix(2, 0) = metric.g13()[i];

    if (const auto det = bout::invert3x3(matrix); det.has_value()) {
      throw BoutException("ERROR: metric tensor is singular at ({}, {}), determinant: {}",
                          i.x(), i.y(), det.value());
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
  return InverseMetric(g_11, g_22, g_33, g_12, g_13, g_23);
}
} // namespace

auto CovariantMetricTensor::inverse(const std::string& region, bool communicate)
    -> ContravariantMetricTensor {
  auto result = inverse_impl<ContravariantMetricTensor>(*this, region);
  if (communicate) {
    result.communicate();
  }
  return result;
}

auto ContravariantMetricTensor::inverse(const std::string& region, bool communicate)
    -> CovariantMetricTensor {
  auto result = inverse_impl<CovariantMetricTensor>(*this, region);
  if (communicate) {
    result.communicate();
  }
  return result;
}

template <class F>
void MetricTensor::normaliseMetric(const MetricNormaliser& norm, const F& op) {
  if (norm.g.has_value()) {
    op(g11_m, norm.g, norm.g_mul);
    op(g22_m, norm.g, norm.g_mul);
    op(g33_m, norm.g, norm.g_mul);
    op(g12_m, norm.g, norm.g_mul);
    op(g13_m, norm.g, norm.g_mul);
    op(g23_m, norm.g, norm.g_mul);
  } else {
    op(g11_m, norm.g11, norm.g11_mul);
    op(g22_m, norm.g22, norm.g22_mul);
    op(g33_m, norm.g33, norm.g33_mul);
    op(g12_m, norm.g12, norm.g12_mul);
    op(g13_m, norm.g13, norm.g13_mul);
    op(g23_m, norm.g23, norm.g23_mul);
  }
}

void ContravariantMetricTensor::normaliseMetric(const MetricNormaliser& norm) {
  MetricTensor::normaliseMetric(norm, [](FieldMetric& f, auto fac, bool fac_mul) {
    if (fac.has_value()) {
      if (f.hasParallelSlices()) {
        if (fac_mul) {
          f.asField3DParallel() /= fac.value();
        } else {
          f.asField3DParallel() *= fac.value();
        }
      } else {
        if (fac_mul) {
          f /= fac.value();
        } else {
          f *= fac.value();
        }
      }
    }
  });
}
void CovariantMetricTensor::normaliseMetric(const MetricNormaliser& norm) {
  MetricTensor::normaliseMetric(norm, [](FieldMetric& f, auto fac, bool fac_mul) {
    if (fac.has_value()) {
      if (f.hasParallelSlices()) {
        if (fac_mul) {
          f.asField3DParallel() *= fac.value();
        } else {
          f.asField3DParallel() /= fac.value();
        }
      } else {
        if (fac_mul) {
          f *= fac.value();
        } else {
          f /= fac.value();
        }
      }
    }
  });
}

void MetricTensor::communicate() {
  g11_m.getMesh()->communicate_no_slices(g11_m, g22_m, g33_m, g12_m, g13_m, g23_m);
}
