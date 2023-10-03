
#include "CovariantMetricTensor.h"
#include "ContravariantMetricTensor.h"
#include <bout/coordinates.hxx>
#include <bout/output.hxx>

CovariantMetricTensor::CovariantMetricTensor(
    const FieldMetric g_11, const FieldMetric g_22, const FieldMetric g_33,
    const FieldMetric g_12, const FieldMetric g_13, const FieldMetric g_23)
    : covariant_components({FieldMetric(std::move(g_11)), FieldMetric(std::move(g_22)),
                            FieldMetric(std::move(g_33)), FieldMetric(std::move(g_12)),
                            FieldMetric(std::move(g_13)), FieldMetric(std::move(g_23))}) {

  Allocate(); // Make sure metric elements are allocated //  ; TODO: Required?
}

CovariantMetricTensor::CovariantMetricTensor(const Array<BoutReal> g_11,
                                             const Array<BoutReal> g_22,
                                             const Array<BoutReal> g_33,
                                             const Array<BoutReal> g_12,
                                             const Array<BoutReal> g_13,
                                             const Array<BoutReal> g_23, Mesh* mesh)
    : covariant_components({FieldMetric(g_11, mesh), FieldMetric(g_22, mesh),
                            FieldMetric(g_33, mesh), FieldMetric(g_12, mesh),
                            FieldMetric(g_13, mesh), FieldMetric(g_23, mesh)}) {

  Allocate(); // Make sure metric elements are allocated //  ; TODO: Required?
}

CovariantMetricTensor::CovariantComponents
CovariantMetricTensor::getCovariantMetricTensor() {
  return CovariantComponents{covariant_components.g_11, covariant_components.g_22,
                             covariant_components.g_33, covariant_components.g_12,
                             covariant_components.g_13, covariant_components.g_23};
}

void CovariantMetricTensor::setCovariantMetricTensor(
    CovariantMetricTensor& metric_tensor) {

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

  // Perform inversion of g_{ij} to get g^{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  BOUT_FOR_SERIAL(i, covariant_components.g_11.getRegion(region)) {
    a(0, 0) = covariant_components.g_11[i];
    a(1, 1) = covariant_components.g_22[i];
    a(2, 2) = covariant_components.g_33[i];

    a(0, 1) = a(1, 0) = covariant_components.g_12[i];
    a(1, 2) = a(2, 1) = covariant_components.g_23[i];
    a(0, 2) = a(2, 0) = covariant_components.g_13[i];

    if (invert3x3(a)) {
      output_error.write("\tERROR: metric tensor is singular at ({:d}, {:d})\n", i.x(),
                         i.y());
      return 1;
    }
  }

  ContravariantMetricTensor const contravariantMetricTensor =
      ContravariantMetricTensor(a(0, 0), a(1, 1), a(2, 2), a(0, 1), a(0, 2), a(1, 2));

  auto const contravariant_components =
      contravariantMetricTensor.getContravariantMetricTensor();

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

  output_info.write("\tMaximum error in diagonal inversion is {:e}\n", maxerr);

  maxerr = BOUTMAX(max(abs(covariant_components.g_11 * contravariant_components.g12
                           + covariant_components.g_12 * contravariant_components.g22
                           + covariant_components.g_13 * contravariant_components.g23)),
                   max(abs(covariant_components.g_11 * contravariant_components.g13
                           + covariant_components.g_12 * contravariant_components.g23
                           + covariant_components.g_13 * contravariant_components.g33)),
                   max(abs(covariant_components.g_12 * contravariant_components.g13
                           + covariant_components.g_22 * contravariant_components.g23
                           + covariant_components.g_23 * contravariant_components.g33)));

  output_info.write("\tMaximum error in off-diagonal inversion is {:e}\n", maxerr);
  return 0;
}

void CovariantMetricTensor::checkCovariant(int ystart) {
  // Diagonal metric components should be finite
  bout::checkFinite(covariant_components.g_11, "g_11", "RGN_NOCORNERS");
  bout::checkFinite(covariant_components.g_22, "g_22", "RGN_NOCORNERS");
  bout::checkFinite(covariant_components.g_33, "g_33", "RGN_NOCORNERS");
  if (covariant_components.g_11.hasParallelSlices()
      && &covariant_components.g_11.ynext(1) != &covariant_components.g_11) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(covariant_components.g_11.ynext(sign * dy), "g_11.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(covariant_components.g_22.ynext(sign * dy), "g_22.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(covariant_components.g_33.ynext(sign * dy), "g_33.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
  // Diagonal metric components should be positive
  bout::checkPositive(covariant_components.g_11, "g_11", "RGN_NOCORNERS");
  bout::checkPositive(covariant_components.g_22, "g_22", "RGN_NOCORNERS");
  bout::checkPositive(covariant_components.g_33, "g_33", "RGN_NOCORNERS");
  if (covariant_components.g_11.hasParallelSlices()
      && &covariant_components.g_11.ynext(1) != &covariant_components.g_11) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkPositive(covariant_components.g_11.ynext(sign * dy), "g_11.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(covariant_components.g_22.ynext(sign * dy), "g_22.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkPositive(covariant_components.g_33.ynext(sign * dy), "g_33.ynext",
                            fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }

  // Off-diagonal metric components should be finite
  bout::checkFinite(covariant_components.g_12, "g_12", "RGN_NOCORNERS");
  bout::checkFinite(covariant_components.g_13, "g_13", "RGN_NOCORNERS");
  bout::checkFinite(covariant_components.g_23, "g_23", "RGN_NOCORNERS");
  if (covariant_components.g_23.hasParallelSlices()
      && &covariant_components.g_23.ynext(1) != &covariant_components.g_23) {
    for (int dy = 1; dy <= ystart; ++dy) {
      for (const auto sign : {1, -1}) {
        bout::checkFinite(covariant_components.g_12.ynext(sign * dy), "g_12.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(covariant_components.g_13.ynext(sign * dy), "g_13.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
        bout::checkFinite(covariant_components.g_23.ynext(sign * dy), "g_23.ynext",
                          fmt::format("RGN_YPAR_{:+d}", sign * dy));
      }
    }
  }
}

void CovariantMetricTensor::Allocate() { //  ; TODO: Required?
  covariant_components.g_11.allocate();
  covariant_components.g_22.allocate();
  covariant_components.g_33.allocate();
  covariant_components.g_12.allocate();
  covariant_components.g_13.allocate();
  covariant_components.g_23.allocate();
}
