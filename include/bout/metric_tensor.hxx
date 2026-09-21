#ifndef BOUT_METRIC_TENSOR_HXX
#define BOUT_METRIC_TENSOR_HXX

#include <bout/bout_types.hxx>
#include <bout/build_defines.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>

#include <string>

namespace bout {
#if BOUT_USE_METRIC_3D
using FieldMetric = Field3D;
using FieldMetricParallel = Field3DParallel;
#else
using FieldMetric = Field2D;
using FieldMetricParallel = Field2D;
#endif
} // namespace bout

class Coordinates;
struct MetricNormaliser;

class MetricTensor {
public:
  friend class Coordinates;

#if BOUT_USE_METRIC_3D
  using Metric2DSlice = const BoutReal*;
#else
  using Metric2DSlice = const BoutReal&;
#endif
  using FieldMetric = bout::FieldMetric;

  MetricTensor(const MetricTensor&) = default;
  MetricTensor(MetricTensor&&) = default;
  MetricTensor& operator=(const MetricTensor&) = default;
  MetricTensor& operator=(MetricTensor&&) = default;
  MetricTensor(FieldMetric g11, FieldMetric g22, FieldMetric g33, FieldMetric g12,
               FieldMetric g13, FieldMetric g23);

  MetricTensor(BoutReal g11, BoutReal g22, BoutReal g33, BoutReal g12, BoutReal g13,
               BoutReal g23, Mesh* mesh);
  virtual ~MetricTensor() = default;

  /// Check that tensors are positive (if expected) and finite (always)
  void check(int ystart);

  const FieldMetric& g11() const { return g11_m; }
  const FieldMetric& g22() const { return g22_m; }
  const FieldMetric& g33() const { return g33_m; }
  const FieldMetric& g12() const { return g12_m; }
  const FieldMetric& g13() const { return g13_m; }
  const FieldMetric& g23() const { return g23_m; }

  const BoutReal& g11(int x, int y, int z) const { return g11_m(x, y, z); }
  const BoutReal& g22(int x, int y, int z) const { return g22_m(x, y, z); }
  const BoutReal& g33(int x, int y, int z) const { return g33_m(x, y, z); }
  const BoutReal& g12(int x, int y, int z) const { return g12_m(x, y, z); }
  const BoutReal& g13(int x, int y, int z) const { return g13_m(x, y, z); }
  const BoutReal& g23(int x, int y, int z) const { return g23_m(x, y, z); }

  Metric2DSlice g11(int x, int y) const { return g11_m(x, y); }
  Metric2DSlice g22(int x, int y) const { return g22_m(x, y); }
  Metric2DSlice g33(int x, int y) const { return g33_m(x, y); }
  Metric2DSlice g12(int x, int y) const { return g12_m(x, y); }
  Metric2DSlice g13(int x, int y) const { return g13_m(x, y); }
  Metric2DSlice g23(int x, int y) const { return g23_m(x, y); }

  /// Transforms the MetricTensor by applying the given function to every component
  template <class F>
  void map(F function) {
    g11_m = function(g11_m);
    g22_m = function(g22_m);
    g33_m = function(g33_m);
    g12_m = function(g12_m);
    g13_m = function(g13_m);
    g23_m = function(g23_m);
  }

  void communicate();

  template <class F>
  void normaliseMetric(const MetricNormaliser& norm, const F& op);

private:
  FieldMetric g11_m, g22_m, g33_m, g12_m, g13_m, g23_m;
};

class CovariantMetricTensor;
class ContravariantMetricTensor;

class CovariantMetricTensor : public MetricTensor {
public:
  using MetricTensor::MetricTensor;

  auto inverse(const std::string& region = "RGN_ALL", bool communicate = true)
      -> ContravariantMetricTensor;

  void normaliseMetric(const MetricNormaliser& norm);
};

class ContravariantMetricTensor : public MetricTensor {
public:
  using MetricTensor::MetricTensor;

  auto inverse(const std::string& region = "RGN_ALL", bool communicate = true)
      -> CovariantMetricTensor;

  void normaliseMetric(const MetricNormaliser& norm);
};

#endif //BOUT_METRIC_TENSOR_HXX
