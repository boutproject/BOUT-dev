
#ifndef BOUT_METRIC_TENSOR_HXX
#define BOUT_METRIC_TENSOR_HXX

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include <bout/bout_types.hxx>
#include <functional>

class MetricTensor {

public:
#if BOUT_USE_METRIC_3D
  using FieldMetric = Field3D;
#else
  using FieldMetric = Field2D;
#endif

  MetricTensor(FieldMetric g11, FieldMetric g22, FieldMetric g33, FieldMetric g12,
               FieldMetric g13, FieldMetric g23);

  MetricTensor(BoutReal g11, BoutReal g22, BoutReal g33, BoutReal g12, BoutReal g13,
               BoutReal g23, Mesh* mesh);

  // check that tensors are positive (if expected) and finite (always)
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

#if BOUT_USE_METRIC_3D
  const BoutReal* g11(int x, int y) const { return g11_m(x, y); }
  const BoutReal* g22(int x, int y) const { return g22_m(x, y); }
  const BoutReal* g33(int x, int y) const { return g33_m(x, y); }
  const BoutReal* g12(int x, int y) const { return g12_m(x, y); }
  const BoutReal* g13(int x, int y) const { return g13_m(x, y); }
  const BoutReal* g23(int x, int y) const { return g23_m(x, y); }
#else
  const BoutReal& g11(int x, int y) const { return g11_m(x, y); }
  const BoutReal& g22(int x, int y) const { return g22_m(x, y); }
  const BoutReal& g33(int x, int y) const { return g33_m(x, y); }
  const BoutReal& g12(int x, int y) const { return g12_m(x, y); }
  const BoutReal& g13(int x, int y) const { return g13_m(x, y); }
  const BoutReal& g23(int x, int y) const { return g23_m(x, y); }
#endif

  void setMetricTensor(const MetricTensor& metric_tensor) {

    g11_m = metric_tensor.g11();
    g22_m = metric_tensor.g22();
    g33_m = metric_tensor.g33();
    g12_m = metric_tensor.g12();
    g13_m = metric_tensor.g13();
    g23_m = metric_tensor.g23();
  }

  void setLocation(const CELL_LOC location) {
    g11_m.setLocation(location);
    g22_m.setLocation(location);
    g33_m.setLocation(location);
    g12_m.setLocation(location);
    g13_m.setLocation(location);
    g23_m.setLocation(location);
  }

  MetricTensor inverse(const std::string& region = "RGN_ALL");

  // Transforms the MetricTensor by applying the given function to every component
  void map(const std::function<const FieldMetric(const FieldMetric)>& function);

  MetricTensor applyToComponents(
      const std::function<const FieldMetric(const FieldMetric)>& function) const;

  void communicate(Mesh* mesh);

private:
  FieldMetric g11_m, g22_m, g33_m, g12_m, g13_m, g23_m;
};

class CovariantMetricTensor : public MetricTensor {

public:
  CovariantMetricTensor(FieldMetric g11, FieldMetric g22, FieldMetric g33,
                        FieldMetric g12, FieldMetric g13, FieldMetric g23)
      : MetricTensor(std::move(g11), std::move(g22), std::move(g33), std::move(g12),
                     std::move(g13), std::move(g23)){};

  CovariantMetricTensor(const BoutReal g11, const BoutReal g22, const BoutReal g33,
                        const BoutReal g12, const BoutReal g13, const BoutReal g23,
                        Mesh* mesh)
      : MetricTensor(g11, g22, g33, g12, g13, g23, mesh){};
};

class ContravariantMetricTensor : public MetricTensor {

public:
  ContravariantMetricTensor(FieldMetric g_11, FieldMetric g_22, FieldMetric g_33,
                            FieldMetric g_12, FieldMetric g_13, FieldMetric g_23)
      : MetricTensor(std::move(g_11), std::move(g_22), std::move(g_33), std::move(g_12),
                     std::move(g_13), std::move(g_23)){};

  ContravariantMetricTensor(const BoutReal g_11, const BoutReal g_22, const BoutReal g_33,
                            const BoutReal g_12, const BoutReal g_13, const BoutReal g_23,
                            Mesh* mesh)
      : MetricTensor(g_11, g_22, g_33, g_12, g_13, g_23, mesh){};
};

#endif //BOUT_METRIC_TENSOR_HXX
