
#ifndef BOUT_METRICTENSOR_HXX
#define BOUT_METRICTENSOR_HXX

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include <bout/bout_types.hxx>

class MetricTensor {

public:
#if BOUT_USE_METRIC_3D
  using FieldMetric = Field
#else
  using FieldMetric = Field2D;
#endif

  MetricTensor(FieldMetric g11, FieldMetric g22, FieldMetric g33, FieldMetric g12,
               FieldMetric g13, FieldMetric g23);

  MetricTensor(BoutReal g11, BoutReal g22, BoutReal g33, BoutReal g12, BoutReal g13,
               BoutReal g23, Mesh* mesh);

  // check that tensors are positive (if expected) and finite (always)
  void check(int ystart);

  const FieldMetric& g11() const;
  const FieldMetric& g22() const;
  const FieldMetric& g33() const;
  const FieldMetric& g12() const;
  const FieldMetric& g13() const;
  const FieldMetric& g23() const;

  void setMetricTensor(const MetricTensor& metric_tensor);

  void setLocation(CELL_LOC location);

  MetricTensor inverse(const std::string& region = "RGN_ALL");

  // Transforms the MetricTensor by applying the given function to every component
  void map(const std::function<const Field2D(const FieldMetric)>& function);

  MetricTensor applyToComponents(
      const std::function<const FieldMetric(const FieldMetric)>& function) const;

protected:
  FieldMetric g11_, g22_, g33_, g12_, g13_, g23_;
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

#endif //BOUT_METRICTENSOR_HXX
