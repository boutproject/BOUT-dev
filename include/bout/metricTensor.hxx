
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

  void Allocate();

  void setLocation(CELL_LOC location);

  MetricTensor inverse(const std::string& region = "RGN_ALL");

  // Transforms the MetricTensor by applying the given function to every component
  void map(const std::function<const Field2D(const FieldMetric)>& function);

  MetricTensor applyToComponents(
      const std::function<const FieldMetric(const FieldMetric)>& function) const;

protected:
  FieldMetric g11_, g22_, g33_, g12_, g13_, g23_;
};

#endif //BOUT_METRICTENSOR_HXX
