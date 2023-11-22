
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

  MetricTensor(const FieldMetric& g11, const FieldMetric& g22, const FieldMetric& g33,
               const FieldMetric& g12, const FieldMetric& g13, const FieldMetric& g23);

  MetricTensor(const BoutReal g11, const BoutReal g22, const BoutReal g33,
               const BoutReal g12, const BoutReal g13, const BoutReal g23, Mesh* mesh);

  // check that tensors are positive (if expected) and finite (always)
  void check(int ystart);

  const FieldMetric& Getg11() const;
  const FieldMetric& Getg22() const;
  const FieldMetric& Getg33() const;
  const FieldMetric& Getg12() const;
  const FieldMetric& Getg13() const;
  const FieldMetric& Getg23() const;

  void setMetricTensor(const MetricTensor& metric_tensor);

  void Allocate();

  void setLocation(const CELL_LOC location);

  MetricTensor oppositeRepresentation(const CELL_LOC location, Mesh* mesh,
                                      const std::string& region = "RGN_ALL");

  // Transforms the MetricTensor by applying the given function to every component
  void map(const std::function<const Field2D(const FieldMetric)> function);

  MetricTensor applyToComponents(
      const std::function<const FieldMetric(const FieldMetric)> function) const;

protected:
  FieldMetric g11, g22, g33, g12, g13, g23;

  std::vector<FieldMetric> getComponents() const;
};

#endif //BOUT_METRICTENSOR_HXX
