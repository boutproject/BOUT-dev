
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

  MetricTensor(BoutReal g11, BoutReal g22, BoutReal g33, BoutReal g12, BoutReal g13,
               BoutReal g23, Mesh* mesh);

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

  void setLocation(CELL_LOC location);

  MetricTensor oppositeRepresentation(CELL_LOC location,
                                      const std::string& region = "RGN_ALL");

  // Transforms the MetricTensor by applying the given function to every component
  void map(std::function<const Field2D(const FieldMetric)> function);

  MetricTensor
  applyToComponents(std::function<const FieldMetric(const FieldMetric)> function) const;

protected:
  FieldMetric g11, g22, g33, g12, g13, g23;
};

#endif //BOUT_METRICTENSOR_HXX