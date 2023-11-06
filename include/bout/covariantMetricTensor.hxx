
#ifndef BOUT_COVARIANTMETRICTENSOR_HXX
#define BOUT_COVARIANTMETRICTENSOR_HXX

#include "field2d.hxx"
#include "bout/field3d.hxx"

class ContravariantMetricTensor;
class CovariantMetricTensor {

public:
#if BOUT_USE_METRIC_3D
  using FieldMetric = Field
#else
  using FieldMetric = Field2D;
#endif

  CovariantMetricTensor(const FieldMetric& g_11, const FieldMetric& g_22,
                        const FieldMetric& g_33, const FieldMetric& g_12,
                        const FieldMetric& g_13, const FieldMetric& g_23);

  CovariantMetricTensor(const BoutReal g_11, const BoutReal g_22, const BoutReal g_33,
                        const BoutReal g_12, const BoutReal g_13, const BoutReal g_23,
                        Mesh* mesh);

  // check that covariant tensors are positive (if expected) and finite (always)
  void checkCovariant(int ystart);

  const FieldMetric& Getg_11() const;
  const FieldMetric& Getg_22() const;
  const FieldMetric& Getg_33() const;
  const FieldMetric& Getg_12() const;
  const FieldMetric& Getg_13() const;
  const FieldMetric& Getg_23() const;

  void setCovariantMetricTensor(const CovariantMetricTensor& metric_tensor);

  void Allocate();

  void setLocation(const CELL_LOC location);

  /// Invert contravariant metric to get covariant components
  void calcCovariant(ContravariantMetricTensor contravariantMetricTensor,
                     CELL_LOC location, const std::string& region = "RGN_ALL");

private:
  FieldMetric g_11, g_22, g_33, g_12, g_13, g_23;
};

#endif //BOUT_COVARIANTMETRICTENSOR_HXX
