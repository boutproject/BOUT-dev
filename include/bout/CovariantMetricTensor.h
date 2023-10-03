
#ifndef BOUT_COVARIANTMETRICTENSOR_H
#define BOUT_COVARIANTMETRICTENSOR_H

#include "field2d.hxx"
#include "bout/field3d.hxx"

class CovariantMetricTensor {

public:
#if BOUT_USE_METRIC_3D
  using FieldMetric = Field
#else
  using FieldMetric = Field2D;
#endif

      struct CovariantComponents {
    FieldMetric g_11, g_22, g_33, g_12, g_13, g_23;
  };

  CovariantMetricTensor(const FieldMetric g_11, const FieldMetric g_22,
                        const FieldMetric g_33, const FieldMetric g_12,
                        const FieldMetric g_13, const FieldMetric g_23);

  CovariantMetricTensor(const Array<BoutReal> g_11, const Array<BoutReal> g_22,
                        const Array<BoutReal> g_33, const Array<BoutReal> g_12,
                        const Array<BoutReal> g_13, const Array<BoutReal> g_23,
                        Mesh* mesh);

  /// Invert covariant metric to get contravariant components
  int calcContravariant(const std::string& region = "RGN_ALL");

  // check that covariant tensors are positive (if expected) and finite (always)
  void checkCovariant(int ystart);

  void setCovariantMetricTensor(CovariantMetricTensor metric_tensor);

  CovariantComponents getCovariantMetricTensor();

  void Allocate();

private:
  CovariantComponents covariant_components;
};

#endif //BOUT_COVARIANTMETRICTENSOR_H
