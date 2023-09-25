
#ifndef BOUT_CONTRAVARIANTMETRICTENSOR_H
#define BOUT_CONTRAVARIANTMETRICTENSOR_H

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/paralleltransform.hxx"
#include "bout/utils.hxx"
#include <bout/bout_types.hxx>

class ContravariantMetricTensor {

public:
#if BOUT_USE_METRIC_3D
  using FieldMetric = Field
#else
  using FieldMetric = Field2D;
#endif

  ContravariantMetricTensor(const FieldMetric g11, const FieldMetric g22,
                            const FieldMetric g33, const FieldMetric g12,
                            const FieldMetric g13, const FieldMetric g23)
      : g11(g11),
        g22(g22), g33(g33), g12(g12), g13(g13), g23(g23) {
  }

  /// Invert contravariant metric to get covariant components
  int calcCovariant(const std::string& region = "RGN_ALL");

private:
  FieldMetric g11, g22, g33, g12, g13, g23;

  // check that contravariant tensors are positive (if expected) and finite (always)
  void checkContravariant();

public:
  void setContravariantMetricTensor(const ContravariantMetricTensor& metric_tensor);
  ContravariantMetricTensor getContravariantMetricTensor() const;
};

#endif //BOUT_CONTRAVARIANTMETRICTENSOR_H
