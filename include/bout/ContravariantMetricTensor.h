
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

      struct ContravariantComponents {
    FieldMetric g11, g22, g33, g12, g13, g23;
  };

  ContravariantMetricTensor(const FieldMetric g11, const FieldMetric g22,
                            const FieldMetric g33, const FieldMetric g12,
                            const FieldMetric g13, const FieldMetric g23);

  ContravariantMetricTensor(const Array<BoutReal> g11, const Array<BoutReal> g22,
                            const Array<BoutReal> g33, const Array<BoutReal> g12,
                            const Array<BoutReal> g13, const Array<BoutReal> g23,
                            Mesh* mesh);

  /// Invert contravariant metric to get covariant components
  int calcCovariant(CELL_LOC location, const std::string& region = "RGN_ALL");

  // check that contravariant tensors are positive (if expected) and finite (always)
  void checkContravariant(int ystart);

  void setContravariantMetricTensor(CELL_LOC location,
                                    const ContravariantMetricTensor& metric_tensor);

  ContravariantComponents getContravariantMetricTensor() const;

  void Allocate();

private:
  ContravariantComponents contravariant_components;
};

#endif //BOUT_CONTRAVARIANTMETRICTENSOR_H
