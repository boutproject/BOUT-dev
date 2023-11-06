
#ifndef BOUT_CONTRAVARIANTMETRICTENSOR_HXX
#define BOUT_CONTRAVARIANTMETRICTENSOR_HXX

#include "covariantMetricTensor.hxx"
#include "metricTensor.hxx"

class ContravariantMetricTensor : public MetricTensor {

public:
  ContravariantMetricTensor(const FieldMetric& g11, const FieldMetric& g22,
                            const FieldMetric& g33, const FieldMetric& g12,
                            const FieldMetric& g13, const FieldMetric& g23);
  ContravariantMetricTensor(const BoutReal g11, const BoutReal g22, const BoutReal g33,
                            const BoutReal g12, const BoutReal g13, const BoutReal g23,
                            Mesh* mesh);

  /// Invert covariant metric to get contravariant components
  void CalculateOppositeRepresentation(MetricTensor& covariantMetricTensor,
                                       CELL_LOC location,
                                       const std::string& region = "RGN_ALL") override;
};

#endif //BOUT_CONTRAVARIANTMETRICTENSOR_HXX
