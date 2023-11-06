
#ifndef BOUT_COVARIANTMETRICTENSOR_HXX
#define BOUT_COVARIANTMETRICTENSOR_HXX

#include "metricTensor.hxx"

class ContravariantMetricTensor;
class CovariantMetricTensor : public MetricTensor {

public:
  CovariantMetricTensor(const FieldMetric& g_11, const FieldMetric& g_22,
                        const FieldMetric& g_33, const FieldMetric& g_12,
                        const FieldMetric& g_13, const FieldMetric& g_23);
  CovariantMetricTensor(const BoutReal g_11, const BoutReal g_22, const BoutReal g_33,
                        const BoutReal g_12, const BoutReal g_13, const BoutReal g_23,
                        Mesh* mesh);

  /// Invert contravariant metric to get covariant components
  void CalculateOppositeRepresentation(MetricTensor& contravariantMetricTensor,
                                       const CELL_LOC location,
                                       const std::string& region = "RGN_ALL") override;
};

#endif //BOUT_COVARIANTMETRICTENSOR_HXX
