
#ifndef BOUT_GVALUES_HXX
#define BOUT_GVALUES_HXX

#include "bout/metricTensor.hxx"

using FieldMetric = MetricTensor::FieldMetric;

/// `GValues` needs renaming, when we know what the name should be
class GValues {

public:
  GValues(FieldMetric G1, FieldMetric G2, FieldMetric G3);

  explicit GValues(const Coordinates& coordinates);

  const FieldMetric& G1() const;
  const FieldMetric& G2() const;
  const FieldMetric& G3() const;

  void setG1(const FieldMetric& G1);
  void setG2(const FieldMetric& G2);
  void setG3(const FieldMetric& G3);

private:
  FieldMetric G1_, G2_, G3_;
};

#endif //BOUT_GVALUES_HXX
