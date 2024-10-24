
#ifndef BOUT_GVALUES_HXX
#define BOUT_GVALUES_HXX

#include "bout/metric_tensor.hxx"

using FieldMetric = MetricTensor::FieldMetric;

/// `GValues` needs renaming, when we know what the name should be
class GValues {
public:
  explicit GValues(const Coordinates& coordinates);

  const FieldMetric& G1() const { return G1_m; }
  const FieldMetric& G2() const { return G2_m; }
  const FieldMetric& G3() const { return G3_m; }

  template <class F>
  void map(F function) {
    G1_m = function(G1_m);
    G2_m = function(G2_m);
    G3_m = function(G3_m);
  }

private:
  FieldMetric G1_m, G2_m, G3_m;
};

#endif //BOUT_GVALUES_HXX
