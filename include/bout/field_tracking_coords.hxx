#pragma once

#include <bout/field_tracking.hxx>
#include <bout/field_metric.hxx>

namespace{
#if BOUT_USE_METRIC_3D
  using FieldMetric = Field3D;
#else
  using FieldMetric = Field2D;
#endif
}

// TODO: actually do things
class TrackingFieldNoOp : public FieldTracking<FieldMetric, TrackingFieldNoOp> {
public:
  using FieldTracking::FieldTracking;
  //using FieldMetric::FieldMetric;
  void onChange(){};
};
