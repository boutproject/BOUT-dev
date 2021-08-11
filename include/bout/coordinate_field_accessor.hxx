#pragma once
#ifndef COORDINATE_FIELD_ACCESSOR_H__
#define COORDINATE_FIELD_ACCESSOR_H__

#include "build_config.hxx"
#include "coordinates.hxx"

/// Thin wrapper around coordinate component data
/// to provide a high performance interface
///
/// This provides an indexing operation, which accounts for the
/// difference between 2D and 3D fields
struct CoordinateFieldAccessor {
  CoordinateFieldAccessor() = delete;

  /// Construct from a Coordinates component
  explicit CoordinateFieldAccessor(Coordinates::FieldMetric& f) {
    ASSERT0(f.isAllocated());

    data = &f(0, 0, 0); // Underlying data. Note: 3 indices always works
    mesh_nz = f.getMesh()->LocalNz;
  }

  BoutReal& BOUT_HOST_DEVICE operator[](int index3D) {
#if BOUT_USE_METRIC_3D
    return data[index3D]; // A Field3D, so just use the index
#else
    return data[index3D / mesh_nz]; // Field2D, so convert index
#endif
  }

  BoutReal* data; ///< The pointer to the underlying data
  int mesh_nz;    ///< Size of the Z coordinate in the mesh. Used to convert indices
};

#endif // COORDINATE_FIELD_ACCESSOR_H__
