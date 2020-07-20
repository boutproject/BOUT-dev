#pragma once
#ifndef FIELD_ACCESSOR_H__
#define FIELD_ACCESSOR_H__

#include "../bout_types.hxx"
#include "../field3d.hxx"

template<CELL_LOC location = CELL_CENTRE>
struct FieldAccessor {
  explicit FieldAccessor(Field3D &f) {
    ASSERT0(f.getLocation() == location);

    // Get raw pointers
    coords = f.getCoordinates();

    ASSERT0(f.isAllocated());
    
    // Underlying data array
    data = &f(0,0,0);
    
    if (f.hasParallelSlices()) {
      yup = &f.yup();
      ydown = &f.ydown();
    }
  }

  BoutReal& operator[](const Ind3D &d) {
    return data[d.ind];
  }
  const BoutReal& operator[](const Ind3D &d) const {
    return data[d.ind];
  }
  
  /// These are for fast (inner loop) access to yup/ydown
  /// Used in single_index_ops.hxx functions
  Field3D *yup {nullptr}, *ydown {nullptr};

  Coordinates* coords {nullptr};

  /// Internal data array
  BoutReal* data;
};

#endif // FIELD_ACCESSOR_H__
