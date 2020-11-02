#pragma once
#ifndef FIELD_ACCESSOR_H__
#define FIELD_ACCESSOR_H__

#include "../bout_types.hxx"
#include "../field3d.hxx"
#include "../field.hxx"
#include "../field2d.hxx"



template<CELL_LOC location = CELL_CENTRE>
struct FieldAccessor {
  explicit FieldAccessor(Field3D &f) {
    ASSERT0(f.getLocation() == location);


    coords = f.getCoordinates();
    
    ASSERT0(f.isAllocated());
  
    data = &f(0,0,0);
    
    if (f.hasParallelSlices()) {
      yup = &f.yup();
      ydown = &f.ydown();
    }
  }


  BoutReal& __host__ __device__ operator[](const Ind3D &d) {
    return data[d.ind];
  }
  const BoutReal& __host__ __device__ operator[](const Ind3D &d) const {
    return data[d.ind];
  }
  

 Field3D   *yup {nullptr};

 Field3D   *ydown {nullptr};

 Coordinates*  coords {nullptr};

  
  BoutReal* data;
};

#endif 
