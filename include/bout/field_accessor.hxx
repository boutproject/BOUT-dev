// GPU version Field-Accessor, updated by Dr. Yining Qin, Oct.27, 2020

#pragma once
#ifndef FIELD_ACCESSOR_H__
#define FIELD_ACCESSOR_H__

#include "../bout_types.hxx"
#include "../field.hxx"
#include "../field2d.hxx"
#include "../field3d.hxx"
#include "coordinates.hxx"

#if defined(BOUT_USE_CUDA) && defined(__CUDACC__)
#define BOUT_HOST_DEVICE __host__ __device__
#define BOUT_HOST __host__
#define BOUT_DEVICE __device__
#else
#define BOUT_HOST_DEVICE
#define BOUT_HOST
#define BOUT_DEVICE
#endif

template <CELL_LOC location = CELL_CENTRE>
struct FieldAccessor {
  explicit FieldAccessor(Field3D& f) {
    ASSERT0(f.getLocation() == location);
    coords = f.getCoordinates();

    ASSERT0(f.isAllocated());

    data = &f(0, 0, 0);

    //----- Field 3d data -> array for GPU
    f_data = f(0, 0);

    // Field size for GPU
    f_nx = f.getNx();
    f_ny = f.getNy();
    f_nz = f.getNz();

    //------ Field2D data to array for GPU
    f2d_dx = &f.getCoordinates()->dx(0, 0);
    f2d_dy = &f.getCoordinates()->dy(0, 0);
    f2d_dz = f.getCoordinates()->dz;
    f2d_J = &f.getCoordinates()->J(0, 0);

    f2d_G1 = &f.getCoordinates()->G1(0, 0);
    f2d_G3 = &f.getCoordinates()->G3(0, 0);
    f2d_g11 = &f.getCoordinates()->g11(0, 0);
    f2d_g12 = &f.getCoordinates()->g12(0, 0);
    f2d_g13 = &f.getCoordinates()->g13(0, 0);
    f2d_g22 = &f.getCoordinates()->g22(0, 0);
    f2d_g23 = &f.getCoordinates()->g23(0, 0);
    f2d_g33 = &f.getCoordinates()->g33(0, 0);

    f2d_g_11 = &f.getCoordinates()->g_11(0, 0);
    f2d_g_12 = &f.getCoordinates()->g_12(0, 0);
    f2d_g_13 = &f.getCoordinates()->g_13(0, 0);
    f2d_g_22 = &f.getCoordinates()->g_22(0, 0);
    f2d_g_23 = &f.getCoordinates()->g_23(0, 0);
    f2d_g_33 = &f.getCoordinates()->g_33(0, 0);

    //---------------------------------------

    if (f.hasParallelSlices()) {
      yup = &f.yup();
      ydown = &f.ydown();

      //----- Field3D data -> array for GPU
      f_yup = static_cast<BoutReal*>(f.yup()(0, 0));
      f_ydown = static_cast<BoutReal*>(f.ydown()(0, 0));
      // ddt() array data for GPU
      f_ddt = static_cast<BoutReal*>(f.timeDeriv()->operator()(0, 0));

      // set field region index for GPU
      // auto indices = f.getRegion("RGN_NOBNDRY").getIndices(); //set index of region
      // Ind3D *ob_i = &(indices[0]);
      //---------------------
    }
  }

  BoutReal& BOUT_HOST_DEVICE operator[](const Ind3D& d) { return data[d.ind]; }
  const BoutReal& BOUT_HOST_DEVICE operator[](const Ind3D& d) const {
    return data[d.ind];
  }

  Field3D* yup{nullptr};

  Field3D* ydown{nullptr};

  Coordinates* coords{nullptr};

  BoutReal* data;

  // Metric tensor (Coordinates) data
  // Note: The data size depends on Coordinates::FieldMetric
  //       and could be Field2D or Field3D
  BoutReal* f2d_dx = nullptr;
  BoutReal* f2d_dy = nullptr;
  BoutReal f2d_dz = 0.0;

  BoutReal* f2d_J = nullptr;

  BoutReal* f2d_G1 = nullptr;
  BoutReal* f2d_G3 = nullptr;
  BoutReal* f2d_g11 = nullptr;
  BoutReal* f2d_g12 = nullptr;
  BoutReal* f2d_g13 = nullptr;
  BoutReal* f2d_g22 = nullptr;
  BoutReal* f2d_g23 = nullptr;
  BoutReal* f2d_g33 = nullptr;

  BoutReal* f2d_g_11 = nullptr;
  BoutReal* f2d_g_12 = nullptr;
  BoutReal* f2d_g_13 = nullptr;
  BoutReal* f2d_g_22 = nullptr;
  BoutReal* f2d_g_23 = nullptr;
  BoutReal* f2d_g_33 = nullptr;

  BoutReal* f_yup = nullptr; ///< Pointer to the Field2D yup data
  BoutReal* f_ydown = nullptr; ///< Pointer to the Field2D ydown data
  BoutReal* f_data = nullptr; ///< Pointer to the Field2D data 

  BoutReal* f_ddt = nullptr;

  int* ind_3D_index = nullptr;

  int f_nx = 0;
  int f_ny = 0;
  int f_nz = 0;
};

#endif
