// GPU version Field-Accessor, updated by Dr. Yining Qin, Oct.27, 2020

#pragma once
#ifndef FIELD2D_ACCESSOR_H__
#define FIELD2D_ACCESSOR_H__

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
struct Field2DAccessor {
  explicit Field2DAccessor(Field2D& f) {
    ASSERT0(f.getLocation() == location);
    ASSERT0(f.isAllocated());


    data = &f(0, 0);
    //----- Field 2d data -> array
    f_data = &f(0, 0);

    // Field size
    f_nx = f.getNx();
    f_ny = f.getNy();
    f_nz = f.getNz();

    // Get coordinate system information
    coords = f.getCoordinates();

    //------ Field2D data to array
    f2d_dx = &coords->dx(0, 0);
    f2d_dy = &coords->dy(0, 0);
    f2d_dz = &coords->dz(0, 0);

    f2d_J = &coords->J(0, 0);

    f2d_G1 = &coords->G1(0, 0);
    f2d_G3 = &coords->G3(0, 0);
    f2d_g11 = &coords->g11(0, 0);
    f2d_g12 = &coords->g12(0, 0);
    f2d_g13 = &coords->g13(0, 0);
    f2d_g22 = &coords->g22(0, 0);
    f2d_g23 = &coords->g23(0, 0);
    f2d_g33 = &coords->g33(0, 0);

    f2d_g_11 = &coords->g_11(0, 0);
    f2d_g_12 = &coords->g_12(0, 0);
    f2d_g_13 = &coords->g_13(0, 0);
    f2d_g_22 = &coords->g_22(0, 0);
    f2d_g_23 = &coords->g_23(0, 0);
    f2d_g_33 = &coords->g_33(0, 0);

    //---------------------------------------

    if (f.hasParallelSlices()) {
      yup = &f.yup();
      ydown = &f.ydown();

      // set field region index for GPU
      // auto indices = f.getRegion("RGN_NOBNDRY").getIndices(); //set index of region
      // Ind3D *ob_i = &(indices[0]);
      //---------------------
    }
  }

  BoutReal& BOUT_HOST_DEVICE operator[](const Ind2D& d) { return data[d.ind]; }
  const BoutReal& BOUT_HOST_DEVICE operator[](const Ind2D& d) const {
    return data[d.ind];
  }

  Field2D* yup{nullptr};

  Field2D* ydown{nullptr};

  Coordinates* coords{nullptr};

  BoutReal* data;

  // Metric tensor (Coordinates) data
  // Note: The data size depends on Coordinates::FieldMetric
  //       and could be Field2D or Field3D
  BoutReal* f2d_dx = nullptr;
  BoutReal* f2d_dy = nullptr;
  BoutReal* f2d_dz = nullptr;

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

  BoutReal* f_yup = nullptr;
  BoutReal* f_ydown = nullptr;
  BoutReal* f_data = nullptr;

  BoutReal* f_ddt = nullptr;

  int* ind_2D_index = nullptr;

  int f_nx = 0;
  int f_ny = 0;
  int f_nz = 0;
  //--------------------------------
};

#endif
