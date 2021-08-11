// GPU version Field-Accessor, updated by Dr. Yining Qin, Oct.27, 2020

#pragma once
#ifndef FIELD_ACCESSOR_H__
#define FIELD_ACCESSOR_H__

#include "../bout_types.hxx"
#include "../field.hxx"
#include "../field2d.hxx"
#include "../field3d.hxx"
#include "coordinates.hxx"
#include "build_config.hxx"


template <CELL_LOC location = CELL_CENTRE>
struct FieldAccessor {
  /// Remove default constructor
  FieldAccessor() = delete;

  /// Constructor from Field3D
  explicit FieldAccessor(Field3D& f) {
    ASSERT0(f.getLocation() == location);
    ASSERT0(f.isAllocated());

    data = &f(0, 0, 0);

    //----- Field 3d data -> array
    f_data = f(0, 0);

    // Field size
    f_nx = f.getNx();
    f_ny = f.getNy();
    f_nz = f.getNz();

    // Get coordinate system information
    coords = f.getCoordinates();

    // Get pointers to the underlying arrays
    // Note that data may be 2D or 3D, depending on Coordinate::FieldMetric
    // In either case we can use three indices to access the data
    dx = &coords->dx(0, 0, 0);
    dy = &coords->dy(0, 0, 0);
    dz = &coords->dz(0, 0, 0);

    J = &coords->J(0, 0, 0);

    G1 = &coords->G1(0, 0, 0);
    G3 = &coords->G3(0, 0, 0);
    g11 = &coords->g11(0, 0, 0);
    g12 = &coords->g12(0, 0, 0);
    g13 = &coords->g13(0, 0, 0);
    g22 = &coords->g22(0, 0, 0);
    g23 = &coords->g23(0, 0, 0);
    g33 = &coords->g33(0, 0, 0);

    g_11 = &coords->g_11(0, 0, 0);
    g_12 = &coords->g_12(0, 0, 0);
    g_13 = &coords->g_13(0, 0, 0);
    g_22 = &coords->g_22(0, 0, 0);
    g_23 = &coords->g_23(0, 0, 0);
    g_33 = &coords->g_33(0, 0, 0);

    if (f.hasParallelSlices()) {
      yup = &f.yup();
      ydown = &f.ydown();

      // Field3D data -> array
      f_yup = &f.yup()(0, 0, 0);
      f_ydown = &(f.ydown()(0, 0, 0));
    }

    // ddt() array data
    f_ddt = &(f.timeDeriv()->operator()(0, 0, 0));
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
  BoutReal* dx = nullptr;
  BoutReal* dy = nullptr;
  BoutReal* dz = nullptr;

  BoutReal* J = nullptr; ///< Coordinate system Jacobian

  BoutReal* G1 = nullptr;
  BoutReal* G3 = nullptr;
  BoutReal* g11 = nullptr;
  BoutReal* g12 = nullptr;
  BoutReal* g13 = nullptr;
  BoutReal* g22 = nullptr;
  BoutReal* g23 = nullptr;
  BoutReal* g33 = nullptr;

  BoutReal* g_11 = nullptr;
  BoutReal* g_12 = nullptr;
  BoutReal* g_13 = nullptr;
  BoutReal* g_22 = nullptr;
  BoutReal* g_23 = nullptr;
  BoutReal* g_33 = nullptr;

  BoutReal* f_yup = nullptr;   ///< Pointer to the Field2D yup data
  BoutReal* f_ydown = nullptr; ///< Pointer to the Field2D ydown data
  BoutReal* f_data = nullptr;  ///< Pointer to the Field2D data

  BoutReal* f_ddt = nullptr;

  int f_nx = 0;
  int f_ny = 0;
  int f_nz = 0;
};

#endif
