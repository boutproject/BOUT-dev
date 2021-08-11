// GPU version Field-Accessor, updated by Dr. Yining Qin, Oct.27, 2020

#pragma once
#ifndef FIELD_ACCESSOR_H__
#define FIELD_ACCESSOR_H__

#include "../bout_types.hxx"
#include "../field.hxx"
#include "../field2d.hxx"
#include "../field3d.hxx"
#include "coordinate_field_accessor.hxx"
#include "coordinates.hxx"
#include "build_config.hxx"

template <CELL_LOC location = CELL_CENTRE, class FieldType = Field3D>
struct FieldAccessor {
  /// Remove default constructor
  FieldAccessor() = delete;

  /// Constructor from Field3D
  explicit FieldAccessor(FieldType& f)
      : coords(f.getCoordinates()), dx(coords->dx), dy(coords->dy), dz(coords->dz),
        J(coords->J), G1(coords->G1), G3(coords->G3), g11(coords->g11), g12(coords->g12),
        g13(coords->g13), g22(coords->g22), g23(coords->g23), g33(coords->g33),
        g_11(coords->g_11), g_12(coords->g_12), g_13(coords->g_13), g_22(coords->g_22),
        g_23(coords->g_23), g_33(coords->g_33) {
    ASSERT0(f.getLocation() == location);
    ASSERT0(f.isAllocated());

    data = &f(0, 0, 0);

    //----- Field 3d data -> array
    f_data = f(0, 0);

    // Field size
    f_nx = f.getNx();
    f_ny = f.getNy();
    f_nz = f.getNz();

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


  Field3D* yup{nullptr};

  Field3D* ydown{nullptr};

  Coordinates* coords{nullptr};

  BoutReal* data;

  // Metric tensor (Coordinates) data
  // Note: The data size depends on Coordinates::FieldMetric
  //       and could be Field2D or Field3D
  BoutReal* f_yup = nullptr;   ///< Pointer to the Field2D yup data
  BoutReal* f_ydown = nullptr; ///< Pointer to the Field2D ydown data
  BoutReal* f_data = nullptr;  ///< Pointer to the Field2D data

  BoutReal* f_ddt = nullptr;

  int f_nx = 0;
  int f_ny = 0;
  int f_nz = 0;

  CoordinateFieldAccessor dx, dy, dz; /// Grid spacing
  CoordinateFieldAccessor J;          ///< Coordinate system Jacobian

  CoordinateFieldAccessor G1, G3;

  CoordinateFieldAccessor g11, g12, g13, g22, g23, g33;

  CoordinateFieldAccessor g_11, g_12, g_13, g_22, g_23, g_33;

};

/// Define a shorthand for 2D fields
template <CELL_LOC location = CELL_CENTRE>
using Field2DAccessor = FieldAccessor<location, Field2D>;
#endif
