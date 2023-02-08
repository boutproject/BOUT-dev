/// FieldAccessor
///
/// Provides quick but unsafe access to field and coordinate system data
///

#pragma once
#ifndef FIELD_ACCESSOR_H__
#define FIELD_ACCESSOR_H__

#include "build_config.hxx"
#include "coordinates.hxx"
#include "coordinates_accessor.hxx"
#include "bout/bout_types.hxx"
#include "bout/field.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"

/// Simple wrapper around a BoutReal* 1D array
///
/// This is used to provide subscript operator [] for Ind3D
struct BoutRealArray {
  BoutReal* data;

  BoutRealArray() = delete; ///< No default constructor

  /// Wrap a BoutReal pointer.
  /// This does not take ownership
  explicit BoutRealArray(BoutReal* data) : data(data) {}

  BOUT_HOST_DEVICE inline BoutReal& operator[](int ind) { return data[ind]; }

  BOUT_HOST_DEVICE inline BoutReal& operator[](const Ind3D& ind) { return data[ind.ind]; }

  BOUT_HOST_DEVICE inline const BoutReal& operator[](int ind) const { return data[ind]; }

  BOUT_HOST_DEVICE inline const BoutReal& operator[](const Ind3D& ind) const {
    return data[ind.ind];
  }

  /// Cast operators, so can be assigned to a raw pointer
  /// Note: Not explicit, so can be cast implicitly
  operator BoutReal*() { return data; }

  operator const BoutReal*() const { return data; }
};

/// Thin wrapper around field data, for fast but unsafe access
///
/// @tparam location   Cell location of the data. This will be checked on construction
/// @tparam FieldType   Either Field3D (default) or Field2D
///
template <CELL_LOC location = CELL_CENTRE, class FieldType = Field3D>
struct FieldAccessor {
  /// Remove default constructor
  FieldAccessor() = delete;

  /// Constructor from Field3D
  ///
  /// @param[in] f    The field to access. Must already be allocated
  explicit FieldAccessor(FieldType& f) : coords(f.getCoordinates()) {
    ASSERT0(f.getLocation() == location);
    ASSERT0(f.isAllocated());

    data = BoutRealArray{&f(0, 0, 0)};

    // Field size
    nx = f.getNx();
    ny = f.getNy();
    nz = f.getNz();

    // Mesh z size, for index conversion
    mesh_nz = f.getMesh()->LocalNz;

    if (f.hasParallelSlices()) {
      // Get arrays from yup and ydown fields
      yup = BoutRealArray{&(f.yup()(0, 0, 0))};
      ydown = BoutRealArray{&(f.ydown()(0, 0, 0))};
    }

    // ddt() array data
    ddt = BoutRealArray{&(f.timeDeriv()->operator()(0, 0, 0))};
  }

  /// Provide shorthand for access to field data.
  /// Does not convert between 3D and 2D indices,
  /// so fa[i] is equivalent to fa.data[i].
  ///
  BOUT_HOST_DEVICE inline const BoutReal& operator[](int ind) const { return data[ind]; }

  BOUT_HOST_DEVICE inline const BoutReal& operator[](const Ind3D& ind) const {
    return data[ind.ind];
  }

  // Pointers to the field data arrays
  // These are wrapped in BoutRealArray types so they can be indexed with Ind3D or int

  BoutRealArray data{nullptr}; ///< Pointer to the Field data
  BoutRealArray ddt{nullptr};  ///< Time-derivative data

  BoutRealArray yup{nullptr};   ///< Pointer to the Field yup data
  BoutRealArray ydown{nullptr}; ///< Pointer to the Field ydown data

  CoordinatesAccessor coords; ///< Provides access to Coordinates data

  // Field size
  int nx = 0;
  int ny = 0;
  int nz = 0;

  // Mesh Z size. Used to convert 3D to 2D indices
  int mesh_nz;
};

/// Define a shorthand for 2D fields
template <CELL_LOC location = CELL_CENTRE>
using Field2DAccessor = FieldAccessor<location, Field2D>;

/// Syntactic sugar for time derivative of a field
///
/// Usage:
///
///   ddt(fa)[i] =
///
///  where fa is a FieldAccessor, and i is an int
///
template <CELL_LOC location, class FieldType>
BOUT_HOST_DEVICE inline BoutRealArray& ddt(const FieldAccessor<location, FieldType>& fa) {
  // Note: FieldAccessor captured by value is const in RAJA kernel.
  //       Need to cast to non-const so that ddt() data can be assigned to
  return const_cast<BoutRealArray&>(fa.ddt);
}

#endif
