#pragma once
#ifndef COORDINATES_ACCESSOR_H__
#define COORDINATES_ACCESSOR_H__

#include "array.hxx"
#include "build_config.hxx"
#include "coordinates.hxx"
#include "macro_for_each.hxx"

/// Provide (hopefully) fast access to Coordinates data
/// e.g. grid spacing, metric tensors etc.
///
/// Contains only two member variables:
///   - data     a BoutReal pointer
///   - mesh_nz  an int
///
/// data is striped, so that quantities at a given
/// grid cell are packed together.
///
/// Example
///
///   auto coord_acc = CoordinatesAccessor(mesh->getCoordinates());
///   coord_acc.dx(index)  -> BoutReal at cell index
///
/// Notes
///
///  * Data from Coordinates is copied into an array which
///    is cached. CoordinatesAccessors created with the same
///    Coordinates pointer will re-use the same array without
///    copying the data again.
///    -> If Coordinates data is changed, the cache should be cleared
///    by calling CoordinatesAccessor::clear()
struct CoordinatesAccessor {
  CoordinatesAccessor() = delete;

  /// Constructor from Coordinates
  /// Copies data from coords, doesn't modify it
  explicit CoordinatesAccessor(const Coordinates* coords);

  /// Clear the cache of Coordinates data
  ///
  /// By default this clears everything; if only a specific
  /// Coordinates should be removed then that can be specified.
  ///
  /// Returns the number of data arrays removed
  /// (mainly to assist in testing)
  static std::size_t clear(const Coordinates* coords = nullptr);

  /// Offsets of each coordinates variable into the striped array
  enum class Offset {
    dx,
    dy,
    dz, // Grid spacing
    d1_dx,
    d1_dy,
    d1_dz, // Grid spacing non-uniformity
    J,     // Jacobian
    B,
    Byup,
    Bydown, // Magnetic field magnitude
    G1,
    G3, // Metric derivatives
    g11,
    g12,
    g13,
    g22,
    g23,
    g33, // Contravariant metric tensor (g^{ij})
    g_11,
    g_12,
    g_13,
    g_22,
    g_23,
    g_33, // Covariant metric tensor
    end
  };

  /// The number of values for each grid point
  /// Note: This might be > end to align to memory boundaries
  ///       e.g. 32-byte boundaries -> multiples of 32 / 8 = 4 desirable
  static constexpr int stripe_size = static_cast<int>(Offset::end);

  static_assert(stripe_size >= static_cast<int>(Offset::end),
                "stripe_size must fit all Coordinates values");

  /// Underlying data pointer.
  /// This array includes all Coordinates fields interleaved
  BoutReal* data;
  int mesh_nz; ///< For converting from 3D to 2D index

  /// Lookup value in data array, based on the cell index
  /// and the variable offset
  BOUT_HOST_DEVICE inline BoutReal lookup(int index, int offset) const {
#if BOUT_USE_METRIC_3D
    const int ind = index; // Use 3D index
#else
    const int ind = index / mesh_nz; // Convert to a 2D index
#endif
    return data[stripe_size * ind + offset];
  }

  /// Create functions to access data e.g. dx(index)
  /// by casting the Offset to an int and passing to lookup function
  ///
  /// e.g. COORD_FN1(dx) defines a function
  ///     BoutReal dx(int index) const {...}
#define COORD_FN1(symbol)                                    \
  BOUT_HOST_DEVICE inline BoutReal symbol(int index) const { \
    return lookup(index, static_cast<int>(Offset::symbol));  \
  }

  /// Generates lookup functions for each symbol
  /// Macro should accept up to 10 arguments
  ///
  /// e.g. COORD_FN(dx, dy) produces two functions
  ///     BoutReal dx(int index) const {...}
  ///     BoutReal dy(int index) const {...}
#define COORD_FN(...) MACRO_FOR_EACH(COORD_FN1, __VA_ARGS__)

  COORD_FN(dx, dy, dz);
  COORD_FN(d1_dx, d1_dy, d1_dz);
  COORD_FN(J);
  COORD_FN(B, Byup, Bydown);
  COORD_FN(G1, G3);
  COORD_FN(g11, g12, g13, g22, g23, g33);
  COORD_FN(g_11, g_12, g_13, g_22, g_23, g_33);
};

#endif // COORDINATES_ACCESSOR_H__
