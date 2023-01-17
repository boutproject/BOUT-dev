/*!
 * RAJA library utilities and wrappers
 *
 * Defines:
 *  - RajaForAll   A class which handles Region indices,
 *                 and wraps RAJA::forall
 *  - BOUT_FOR_RAJA   A macro which uses RajaForAll when BOUT_HAS_RAJA
 *                    Falls back to BOUT_FOR when RAJA not enabled
 *
 * Notes:
 *
 *  - DISABLE_RAJA can be defined to 1 before including this header,
 *    to locally disable RAJA for testing.
 */

#pragma once
#ifndef RAJALIB_H
#define RAJALIB_H

#include "bout/region.hxx"

#if BOUT_HAS_RAJA and !DISABLE_RAJA

#include "RAJA/RAJA.hpp" // using RAJA lib

/// Wrapper around RAJA::forall
/// Enables computations to be done on CPU or GPU (CUDA).
///
/// Must be constructed with a `Region`. When passed a lambda
/// function via the `<<` operator, the lambda function will be
/// called with the `Region` indices (index.ind).
///
/// Usage:
///
///   RajaForAll(f.getRegion("RGN_NOBNDRY")) << [=](int id) { /* ... */ };
///
///  where `f` is a Field
///
/// For multiple loops, the `RajaForAll` object can be created once:
///
///   RajaForAll  raja_nobndry(f.getRegion("RGN_NOBNDRY"));
///
/// and then passed lambda functions multiple times:
///
///   raja_nobndry << [=](int id) { /* ... */ };
///
struct RajaForAll {
  RajaForAll() = delete; ///< Need a Region to construct

  /// Construct by specifying a Region to iterate over.
  /// Converts the range into a form which can be used in a RAJA loop.
  ///
  /// @tparam IndType     Ind2D or Ind3D. The lambda function will
  ///                     be called with 1D (flattened) indices of this type.
  ///
  /// @param region    The region to iterate over
  ///
  template <typename IndType>
  RajaForAll(const Region<IndType>& region) {
    auto indices = region.getIndices(); // A std::vector of Ind3D objects
    _ob_i_ind.reallocate(indices.size());
    // Copy indices into Array
    for (auto i = 0; i < indices.size(); i++) {
      _ob_i_ind[i] = indices[i].ind;
    }
  }

  /// Pass a lambda function to RAJA::forall
  /// Iterates over the range passed to the constructor
  ///
  /// Returns a reference to `this`, so that `<<` can be chained.
  ///
  /// @tparam F  The function type. Expected to take an `int` input
  ///            e.g. Lambda(int) -> void
  ///
  /// @param f   Lambda function to call each iteration
  ///
  template <typename F>
  const RajaForAll& operator<<(F f) const {
    // Get the raw pointer to use on the device
    // Note: must be a local variable
    const int* _ob_i_ind_raw = &_ob_i_ind[0];
    RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, _ob_i_ind.size()),
                           [=] RAJA_DEVICE(int id) {
                             // Look up index and call user function
                             f(_ob_i_ind_raw[id]);
                           });
    return *this;
  }

private:
  Array<int> _ob_i_ind; ///< Holds the index array
};

/// Create a variable which shadows another (has the same name)
#define SHADOW_ARG(var) var = var

/// Transform a list of variables into a list of var=var assignments
/// Useful for capturing class members in lambda function arguments.
#define CAPTURE(...) MACRO_FOR_EACH_ARG(SHADOW_ARG, __VA_ARGS__)

/// Iterate an index over a region
///
/// If BOUT_HAS_RAJA is true and DISABLE_RAJA is false, then this macro
/// uses RAJA (via RajaForAll) to perform the iteration.
///
/// Usage:
///
///   BOUT_FOR_RAJA(i, f.region("RGN_NOBNDRY")) {
///     /* ... */
///   };   //<- Note semicolon!
///
/// Note: Needs to be closed with `};` because it's a lambda function
///
/// Extra arguments can be passed after the region, and will be added
/// to the lambda capture. The intended use for this is to capture
/// class member variables, which can't be used directly in a RAJA CUDA
/// loop.
///
/// Usage:
///
///   BOUT_FOR_RAJA(i, region, CAPTURE(var1, var2)) {
///     /* ... */
///   };
///
/// which will have a lambda capture [=, var1=var1, var2=var2]
/// to create variables which shadow the class members.
///
#define BOUT_FOR_RAJA(index, region, ...) \
  RajaForAll(region) << [ =, ##__VA_ARGS__ ] RAJA_DEVICE(int index)

#else // BOUT_HAS_RAJA

#warning RAJA not enabled. BOUT_FOR_RAJA falling back to BOUT_FOR.

/// If no RAJA, BOUT_FOR_RAJA reverts to BOUT_FOR
/// Note: Redundant ';' after closing brace should be ignored by compiler
///       Ignores any additional arguments
#define BOUT_FOR_RAJA(index, region, ...) BOUT_FOR(index, region)

/// If not using RAJA, CAPTURE doesn't do anything
#define CAPTURE(...)

#endif // BOUT_HAS_RAJA
#endif // RAJALIB_H
