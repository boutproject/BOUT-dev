#pragma once

#include "bout/boundary_region.hxx"
#include "bout/boundary_region_iter.hxx"
#include "bout/boutexception.hxx"
#include "bout/field_data.hxx"
#include "bout/mask.hxx"
#include "bout/region.hxx"

#include <concepts>
#include <memory>
#include <vector>

class Options;

namespace bout {
namespace boundary {
/// This class allows to simplify iterating over y-boundaries.
///
/// It makes it easier to write code for FieldAligned boundaries, but if a bit
/// care is taken the code also works with FluxCoordinateIndependent code.
///
/// An example how to replace old code is given here:
/// ../../manual/sphinx/user_docs/boundary_options.rst
class YBoundary {
public:
  YBoundary(YBndryType type, Options* options_ptr, const Mesh& mesh);

  /// Iterate over the boundary.
  /// This function takes a lamda / templated function, that applies the boundary on the given point.
  /// The function must take a `const BoundaryIterator auto& point` as its only argument.
  /// See also the documentation at ../../manual/sphinx/user_docs/boundary_options.rst
  template <class F>
    requires std::regular_invocable<F&, const BoundaryRegionY::Iterator&>
             || std::regular_invocable<F&, const BoundaryRegionFCI::Point&>
  void iter(F func) const {
    for (const auto& region : boundary_regions) {
      for (const auto& point : *region) {
        func(point);
      }
    }
    for (const auto& region : boundary_regions_par) {
      for (const auto& point : *region) {
        func(point);
      }
    }
  }

  /// Return true if this boundary in the lower direction contains point \p ind
  bool contains_low(const Ind3D& ind) const { return _contains_low[ind]; }
  /// Return true if this boundary in the upper direction contains point \p ind
  bool contains_high(const Ind3D& ind) const { return _contains_high[ind]; }

  /// Return true if this boundary in direction \p dir contains point \p ind
  template <int dir>
  bool contains(const Ind3D& ind) const {
    static_assert(dir == 1 || dir == -1);
    if constexpr (dir == 1) {
      return _contains_high[ind];
    }
    if constexpr (dir == -1) {
      return _contains_low[ind];
    }
  }

  /// Return true if this boundary in direction \p dir contains point \p ind
  bool contains(int dir, const Ind3D& ind) const {
    if (dir == 1) {
      return contains<+1>(ind);
    }
    if (dir == -1) {
      return contains<-1>(ind);
    }
    throw BoutException("only dir == 1 and dir == -1 are implemented, not {}", dir);
  }

private:
  std::vector<std::shared_ptr<BoundaryRegionFCI>> boundary_regions_par;
  std::vector<std::shared_ptr<BoundaryRegionY>> boundary_regions;

  BoutMask _contains_low;
  BoutMask _contains_high;
};

} // namespace boundary
} // namespace bout
