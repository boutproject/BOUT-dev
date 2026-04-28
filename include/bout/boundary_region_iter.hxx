#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <utility>

#include <bout/bout_enum_class.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/parallel_boundary_region.hxx>
#include <bout/region.hxx>
#include <bout/sys/parallel_stencils.hxx>

namespace bout {
namespace boundary {

/// Physical type of y boundary
enum class BndryType : std::int8_t {
  sheath,
  not_sheath_par,
  core,
  sol,
  sol_perp,
  sol_par,
  all,
  num
};

/// Types of free boundary condition
/// ================== ===================================================
/// Name               Description
/// ================== ===================================================
/// ``limited``        use exponential if decreasing, otherwise Neumanm
/// ``exponential``    use exponential extrapolation
/// ``linear``         use linear extrapolation
/// ================== ===================================================
enum class BoundaryFreeExtrapolation : std::int8_t { limited, exponential, linear };
//BOUT_ENUM_CLASS(BoundaryFreeExtrapolation, limited, exponential, linear);

template <typename impl>
class BoundaryRegionIterBase { /// get the index at the last point in domain
public:
  Ind3D ind() const { return static_cast<const impl*>(this)->_ind(); }
  /// get the length from the point in the domain to the boundary in index
  /// space. It is in the range [0, 1]
  BoutReal length() const { return static_cast<const impl*>(this)->_length(); }
  /// Lower bound of how many points are between the first point in the domain
  /// and the boundary in the other direction.
  signed char valid() const { return static_cast<const impl*>(this)->_valid(); }
  /// Get the width of the boundary at the current point
  int boundary_width() const { return static_cast<const impl*>(this)->_boundary_width(); }
  /// Is this the lower boundary?
  bool is_lower() const { return static_cast<const impl*>(this)->_is_lower(); }

  /*
   *         FIELD3D ACCESSORS
   */

  /// get the value at a given offset `off` of a field `f`.
  /// off = -1 is the second point in the boundary
  /// off = 0 is the first point in the boundary
  /// off = 1 is the last point in the domain
  /// off = 2 is the second to last point in the domain
  template <bool check = true>
  BoutReal& getAt(Field3D& f, int off) const {
    return static_cast<const impl*>(this)->template _getAt<check>(f, off);
  }
  /// get the value at a given offset `off` of a field `f`.
  template <bool check = true>
  BoutReal& getAt(const Field3D& f, int off) const {
    return static_cast<const impl*>(this)->template _getAt<check>(f, off);
  }

  /// Get the first point in the boundary
  const BoutReal& next(const Field3D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 0);
  }
  /// Get the first point in the boundary
  BoutReal& next(Field3D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 0);
  }
  /// Get the last point in the domain
  const BoutReal& current(const Field3D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 1);
  }
  /// Get the last point in the domain
  BoutReal& current(Field3D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 1);
  }
  /// Get the second to last point in the domain - this may not be valid and thus throw
  const BoutReal& prev(const Field3D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 2);
  }
  /*
   *         FIELD2D ACCESSORS
   */

  /// get the value at a given offset `off` of a field `f`.
  /// off = -1 is the second point in the boundary
  /// off = 0 is the first point in the boundary
  /// off = 1 is the last point in the domain
  /// off = 2 is the second to last point in the domain
  template <bool check = true>
  BoutReal& getAt(Field2D& f, int off) const {
    return static_cast<const impl*>(this)->template _getAt<check>(f, off);
  }
  /// get the value at a given offset `off` of a field `f`.
  template <bool check = true>
  BoutReal& getAt(const Field2D& f, int off) const {
    return static_cast<const impl*>(this)->template _getAt<check>(f, off);
  }

  /// Get the first point in the boundary
  const BoutReal& next(const Field2D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 0);
  }
  /// Get the first point in the boundary
  BoutReal& next(Field2D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 0);
  }
  /// Get the last point in the domain
  const BoutReal& current(const Field2D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 1);
  }
  /// Get the last point in the domain
  BoutReal& current(Field2D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 1);
  }
  /// Get the second to last point in the domain - this may not be valid and thus throw
  const BoutReal& prev(const Field2D& f) const {
    return static_cast<const impl*>(this)->_getAt(f, 2);
  }

  /*
   *         FUNCTIONS ACCESSORS
   */

  /// get the value at a given offset `off` of a field `f`.
  /// off = -1 is the second point in the boundary
  /// off = 0 is the first point in the boundary
  /// off = 1 is the last point in the domain
  /// off = 2 is the second to last point in the domain
  template <bool check = true>
  BoutReal& getAt(const std::function<BoutReal(int yoffset, Ind3D ind)>& func,
                  int off) const {
    return static_cast<const impl*>(this)->template _getAt<check>(func, off);
  }
  /// Get the first point in the boundary
  const BoutReal&
  next(const std::function<BoutReal(int yoffset, Ind3D ind)>& func) const {
    return static_cast<const impl*>(this)->_getAt(func, 0);
  }
  /// Get the last point in the domain
  const BoutReal&
  current(const std::function<BoutReal(int yoffset, Ind3D ind)>& func) const {
    return static_cast<const impl*>(this)->_getAt(func, 1);
  }
  /// Get the second to last point in the domain - this may not be valid and thus throw
  const BoutReal&
  prev(const std::function<BoutReal(int yoffset, Ind3D ind)>& func) const {
    return static_cast<const impl*>(this)->_getAt(func, 2);
  }

  /*
   *     INTERPOLATION and EXTRAPOLATION
   */

  // extrapolate a given field to the boundary
  BoutReal extrapolate_boundary_o1(const Field3D& f) const { return current(f); }
  // extrapolate a given field to the boundary
  BoutReal extrapolate_boundary_o2(const Field3D& f) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_boundary_o1(f);
    }
    return current(f) * (1 + length()) - prev(f) * length();
  }
  /// Extrapolate a given function to the boundary
  BoutReal extrapolate_bounday_o1(
      const std::function<BoutReal(int yoffset, Ind3D ind)>& func) const {
    return current(func);
  }
  /// Extrapolate a given function to the boundary
  BoutReal extrapolate_boundary_o2(
      const std::function<BoutReal(int yoffset, Ind3D ind)>& func) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_boundary_o1(func);
    }
    return current(func) * (1 + length()) - prev(func) * length();
  }

  /// Interpolate a field to the boundary, using the boundary values
  BoutReal interpolate_boundary_o2(const Field3D& f) const {
    return current(f) * (1 - length()) + next(f) * length();
  }
  /// Interpolate a field to the boundary, using the boundary values
  BoutReal interpolate_boundary_o2(
      const std::function<BoutReal(int yoffset, Ind3D ind)>& func) const {
    return current(func) * (1 - length()) + next(func) * length();
  }
  /// Extrapolate to the first boundary value freely
  BoutReal extrapolate_next_o1(const Field3D& f) const { return current(f); }
  /// Extrapolate to the first boundary value freely
  BoutReal extrapolate_next_o2(const Field3D& f) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_next_o1(f);
    }
    return current(f) * 2 - prev(f);
  }

  /// Extrapolate to the first boundary value freely
  BoutReal
  extrapolate_next_o1(const std::function<BoutReal(int yoffset, Ind3D ind)>& func) const {
    return current(func);
  }
  /// Extrapolate to the first boundary value freely
  BoutReal
  extrapolate_next_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& func) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_next_o1(func);
    }
    return current(func) * 2 - prev(func);
  }

  /// extrapolate the gradient into the boundary
  BoutReal extrapolate_grad_o1([[maybe_unused]] const Field3D& f) const { return 0; }
  /// extrapolate the gradient into the boundary
  BoutReal extrapolate_grad_o2(const Field3D& f) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_grad_o1(f);
    }
    return current(f) - next(f);
  }

  BoutReal extrapolate_boundary_free(const Field3D& f,
                                     BoundaryFreeExtrapolation mode) const {
    const auto fac = valid() > 0 ? limitFreeScale(prev(f), current(f), mode)
                                 : (mode == BoundaryFreeExtrapolation::linear ? 0 : 1);
    auto val = current(f);
    BoutReal next = mode == BoundaryFreeExtrapolation::linear ? val + fac : val * fac;
    return val * length() + next * (1 - length());
  }

  /*
   *     APPLY BOUNDARY CONDITIONS
   */

  /// Apply a dirichlet boundary condition
  void dirichlet_o1(Field3D& f, BoutReal value) const {
    for (int i = 0; i < boundary_width(); ++i) {
      getAt(f, i) = value;
    }
  }

  /// Apply a dirichlet boundary condition
  void dirichlet_o2(Field3D& f, BoutReal value) const {
    if (length() < small_value) {
      return dirichlet_o1(f, value);
    }
    for (int i = 0; i < boundary_width(); ++i) {
      getAt(f, i) =
          parallel_stencil::dirichlet_o2(i + 1, current(f), i + 1 - length(), value);
    }
  }

  /// Apply a dirichlet boundary condition
  void dirichlet_o3(Field3D& f, BoutReal value) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return dirichlet_o2(f, value);
    }
    if (length() < small_value) {
      for (int i = 0; i < boundary_width(); ++i) {
        getAt(f, i) =
            parallel_stencil::dirichlet_o2(i + 2, prev(f), i + 1 - length(), value);
      }
    } else {
      for (int i = 0; i < boundary_width(); ++i) {
        getAt(f, i) = parallel_stencil::dirichlet_o3(i + 2, prev(f), i + 1, current(f),
                                                     i + 1 - length(), value);
      }
    }
  }

  /// Ensure the value in the boundary is at least `value`
  void limit_at_least(Field3D& f, BoutReal value) const {
    for (int i = 0; i < boundary_width(); ++i) {
      if (getAt(f, i) < value) {
        getAt(f, i) = value;
      }
    }
  }

  /// Apply neumann boundary condition, where `value` is the gradient in index space

  // neumann_o1 would give second order convergence, given an appropriate one-sided stencil.
  // But in general we do not, and thus for normal C2 stencils, this is 1st order.
  void neumann_o1(Field3D& f, BoutReal value) const {
    for (int i = 0; i < boundary_width(); ++i) {
      getAt(f, i) = current(f) + value * (i + 1);
    }
  }

  /// Apply neumann boundary condition, where `value` is the gradient in index space
  void neumann_o2(Field3D& f, BoutReal value) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return neumann_o1(f, value);
    }
    for (int i = 0; i < boundary_width(); ++i) {
      getAt(f, i) = prev(f) + (2 + i) * value;
    }
  }

  /// Apply neumann boundary condition, where `value` is the gradient in index space
  void neumann_o3(Field3D& f, BoutReal value) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return neumann_o2(f, value);
    }
    for (int i = 0; i < boundary_width(); ++i) {
      getAt(f, i) = parallel_stencil::neumann_o3(i + 1 - length(), value, i + 1,
                                                 current(f), 2, prev(f));
    }
  }

  void set_free(Field3D& f, BoundaryFreeExtrapolation mode) const {
    int fac;
    if (valid() > 0) {
      fac = limitFreeScale(prev(f), current(f), mode);
    } else {
      fac = mode == BoundaryFreeExtrapolation::linear ? 0 : 1;
    }
    auto val = current(f);
    if (mode == BoundaryFreeExtrapolation::linear) {
      for (int i = 0; i < boundary_width(); ++i) {
        val += fac;
        getAt(f, i) = val;
      }
    } else {
      for (int i = 0; i < boundary_width(); ++i) {
        val *= fac;
        getAt(f, i) = val;
      }
    }
  }
  void setSmallValue(BoutReal val) {
    ASSERT2(val > 0);
    ASSERT2(val < 0.5);
    small_value = val;
  }

private:
  BoutReal small_value = 1e-4;
};

namespace {
/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
///
/// exp( 2*log(fc) - log(fm) )
inline BoutReal limitFreeScale(BoutReal fm, BoutReal fc, BoundaryFreeExtrapolation mode) {
  if ((fm < fc) && (mode == BoundaryFreeExtrapolation::limited)) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }

  BoutReal fp = 0;
  switch (mode) {
  case BoundaryFreeExtrapolation::limited:
  case BoundaryFreeExtrapolation::exponential:
    fp = SQ(fc) / fm; // Exponential
    break;
  case BoundaryFreeExtrapolation::linear:
    fp = (2.0 * fc) - fm; // Linear
    break;
  }

#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree {}: {}, {} -> {}",
                        static_cast<int>(mode), fm, fc, fp);
  }
#endif

  return fp;
}
} // namespace

class BoundaryRegionFCI {
public:
  BoundaryRegionFCI(int dir, Mesh* mesh) : _dir(dir), localmesh(mesh) {};
  /// Add a point to the boundary
  void add_point(Ind3D ind, BoutReal x, BoutReal y, BoutReal z, BoutReal length,
                 char valid, signed char offset) {
    if (!bndry_points.empty() && bndry_points.back().index > ind) {
      is_sorted = false;
    }
    bndry_points.emplace_back(ind, bout::parallel_boundary_region::RealPoint{x, y, z},
                              length, valid, offset,
                              static_cast<unsigned char>(std::abs(offset)));
  }
  void add_point(int ix, int iy, int iz, BoutReal x, BoutReal y, BoutReal z,
                 BoutReal length, char valid, signed char offset) {
    add_point(xyz2ind(ix, iy, iz), x, y, z, length, valid, offset);
  }

private:
  friend class BoundaryRegionIterFCI;
  int _dir;
  // Vector of points in the boundary
  bout::parallel_boundary_region::IndicesVec bndry_points;
  Ind3D xyz2ind(int x, int y, int z) const {
    const int ny = localmesh->LocalNy;
    const int nz = localmesh->LocalNz;
    return Ind3D{(x * ny + y) * nz + z, ny, nz};
  }
  bool is_sorted{true};
  Mesh* localmesh;
};

class BoundaryRegionIterFCI : public BoundaryRegionIterBase<BoundaryRegionIterFCI> {
private:
  BoundaryRegionFCI* region;
  size_t pos{0};

public:
  template <bool check = true>
  BoutReal& _getAt(Field3D& f, int off) const {
    ASSERT3(f.hasParallelSlices());
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() + off * region->_dir;
    return f.ynext(_off)[ind().yp(_off)];
  }
  template <bool check = true>
  const BoutReal& _getAt(const Field3D& f, int off) const {
    ASSERT3(f.hasParallelSlices());
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() + off * region->_dir;
    return f.ynext(_off)[ind().yp(_off)];
  }
  template <bool check = true>
  BoutReal& _getAt(Field2D& f, int off) const {
    ASSERT3(f.hasParallelSlices());
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() + off * region->_dir;
    return f.ynext(_off)[ind().yp(_off)];
  }
  template <bool check = true>
  const BoutReal& _getAt(const Field2D& f, int off) const {
    ASSERT3(f.hasParallelSlices());
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() + off * region->_dir;
    return f.ynext(_off)[ind().yp(_off)];
  }
  template <bool check = true>
  BoutReal getAt(const std::function<BoutReal(int yoffset, Ind3D ind)>& f,
                 int off) const {
    if constexpr (check) {
      ASSERT3(valid() > -off - 2);
    }
    auto _off = _offset() + off * region->_dir;
    return f(_off, ind().yp(_off));
  }
  signed char _offset() const { return region->bndry_points[pos].offset; }
  signed char _valid() const { return region->bndry_points[pos].valid; }
  Ind3D& _ind() const { return region->bndry_points[pos].index; }
  signed char _boundary_width() const {
    return region->localmesh->ystart - region->bndry_points[pos].abs_offset;
  }
  const BoutReal& _length() const { return region->bndry_points[pos].length; }
};
} // namespace boundary
} // namespace bout
