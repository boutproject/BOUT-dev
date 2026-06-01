#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <type_traits>
#include <vector>

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/field_data.hxx"
#include "bout/utils.hxx"
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
class BoundaryRegionIterBase {
  BoundaryRegionIterBase() = default;
  /// get the index at the last point in domain
public:
  Ind3D ind() const { return static_cast<const impl*>(this)->_ind(); }
  /// get the length from the point in the domain to the boundary in index
  /// space. It is in the range [0, 1]
  BoutReal length(CELL_LOC loc) const {
    return static_cast<const impl*>(this)->_length(loc);
  }
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
  /// Get the offset from the last point in the domain
  /// For FA this is always ±1, for FCI this can be up to ±MYG, excluding 0
  int offset() const { return static_cast<const impl*>(this)->_offset(); }
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
    return current(f) * (1 + length(f.getLocation())) - prev(f) * length(f.getLocation());
  }
  /// Extrapolate a given function to the boundary
  BoutReal
  extrapolate_bounday_o1(const std::function<BoutReal(int yoffset, Ind3D ind)>& func,
                         [[maybe_unused]] CELL_LOC loc = CELL_CENTRE) const {
    return current(func);
  }
  /// Extrapolate a given function to the boundary
  BoutReal
  extrapolate_boundary_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& func,
                          CELL_LOC loc = CELL_CENTRE) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_boundary_o1(func);
    }
    return current(func) * (1 + length(loc)) - prev(func) * length(loc);
  }

  /// Interpolate a field to the boundary, using the boundary values
  BoutReal interpolate_boundary_o2(const Field3D& f) const {
    return current(f) * (1 - length(f.getLocation())) + next(f) * length(f.getLocation());
  }
  /// Interpolate a field to the boundary, using the boundary values
  BoutReal
  interpolate_boundary_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& func,
                          CELL_LOC loc = CELL_CENTRE) const {
    return current(func) * (1 - length(loc)) + next(func) * length(loc);
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
    BoutReal fac = BoutNaN;
    if (valid() > 0) {
      fac = limitFreeScale(prev(f), current(f), mode);
    } else {
      fac = mode == BoundaryFreeExtrapolation::linear ? 0 : 1;
    }
    auto val = current(f);
    BoutReal next = mode == BoundaryFreeExtrapolation::linear ? val + fac : val * fac;
    return val * length(f.getLocation()) + next * (1 - length(f.getLocation()));
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
    if (length(f.getLocation()) < small_value) {
      return dirichlet_o1(f, value);
    }
    for (int i = 0; i < boundary_width(); ++i) {
      getAt(f, i) = parallel_stencil::dirichlet_o2(
          i + 1, current(f), i + 1 - length(f.getLocation()), value);
    }
  }

  /// Apply a dirichlet boundary condition
  void dirichlet_o3(Field3D& f, BoutReal value) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return dirichlet_o2(f, value);
    }
    if (length(f.getLocation()) < small_value) {
      for (int i = 0; i < boundary_width(); ++i) {
        getAt(f, i) = parallel_stencil::dirichlet_o2(
            i + 2, prev(f), i + 1 - length(f.getLocation()), value);
      }
    } else {
      for (int i = 0; i < boundary_width(); ++i) {
        getAt(f, i) = parallel_stencil::dirichlet_o3(
            i + 2, prev(f), i + 1, current(f), i + 1 - length(f.getLocation()), value);
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
      getAt(f, i) = parallel_stencil::neumann_o3(i + 1 - length(f.getLocation()), value,
                                                 i + 1, current(f), 2, prev(f));
    }
  }

  void set_free(Field3D& f, BoundaryFreeExtrapolation mode) const {
    BoutReal fac = BoutNaN;
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
  friend impl;
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

class BoundaryRegionFCI : public BoundaryRegionBase {
public:
  BoundaryRegionFCI(const std::string& name, const BndryLoc& loc, int dir, Mesh* mesh)
      : BoundaryRegionBase(name, loc, mesh), _dir(dir), localmesh(mesh) {
    isParallel = true;
  };
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
  bool contains(int ix, int iy, int iz) {
    const auto ind = xyz2ind(ix, iy, iz);
    ensureSorted();
    const auto found =
        std::lower_bound(std::begin(bndry_points), std::end(bndry_points), ind);
    return found != std::end(bndry_points) and found->index == ind;
  }
  int dir() const { return _dir; }
  // legacy interface
  void first() override { throw BoutException("Legacy interface is not suppored"); }
  void next() override { throw BoutException("Legacy interface is not suppored"); }
  bool isDone() override { throw BoutException("Legacy interface is not suppored"); }

private:
  friend class BoundaryRegionIterFCI;
  int _dir;
  // Vector of points in the boundary
  bout::parallel_boundary_region::IndicesVec bndry_points;
  Mesh* localmesh;
  bool is_sorted{true};
  void ensureSorted() {
    if (is_sorted) {
      return;
    }
    std::sort(std::begin(bndry_points), std::end(bndry_points));
  }
  Ind3D xyz2ind(int x, int y, int z) const {
    const int ny = localmesh->LocalNy;
    const int nz = localmesh->LocalNz;
    return Ind3D{((x * ny + y) * nz) + z, ny, nz};
  }
};

class BoundaryRegionIterFCI : public BoundaryRegionIterBase<BoundaryRegionIterFCI> {
private:
  // TODO(dave) make non-const?
  const BoundaryRegionFCI* region;
  size_t pos{0};

public:
  BoundaryRegionIterFCI() = delete;
  BoundaryRegionIterFCI(const BoundaryRegionFCI* reg, bool isstart)
      : region(reg), pos(isstart ? 0 : reg->bndry_points.size()) {}
  void setValid(char valid) {
    const_cast<BoundaryRegionFCI*>(region)->bndry_points[pos].valid = valid;
  };
  BoutReal s_x() const { return region->bndry_points[pos].intersection.s_x; };
  BoutReal s_y() const { return region->bndry_points[pos].intersection.s_y; };
  BoutReal s_z() const { return region->bndry_points[pos].intersection.s_z; };
  Mesh* localmesh() const { return region->localmesh; };
  int dir() const { return region->_dir; }
  template <bool check = true>
  BoutReal& _getAt(Field3D& f, int off) const {
    ASSERT3(f.hasParallelSlices());
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() - (off * region->_dir);
    return f.ynext(_off)[_ind().yp(_off)];
  }
  template <bool check = true>
  const BoutReal& _getAt(const Field3D& f, int off) const {
    ASSERT3(f.hasParallelSlices());
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() - (off * region->_dir);
    return f.ynext(_off)[_ind().yp(_off)];
  }
  template <bool check = true>
  BoutReal& _getAt(Field2D& f, int off) const {
    ASSERT3(f.hasParallelSlices());
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() - (off * region->_dir);
    return f.ynext(_off)[_ind().yp(_off)];
  }
  template <bool check = true>
  const BoutReal& _getAt(const Field2D& f, int off) const {
    ASSERT3(f.hasParallelSlices());
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() - (off * region->_dir);
    return f.ynext(_off)[_ind().yp(_off)];
  }
  template <bool check = true>
  BoutReal getAt(const std::function<BoutReal(int yoffset, Ind3D ind)>& f,
                 int off) const {
    if constexpr (check) {
      ASSERT3(valid() > -off - 2);
    }
    auto _off = _offset() + (off * region->_dir);
    return f(_off, _ind().yp(_off));
  }
  signed char _offset() const { return region->bndry_points[pos].offset; }
  signed char _valid() const { return region->bndry_points[pos].valid; }
  Ind3D _ind() const { return region->bndry_points[pos].index; }
  signed char _boundary_width() const {
    return region->localmesh->ystart - region->bndry_points[pos].abs_offset;
  }
  BoutReal _length([[maybe_unused]] CELL_LOC loc) const {
    ASSERT3(loc == CELL_CENTRE);
    return region->bndry_points[pos].length;
  }
  bool operator!=(BoundaryRegionIterFCI lhs) const {
    ASSERT3(region == lhs.region);
    return pos != lhs.pos;
  }
  BoundaryRegionIterFCI& operator++() {
    ++pos;
    return *this;
  }
  // No-op for compatibility
  BoundaryRegionIterFCI& operator*() { return *this; }
};

template <bool isXtemp>
class BoundaryRegionXY : public BoundaryRegionBase {
public:
  BoundaryRegionXY() = delete;
  BoundaryRegionXY(const std::string& name, int dir, Mesh* mesh, Region<Ind3D>&& rgn)
      : BoundaryRegionBase(name, mesh), _dir(dir),
        valid(isXtemp ? mesh->xstart : mesh->ystart) {
    BOUT_FOR_SERIAL(i, rgn) { this->rgn.emplace_back(i); }
    if constexpr (isXtemp) {
      this->isX = true;
      location = dir == 1 ? BNDRY_XOUT : BNDRY_XIN;
    } else {
      this->isY = true;
      location = dir == 1 ? BNDRY_YDOWN : BNDRY_YUP;
    }
  }
  int dir() { return _dir; }
  // legacy interface
  void first() override { throw BoutException("Legacy interface is not suppored"); }
  void next() override { throw BoutException("Legacy interface is not suppored"); }
  bool isDone() override { throw BoutException("Legacy interface is not suppored"); }

private:
  template <bool isIterX>
  friend class BoundaryRegionIterXY;
  int _dir;
  std::vector<Ind3D> rgn;
  signed char valid;
};

template <bool isX>
class BoundaryRegionIterXY : public BoundaryRegionIterBase<BoundaryRegionIterXY<isX>> {
private:
  const BoundaryRegionXY<isX>* region;
  size_t pos{0};

public:
  BoundaryRegionIterXY() = delete;
  BoundaryRegionIterXY(const BoundaryRegionXY<isX>* reg, bool isstart)
      : region(reg), pos(isstart ? 0 : reg->rgn.size()) {}
  int dir() const { return region->_dir; }
  template <bool check = true>
  BoutReal& _getAt(Field3D& f, int off) const {
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() - off * region->_dir;
    if constexpr (isX) {
      return f[_ind().xp(_off)];
    } else {
      return f[_ind().yp(_off)];
    }
  }
  template <bool check = true>
  const BoutReal& _getAt(const Field3D& f, int off) const {
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() - off * region->_dir;
    if constexpr (isX) {
      return f[_ind().xp(_off)];
    } else {
      return f[_ind().yp(_off)];
    }
  }
  template <bool check = true>
  BoutReal& _getAt(Field2D& f, int off) const {
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() - off * region->_dir;
    if constexpr (isX) {
      return f[_ind().xp(_off)];
    } else {
      return f[_ind().yp(_off)];
    }
  }
  template <bool check = true>
  const BoutReal& _getAt(const Field2D& f, int off) const {
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() - off * region->_dir;
    if constexpr (isX) {
      return f[_ind().xp(_off)];
    } else {
      return f[_ind().yp(_off)];
    }
  }
  template <bool check = true>
  BoutReal getAt(const std::function<BoutReal(int yoffset, Ind3D ind)>& f,
                 int off) const {
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = _offset() + off * region->_dir;
    if constexpr (isX) {
      return f(0, _ind().xp(_off));
    } else {
      return f(0, _ind().yp(_off));
    }
  }
  signed char _offset() const { return region->_dir; }
  signed char _valid() const { return region->valid; }
  Ind3D _ind() const { return region->rgn[pos]; }
  signed char _boundary_width() const { return region->localmesh->xstart; }
  BoutReal _length(CELL_LOC loc) const {
    if (loc == CELL_XLOW) {
      if (dir() == 1) {
        return 1;
      }
      return 0;
    }
    return 0.5;
  }
  bool operator!=(BoundaryRegionIterXY<isX> lhs) {
    ASSERT3(region == lhs.region);
    return pos != lhs.pos;
  }
  BoundaryRegionIterXY& operator++() {
    ++pos;
    return *this;
  }
  // No-op for compatibility
  BoundaryRegionIterXY& operator*() { return *this; }
};
using BoundaryRegionX = BoundaryRegionXY<true>;
using BoundaryRegionY = BoundaryRegionXY<false>;
using BoundaryRegionIterX = BoundaryRegionIterXY<true>;
using BoundaryRegionIterY = BoundaryRegionIterXY<false>;

inline BoundaryRegionX* NewBoundaryRegionXIn(const std::string& name, int ymin, int ymax,
                                             Mesh* mesh) {
  auto* pointer = new BoundaryRegionX(
      name, -1, mesh,
      Region<Ind3D>(mesh->xstart, mesh->xstart, ymin, ymax, mesh->zstart, mesh->zend,
                    mesh->LocalNy, mesh->LocalNz, mesh->maxregionblocksize));
  pointer->legacy = new ::BoundaryRegionXIn(name, ymin, ymax, mesh);
  return pointer;
}

inline BoundaryRegionX* NewBoundaryRegionXOut(const std::string& name, int ymin, int ymax,
                                              Mesh* mesh) {
  auto* pointer = new BoundaryRegionX(
      name, 1, mesh,
      Region<Ind3D>(mesh->xend, mesh->xend, ymin, ymax, mesh->zstart, mesh->zend,
                    mesh->LocalNy, mesh->LocalNz, mesh->maxregionblocksize));
  pointer->legacy = new ::BoundaryRegionXOut(name, ymin, ymax, mesh);
  return pointer;
}

inline BoundaryRegionY* NewBoundaryRegionYUp(const std::string& name, int xmin, int xmax,
                                             Mesh* mesh) {
  auto* pointer = new BoundaryRegionY(
      name, -1, mesh,
      Region<Ind3D>(xmin, xmax, mesh->yend, mesh->yend, mesh->zstart, mesh->zend,
                    mesh->LocalNy, mesh->LocalNz, mesh->maxregionblocksize));
  pointer->legacy = new ::BoundaryRegionYUp(name, xmin, xmax, mesh);
  return pointer;
}

inline BoundaryRegionY* NewBoundaryRegionYDown(const std::string& name, int xmin,
                                               int xmax, Mesh* mesh) {
  auto* pointer = new BoundaryRegionY(
      name, 1, mesh,
      Region<Ind3D>(xmin, xmax, mesh->ystart, mesh->ystart, mesh->zstart, mesh->zend,
                    mesh->LocalNy, mesh->LocalNz, mesh->maxregionblocksize));
  pointer->legacy = new ::BoundaryRegionYDown(name, xmin, xmax, mesh);
  return pointer;
}

template <class Func>
void iter_boundary(const BoundaryRegionBase* bndrybase, const Func& func) {
  if (bndrybase->isX) {
    const auto* const bndry = dynamic_cast<const BoundaryRegionX*>(bndrybase);
    return iter_boundary(*bndry, func);
  }
  if (bndrybase->isY) {
    const auto* const bndry = dynamic_cast<const BoundaryRegionY*>(bndrybase);
    return iter_boundary(*bndry, func);
  }
  if (bndrybase->isParallel) {
    const auto* const bndry = dynamic_cast<const BoundaryRegionFCI*>(bndrybase);
    return iter_boundary(*bndry, func);
  }
  throw BoutException("{} is of unknown type - probably a legacy iterator",
                      bndrybase->label);
}

template <class Bndry, class Func,
          typename = std::enable_if_t<std::is_base_of<BoundaryRegionBase, Bndry>::value>>
void iter_boundary(const Bndry& bndry, const Func& func) {
  static_assert(std::is_base_of<BoundaryRegionBase, Bndry>::value,
                "Bndry must derive from BoundaryRegionY");
  for (auto& point : bndry) {
    func(point);
  }
}

} // namespace boundary
} // namespace bout
inline bout::boundary::BoundaryRegionIterFCI
begin(const bout::boundary::BoundaryRegionFCI& reg) {
  return bout::boundary::BoundaryRegionIterFCI(&reg, true);
}
inline bout::boundary::BoundaryRegionIterFCI
end(const bout::boundary::BoundaryRegionFCI& reg) {
  return bout::boundary::BoundaryRegionIterFCI(&reg, false);
}

template <bool isX>
inline bout::boundary::BoundaryRegionIterXY<isX>
begin(const bout::boundary::BoundaryRegionXY<isX>& reg) {
  return bout::boundary::BoundaryRegionIterXY<isX>(&reg, true);
}
template <bool isX>
inline bout::boundary::BoundaryRegionIterXY<isX>
end(const bout::boundary::BoundaryRegionXY<isX>& reg) {
  return bout::boundary::BoundaryRegionIterXY<isX>(&reg, false);
}
