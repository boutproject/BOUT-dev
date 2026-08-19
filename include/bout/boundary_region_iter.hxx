#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include <bout/assert.hxx>
#include <bout/boundary_common.hxx>
#include <bout/boundary_region.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/field_data.hxx>
#include <bout/parallel_boundary_region.hxx>
#include <bout/region.hxx>
#include <bout/sys/parallel_stencils.hxx>
#include <bout/utils.hxx>

namespace bout {
namespace boundary {

/// Helper concept for `BoundaryRegionIterBase::getAt` function accessor overloads
///
/// This is a callable that takes two arguments:
/// - `int yoffset`, the parallel slice offset
/// - `Ind3D ind`, the index of the boundary point
///
/// and returns a `BoutReal`
template <class F>
concept function_accessor =
    std::regular_invocable<F, int, Ind3D>
    and std::is_same_v<std::invoke_result_t<F, int, Ind3D>, BoutReal>;

/// Interface for boundary iterators
template <class Iter>
concept boundary_iterator = requires(Iter point, Field3D f) {
  point.getAt(f, int{});
  point.next(f);
  point.current(f);
  point.prev(f);

  point.ind();
  point.length(CELL_LOC{});
  point.valid();
  point.boundary_width();
  point.is_lower();
  point.offset();

  point.smallValue();
};

/// Common base class for boundary region iterators
///
/// This uses CRTP: boundary region iterators should inherit from this,
/// templated on themselves, and they must implement all methods in
/// `boundary_iterator`, (that is, those methods that call `impl().` in this
/// class)
template <typename Impl>
class BoundaryRegionIterBase {
  BoundaryRegionIterBase() = default;

  /// Get a reference to the derived/implementation class
  const Impl& impl() const { return *static_cast<const Impl*>(this); }

public:
  /// Get the index at the last point in domain
  Ind3D ind() const { return impl()._ind(); }
  /// get the length from the point in the domain to the boundary in index
  /// space. It is in the range [0, 1]
  BoutReal length(CELL_LOC loc) const { return impl()._length(loc); }
  /// Lower bound of how many points are between the first point in the domain
  /// and the boundary in the other direction.
  signed char valid() const { return impl()._valid(); }
  /// Get the width of the boundary at the current point
  int boundary_width() const { return impl()._boundary_width(); }
  /// Is this the lower boundary?
  bool is_lower() const { return impl()._is_lower(); }
  /// Get the offset from the last point in the domain
  /// For FA this is always ±1, for FCI this can be up to ±MYG, excluding 0
  int offset() const { return impl()._offset(); }

  /*
   *         FIELD3D ACCESSORS
   */

  /// Get the value at a given \p offset of a field \p f.
  ///
  /// `offset = -1` is the second point in the boundary
  /// `offset = 0` is the first point in the boundary
  /// `offset = 1` is the last point in the domain
  /// `offset = 2` is the second to last point in the domain
  ///
  /// |---|---|---|--> interior points
  /// -1  0 ^ 1   2
  ///       |
  ///    boundary
  template <bool check = true>
  BoutReal& getAt(Field3D& f, int offset) const {
    return impl().template _getAt<check>(f, offset);
  }
  /// get the value at a given offset `off` of a field `f`.
  template <bool check = true>
  BoutReal getAt(const Field3D& f, int offset) const {
    return impl().template _getAt<check>(f, offset);
  }

  /// Get the first point in the boundary
  const BoutReal& next(const Field3D& f) const { return impl()._getAt(f, 0); }
  /// Get the first point in the boundary
  BoutReal& next(Field3D& f) const { return impl()._getAt(f, 0); }
  /// Get the last point in the domain
  const BoutReal& current(const Field3D& f) const { return impl()._getAt(f, 1); }
  /// Get the last point in the domain
  BoutReal& current(Field3D& f) const { return impl()._getAt(f, 1); }
  /// Get the second to last point in the domain - this may not be valid and thus throw
  const BoutReal& prev(const Field3D& f) const { return impl()._getAt(f, 2); }

  /*
   *         FIELD2D ACCESSORS
   */

  /// Get the value at a given \p offset of a field \p f.
  ///
  /// `offset = -1` is the second point in the boundary
  /// `offset = 0` is the first point in the boundary
  /// `offset = 1` is the last point in the domain
  /// `offset = 2` is the second to last point in the domain
  ///
  /// |---|---|---|--> interior points
  /// -1  0 ^ 1   2
  ///       |
  ///    boundary
  template <bool check = true>
  BoutReal& getAt(Field2D& f, int offset) const {
    return impl().template _getAt<check>(f, offset);
  }
  /// get the value at a given offset `off` of a field `f`.
  template <bool check = true>
  BoutReal getAt(const Field2D& f, int offset) const {
    return impl().template _getAt<check>(f, offset);
  }

  /// Get the first point in the boundary
  const BoutReal& next(const Field2D& f) const { return impl()._getAt(f, 0); }
  /// Get the first point in the boundary
  BoutReal& next(Field2D& f) const { return impl()._getAt(f, 0); }
  /// Get the last point in the domain
  const BoutReal& current(const Field2D& f) const { return impl()._getAt(f, 1); }
  /// Get the last point in the domain
  BoutReal& current(Field2D& f) const { return impl()._getAt(f, 1); }
  /// Get the second to last point in the domain - this may not be valid and thus throw
  const BoutReal& prev(const Field2D& f) const { return impl()._getAt(f, 2); }

  /*
   *         FUNCTIONS ACCESSORS
   */

  /// Apply the function \p func at a given \p offset
  ///
  /// |---|---|---|--> interior points
  /// -1  0 ^ 1   2
  ///       |
  ///    boundary
  template <bool check = true>
  BoutReal getAt(const function_accessor auto& func, int offset) const {
    return impl().template _getAt<check>(func, offset);
  }
  /// Get the first point in the boundary
  BoutReal next(const function_accessor auto& func) const {
    return impl()._getAt(func, 0);
  }
  /// Get the last point in the domain
  BoutReal current(const function_accessor auto& func) const {
    return impl()._getAt(func, 1);
  }
  /// Get the second to last point in the domain - this may not be valid and thus throw
  BoutReal prev(const function_accessor auto& func) const {
    return impl()._getAt(func, 2);
  }

  void setSmallValue(BoutReal val) {
    ASSERT2(val > 0);
    ASSERT2(val < 0.5);
    small_value = val;
  }

  BoutReal smallValue() const { return small_value; }

private:
  BoutReal small_value = 1e-4;
  friend Impl;
};

namespace details {
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
} // namespace details

/// An FCI-aware boundary region
///
/// Uses `BoundaryRegionIterFCI` as its iterator.
///
/// This can't use the legacy iteration methods (`first()`, `next()`, and so on)
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

/// Iterator over a `BoundaryRegionFCI`
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
  bool _is_lower() const { return region->_dir < 0; }

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
  BoutReal _getAt(const function_accessor auto& f, int off) const {
    if constexpr (check) {
      ASSERT3(valid() > -off - 2);
    }
    auto _off = _offset() + (off * region->_dir);
    return f(_off, _ind().yp(_off));
  }
  signed char _offset() const { return region->bndry_points[pos].offset; }
  signed char _valid() const { return region->bndry_points[pos].valid; }
  Ind3D _ind() const { return region->bndry_points[pos].index; }
  int _boundary_width() const {
    return region->localmesh->ystart - region->bndry_points[pos].abs_offset + 1;
  }
  BoutReal _length([[maybe_unused]] CELL_LOC loc) const {
    ASSERT3(loc == CELL_CENTRE);
    return region->bndry_points[pos].length;
  }

  auto operator<=>(const BoundaryRegionIterFCI& rhs) const {
    ASSERT3(region == rhs.region);
    return pos <=> rhs.pos;
  }

  bool operator==(BoundaryRegionIterFCI lhs) const {
    ASSERT3(region == lhs.region);
    return pos == lhs.pos;
  }

  BoundaryRegionIterFCI& operator++() {
    ++pos;
    return *this;
  }
  // No-op for compatibility
  BoundaryRegionIterFCI& operator*() { return *this; }
};

/// Boundary region for field-aligned grids
///
/// FCI grids should use `BoundaryRegionFCI`
///
/// Uses `BoundaryRegionIterXY` as its iterator
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
      location = dir == 1 ? BNDRY_YUP : BNDRY_YDOWN;
    }
  }
  int dir() const { return _dir; }
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

/// Iterator over a `BoundaryRegionXY`
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
  Ind3D ind() const { return _ind(); }

  template <bool check = true>
  BoutReal& _getAt(Field3D& f, int off) const {
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = (1 - off) * region->_dir;
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
    auto _off = (1 - off) * region->_dir;
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
    auto _off = (1 - off) * region->_dir;
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
    auto _off = (1 - off) * region->_dir;
    if constexpr (isX) {
      return f[_ind().xp(_off)];
    } else {
      return f[_ind().yp(_off)];
    }
  }
  template <bool check = true>
  BoutReal _getAt(const function_accessor auto& f, int off) const {
    if constexpr (check) {
      ASSERT3(_valid() > -off - 2);
    }
    auto _off = (1 - off) * region->_dir;
    if constexpr (isX) {
      return f(0, _ind().xp(_off));
    } else {
      return f(0, _ind().yp(_off));
    }
  }
  signed char _offset() const { return region->_dir; }
  signed char _valid() const { return region->valid; }
  Ind3D _ind() const { return region->rgn[pos]; }
  int _boundary_width() const {
    if constexpr (isX) {
      return region->localmesh->xstart;
    }
    return region->localmesh->ystart;
  }
  BoutReal _length(CELL_LOC loc) const {
    if (loc == CELL_XLOW) {
      if (dir() == 1) {
        return 1;
      }
      return 0;
    }
    return 0.5;
  }

  auto operator<=>(const BoundaryRegionIterXY<isX>& rhs) const {
    ASSERT3(region == rhs.region);
    return pos <=> rhs.pos;
  }

  bool operator==(BoundaryRegionIterXY<isX> rhs) const {
    ASSERT3(region == rhs.region);
    return pos == rhs.pos;
  }

  BoundaryRegionIterXY& operator++() {
    ++pos;
    return *this;
  }
  // No-op for compatibility
  BoundaryRegionIterXY& operator*() { return *this; }
};

/// Alias for boundary regions over X specifically
using BoundaryRegionX = BoundaryRegionXY<true>;
/// Alias for boundary regions over Y specifically
using BoundaryRegionY = BoundaryRegionXY<false>;
/// Alias for boundary region iterators over X specifically
using BoundaryRegionIterX = BoundaryRegionIterXY<true>;
/// Alias for boundary region iterators over Y specifically
using BoundaryRegionIterY = BoundaryRegionIterXY<false>;

inline std::shared_ptr<BoundaryRegionX>
NewBoundaryRegionXIn(const std::string& name, int ymin, int ymax, Mesh* mesh) {
  auto pointer = std::make_shared<BoundaryRegionX>(
      name, -1, mesh,
      Region<Ind3D>(mesh->xstart, mesh->xstart, ymin, ymax, mesh->zstart, mesh->zend,
                    mesh->LocalNy, mesh->LocalNz, mesh->maxregionblocksize));
  pointer->legacy = std::make_unique<::BoundaryRegionXIn>(name, ymin, ymax, mesh);
  return pointer;
}

inline std::shared_ptr<BoundaryRegionX>
NewBoundaryRegionXOut(const std::string& name, int ymin, int ymax, Mesh* mesh) {
  auto pointer = std::make_shared<BoundaryRegionX>(
      name, 1, mesh,
      Region<Ind3D>(mesh->xend, mesh->xend, ymin, ymax, mesh->zstart, mesh->zend,
                    mesh->LocalNy, mesh->LocalNz, mesh->maxregionblocksize));
  pointer->legacy = std::make_unique<::BoundaryRegionXOut>(name, ymin, ymax, mesh);
  return pointer;
}

inline std::shared_ptr<BoundaryRegionY>
NewBoundaryRegionYUp(const std::string& name, int xmin, int xmax, Mesh* mesh) {
  auto pointer = std::make_shared<BoundaryRegionY>(
      name, 1, mesh,
      Region<Ind3D>(xmin, xmax, mesh->yend, mesh->yend, mesh->zstart, mesh->zend,
                    mesh->LocalNy, mesh->LocalNz, mesh->maxregionblocksize));
  pointer->legacy = std::make_unique<::BoundaryRegionYUp>(name, xmin, xmax, mesh);
  return pointer;
}

inline std::shared_ptr<BoundaryRegionY>
NewBoundaryRegionYDown(const std::string& name, int xmin, int xmax, Mesh* mesh) {
  auto pointer = std::make_shared<BoundaryRegionY>(
      name, -1, mesh,
      Region<Ind3D>(xmin, xmax, mesh->ystart, mesh->ystart, mesh->zstart, mesh->zend,
                    mesh->LocalNy, mesh->LocalNz, mesh->maxregionblocksize));
  pointer->legacy = std::make_unique<::BoundaryRegionYDown>(name, xmin, xmax, mesh);
  return pointer;
}

template <class Func>
void iter_boundary(const std::shared_ptr<const BoundaryRegionBase>& bndrybase,
                   const Func& func) {
  if (bndrybase->isX) {
    const auto* const bndry = dynamic_cast<const BoundaryRegionX*>(bndrybase.get());
    return iter_boundary(*bndry, func);
  }
  if (bndrybase->isY) {
    const auto* const bndry = dynamic_cast<const BoundaryRegionY*>(bndrybase.get());
    return iter_boundary(*bndry, func);
  }
  if (bndrybase->isParallel) {
    const auto* const bndry = dynamic_cast<const BoundaryRegionFCI*>(bndrybase.get());
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

inline auto begin(const BoundaryRegionFCI& reg) {
  return BoundaryRegionIterFCI(&reg, true);
}
inline auto end(const BoundaryRegionFCI& reg) {
  return BoundaryRegionIterFCI(&reg, false);
}

template <bool isX>
inline auto begin(const BoundaryRegionXY<isX>& reg) {
  return BoundaryRegionIterXY<isX>(&reg, true);
}
template <bool isX>
inline auto end(const BoundaryRegionXY<isX>& reg) {
  return BoundaryRegionIterXY<isX>(&reg, false);
}

/*
 *     INTERPOLATION and EXTRAPOLATION
 */

/// Extrapolate a given field to the boundary
template <boundary_iterator Iter>
BoutReal extrapolate_boundary_o1(const Iter& point, const Field3D& f) {
  return point.current(f);
}

/// Extrapolate a given field to the boundary
template <boundary_iterator Iter>
BoutReal extrapolate_boundary_o2(const Iter& point, const Field3D& f) {
  ASSERT3(point.valid() >= 0);
  if (point.valid() < 1) {
    return extrapolate_boundary_o1(point, f);
  }
  return point.current(f) * (1 + point.length(f.getLocation()))
         - point.prev(f) * point.length(f.getLocation());
}

/// Extrapolate a given function to the boundary
template <boundary_iterator Iter>
BoutReal extrapolate_bounday_o1(const Iter& point, const function_accessor auto& func,
                                [[maybe_unused]] CELL_LOC loc = CELL_CENTRE) {
  return point.current(func);
}

/// Extrapolate a given function to the boundary
template <boundary_iterator Iter>
BoutReal extrapolate_boundary_o2(const Iter& point, const function_accessor auto& func,
                                 CELL_LOC loc = CELL_CENTRE) {
  ASSERT3(point.valid() >= 0);
  if (point.valid() < 1) {
    return extrapolate_boundary_o1(point, func);
  }
  return point.current(func) * (1 + point.length(loc))
         - point.prev(func) * point.length(loc);
}

/// Interpolate a field to the boundary, using the boundary values
template <boundary_iterator Iter>
BoutReal interpolate_boundary_o2(const Iter& point, const Field3D& f) {
  return point.current(f) * (1 - point.length(f.getLocation()))
         + point.next(f) * point.length(f.getLocation());
}
/// Interpolate a field to the boundary, using the boundary values
template <boundary_iterator Iter>
BoutReal interpolate_boundary_o2(const Iter& point, const function_accessor auto& func,
                                 CELL_LOC loc = CELL_CENTRE) {
  return point.current(func) * (1 - point.length(loc))
         + point.next(func) * point.length(loc);
}
/// Extrapolate to the first boundary value freely
template <boundary_iterator Iter>
BoutReal extrapolate_next_o1(const Iter& point, const Field3D& f) {
  return point.current(f);
}
/// Extrapolate to the first boundary value freely
template <boundary_iterator Iter>
BoutReal extrapolate_next_o2(const Iter& point, const Field3D& f) {
  ASSERT3(point.valid() >= 0);
  if (point.valid() < 1) {
    return extrapolate_next_o1(point, f);
  }
  return point.current(f) * 2 - point.prev(f);
}

/// Extrapolate to the first boundary value freely
template <boundary_iterator Iter>
BoutReal extrapolate_next_o1(const Iter& point, const function_accessor auto& func) {
  return point.current(func);
}
/// Extrapolate to the first boundary value freely
template <boundary_iterator Iter>
BoutReal extrapolate_next_o2(const Iter& point, const function_accessor auto& func) {
  ASSERT3(point.valid() >= 0);
  if (point.valid() < 1) {
    return extrapolate_next_o1(func);
  }
  return point.current(func) * 2 - prev(func);
}

/// extrapolate the gradient into the boundary
template <boundary_iterator Iter>
BoutReal extrapolate_grad_o1([[maybe_unused]] const Iter& point,
                             [[maybe_unused]] const Field3D& f) {
  return 0;
}
/// extrapolate the gradient into the boundary
template <boundary_iterator Iter>
BoutReal extrapolate_grad_o2(const Iter& point, const Field3D& f) {
  ASSERT3(point.valid() >= 0);
  if (point.valid() < 1) {
    return extrapolate_grad_o1(point, f);
  }
  return point.current(f) - point.next(f);
}

template <boundary_iterator Iter>
BoutReal extrapolate_boundary_free(const Iter& point, const Field3D& f,
                                   BoundaryFreeExtrapolation mode) {
  BoutReal fac = BoutNaN;
  if (point.valid() > 0) {
    fac = limitFreeScale(point.prev(f), point.current(f), mode);
  } else {
    fac = mode == BoundaryFreeExtrapolation::linear ? 0 : 1;
  }
  const auto val = point.current(f);
  const BoutReal next = mode == BoundaryFreeExtrapolation::linear ? val + fac : val * fac;
  return val * point.length(f.getLocation()) + next * (1 - point.length(f.getLocation()));
}

/*
 *     APPLY BOUNDARY CONDITIONS
 */

/// Apply a dirichlet boundary condition
template <boundary_iterator Iter>
void dirichlet_o1(const Iter& point, Field3D& f, BoutReal value) {
  for (int i = 0; i < point.boundary_width(); ++i) {
    point.getAt(f, -i) = value;
  }
}

/// Apply a dirichlet boundary condition
template <boundary_iterator Iter>
void dirichlet_o2(const Iter& point, Field3D& f, BoutReal value) {
  if (point.length(f.getLocation()) < point.smallValue()) {
    return dirichlet_o1(point, f, value);
  }
  for (int i = 0; i < point.boundary_width(); ++i) {
    point.getAt(f, -i) = parallel_stencil::dirichlet_o2(
        i + 1, point.current(f), i + 1 - point.length(f.getLocation()), value);
  }
}

/// Apply a dirichlet boundary condition
template <boundary_iterator Iter>
void dirichlet_o3(const Iter& point, Field3D& f, BoutReal value) {
  ASSERT3(point.valid() >= 0);
  if (point.valid() < 1) {
    return dirichlet_o2(point, f, value);
  }
  if (point.length(f.getLocation()) < point.smallValue()) {
    for (int i = 0; i < point.boundary_width(); ++i) {
      point.getAt(f, -i) = parallel_stencil::dirichlet_o2(
          i + 2, point.prev(f), i + 1 - point.length(f.getLocation()), value);
    }
  } else {
    for (int i = 0; i < point.boundary_width(); ++i) {
      point.getAt(f, -i) =
          parallel_stencil::dirichlet_o3(i + 2, point.prev(f), i + 1, point.current(f),
                                         i + 1 - point.length(f.getLocation()), value);
    }
  }
}

/// Ensure the value in the boundary is at least `value`
template <boundary_iterator Iter>
void limit_at_least(const Iter& point, Field3D& f, BoutReal value) {
  for (int i = 0; i < point.boundary_width(); ++i) {
    if (point.getAt(f, -i) < value) {
      point.getAt(f, -i) = value;
    }
  }
}

/// Apply neumann boundary condition, where `value` is the gradient in index space

// neumann_o1 would give second order convergence, given an appropriate one-sided stencil.
// But in general we do not, and thus for normal C2 stencils, this is 1st order.
template <boundary_iterator Iter>
void neumann_o1(const Iter& point, Field3D& f, BoutReal value) {
  for (int i = 0; i < point.boundary_width(); ++i) {
    point.getAt(f, -i) = point.current(f) + value * (i + 1);
  }
}

/// Apply neumann boundary condition, where `value` is the gradient in index space
template <boundary_iterator Iter>
void neumann_o2(const Iter& point, Field3D& f, BoutReal value) {
  ASSERT3(point.valid() >= 0);
  if (point.valid() < 1) {
    return neumann_o1(point, f, value);
  }
  for (int i = 0; i < point.boundary_width(); ++i) {
    point.getAt(f, -i) = point.prev(f) + (2 + i) * value;
  }
}

/// Apply neumann boundary condition, where `value` is the gradient in index space
template <boundary_iterator Iter>
void neumann_o3(const Iter& point, Field3D& f, BoutReal value) {
  ASSERT3(point.valid() >= 0);
  if (point.valid() < 1) {
    return neumann_o2(point, f, value);
  }
  for (int i = 0; i < point.boundary_width(); ++i) {
    point.getAt(f, -i) =
        parallel_stencil::neumann_o3(i + 1 - point.length(f.getLocation()), value, i + 1,
                                     point.current(f), 2, point.prev(f));
  }
}

template <boundary_iterator Iter>
void set_free(const Iter& point, Field3D& f, BoundaryFreeExtrapolation mode) {
  BoutReal fac = BoutNaN;
  if (point.valid() > 0) {
    fac = limitFreeScale(point.prev(f), point.current(f), mode);
  } else {
    fac = mode == BoundaryFreeExtrapolation::linear ? 0 : 1;
  }
  auto val = point.current(f);
  if (mode == BoundaryFreeExtrapolation::linear) {
    for (int i = 0; i < point.boundary_width(); ++i) {
      val += fac;
      point.getAt(f, -i) = val;
    }
  } else {
    for (int i = 0; i < point.boundary_width(); ++i) {
      val *= fac;
      point.getAt(f, -i) = val;
    }
  }
}

} // namespace boundary
} // namespace bout
