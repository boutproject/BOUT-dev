#pragma once

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/parallel_boundary_region.hxx"
#include "bout/region.hxx"
#include "bout/sys/parallel_stencils.hxx"
#include "bout/sys/range.hxx"
#include <algorithm>
#include <functional>

class BoundaryRegionIter {
public:
  virtual ~BoundaryRegionIter() = default;
  BoundaryRegionIter(int x, int y, int bx, int by, const Mesh& mesh)
      : _dir(bx + by), _x(x), _y(y), _bx(bx), _by(by), localmesh(&mesh) {
    ASSERT3(_bx * _by == 0);
  }
  bool operator!=(const BoundaryRegionIter& rhs) const { return ind() != rhs.ind(); }

  Ind3D ind() const { return xyz2ind(_x, _y, _z); }
  BoundaryRegionIter& operator++() {
    ASSERT3(_z < nz());
    _z++;
    if (_z == nz()) {
      _z = 0;
      _next();
    }
    return *this;
  }

protected:
  virtual void _next() = 0;

public:
  BoundaryRegionIter& operator*() { return *this; }

  void dirichlet_o2(Field3D& f, BoutReal value) const {
    ynext(f) = bout::parallel_stencil::dirichlet_o2(1, f[ind()], 0.5, value);
  }

  BoutReal extrapolate_grad_o2(const Field3D& f) const { return f[ind()] - yprev(f); }

  BoutReal extrapolate_sheath_o2(const Field3D& f) const {
    return (f[ind()] * 3 - yprev(f)) * 0.5;
  }

  BoutReal extrapolate_next_o2(const Field3D& f) const {
    return (2 * f[ind()]) - yprev(f);
  }

  BoutReal
  extrapolate_next_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& f) const {
    return (2 * f(0, ind())) - f(0, ind().yp(-_by).xp(-_bx));
  }

  BoutReal interpolate_sheath_o2(const Field3D& f) const {
    return (f[ind()] + ynext(f)) * 0.5;
  }

  BoutReal
  interpolate_sheath_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& f) const {
    return (f(0, ind()) + f(0, ind().yp(-_by).xp(-_bx))) * 0.5;
  }

  BoutReal
  extrapolate_sheath_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& f) const {
    return 0.5 * (3 * f(0, ind()) - f(0, ind().yp(-_by).xp(-_bx)));
  }

  BoutReal extrapolate_sheath_free(const Field3D& f, SheathLimitMode mode) const {
    const BoutReal fac =
        bout::parallel_boundary_region::limitFreeScale(yprev(f), ythis(f), mode);
    const BoutReal val = ythis(f);
    const BoutReal next = mode == SheathLimitMode::linear_free ? val + fac : val * fac;
    return 0.5 * (val + next);
  }

  void set_free(Field3D& f, SheathLimitMode mode) const {
    const BoutReal fac =
        bout::parallel_boundary_region::limitFreeScale(yprev(f), ythis(f), mode);
    BoutReal val = ythis(f);
    if (mode == SheathLimitMode::linear_free) {
      for (int i = 1; i <= localmesh->ystart; ++i) {
        val += fac;
        f[ind().yp(_by * i).xp(_bx * i)] = val;
      }
    } else {
      for (int i = 1; i <= localmesh->ystart; ++i) {
        val *= fac;
        f[ind().yp(_by * i).xp(_bx * i)] = val;
      }
    }
  }

  void limitFree(Field3D& f) const {
    const BoutReal fac =
        bout::parallel_boundary_region::limitFreeScale(yprev(f), ythis(f));
    BoutReal val = ythis(f);
    for (int i = 1; i <= localmesh->ystart; ++i) {
      val *= fac;
      f[ind().yp(_by * i).xp(_bx * i)] = val;
    }
  }

  bool is_lower() const {
    ASSERT2(_bx == 0);
    return _by == -1;
  }

  void neumann_o1(Field3D& f, BoutReal grad) const {
    BoutReal val = ythis(f);
    for (int i = 1; i <= localmesh->ystart; ++i) {
      val += grad;
      f[ind().yp(_by * i).xp(_bx * i)] = val;
    }
  }

  void neumann_o2(Field3D& f, BoutReal grad) const {
    BoutReal val = yprev(f) + grad;
    for (int i = 1; i <= localmesh->ystart; ++i) {
      val += grad;
      f[ind().yp(_by * i).xp(_bx * i)] = val;
    }
  }

  void limit_at_least(Field3D& f, BoutReal value) const {
    ynext(f) = std::max(ynext(f), value);
  }

  BoutReal& ynext(Field3D& f) const { return f[ind().yp(_by).xp(_bx)]; }
  const BoutReal& ynext(const Field3D& f) const { return f[ind().yp(_by).xp(_bx)]; }
  BoutReal& yprev(Field3D& f) const { return f[ind().yp(-_by).xp(-_bx)]; }
  const BoutReal& yprev(const Field3D& f) const { return f[ind().yp(-_by).xp(-_bx)]; }
  BoutReal& ythis(Field3D& f) const { return f[ind()]; }
  const BoutReal& ythis(const Field3D& f) const { return f[ind()]; }

  void setYPrevIfValid(Field3D& f, BoutReal val) const { yprev(f) = val; }
  void setAll(Field3D& f, const BoutReal val) const {
    for (int i = -localmesh->ystart; i <= localmesh->ystart; ++i) {
      f[ind().yp(_by * i).xp(_bx * i)] = val;
    }
  }

  static int abs_offset() { return 1; }

  BoutReal& ynext(Field2D& f) const { return f[ind().yp(_by).xp(_bx)]; }
  const BoutReal& ynext(const Field2D& f) const { return f[ind().yp(_by).xp(_bx)]; }
  BoutReal& yprev(Field2D& f) const { return f[ind().yp(-_by).xp(-_bx)]; }
  const BoutReal& yprev(const Field2D& f) const { return f[ind().yp(-_by).xp(-_bx)]; }

  int dir() const { return _dir; }

protected:
  int x() const { return _x; }
  int y() const { return _y; }
  int z() const { return _z; }

  void setx(int x) { _x = x; }
  void sety(int y) { _y = y; }

private:
  int _dir;

  int _z{0};
  int _x;
  int _y;
  int _bx;
  int _by;

  const Mesh* localmesh;
  int nx() const { return localmesh->LocalNx; }
  int ny() const { return localmesh->LocalNy; }
  int nz() const { return localmesh->LocalNz; }

  Ind3D xyz2ind(int x, int y, int z) const {
    return Ind3D{((x * ny() + y) * nz()) + z, ny(), nz()};
  }
};

class BoundaryRegionIterY : public BoundaryRegionIter {
public:
  ~BoundaryRegionIterY() override = default;
  BoundaryRegionIterY(const RangeIterator& r, int y, int dir, bool is_end,
                      const Mesh& mesh)
      : BoundaryRegionIter(r.ind, y, 0, dir, mesh), r(r), is_end(is_end) {}

  bool operator!=(const BoundaryRegionIterY& rhs) {
    ASSERT2(y() == rhs.y());
    if (is_end) {
      if (rhs.is_end) {
        return false;
      }
      return !rhs.r.isDone();
    }
    if (rhs.is_end) {
      return !r.isDone();
    }
    return x() != rhs.x();
  }

  void _next() override {
    ++r;
    setx(r.ind);
  }

private:
  RangeIterator r;
  bool is_end;
};

class NewBoundaryRegionY {
public:
  NewBoundaryRegionY(const Mesh& mesh, bool lower, const RangeIterator& r)
      : mesh(&mesh), lower(lower), r(r) {}
  BoundaryRegionIterY begin(bool begin = true) {
    return BoundaryRegionIterY(r, lower ? mesh->ystart : mesh->yend, lower ? -1 : +1,
                               !begin, *mesh);
  }
  BoundaryRegionIterY end() { return begin(false); }

private:
  const Mesh* mesh;
  bool lower;
  RangeIterator r;
};
