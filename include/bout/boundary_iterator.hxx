#pragma once

#include "bout/assert.hxx"
#include "bout/boundary_common.hxx"
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
