#pragma once

#include "bout/mesh.hxx"
#include "bout/sys/parallel_stencils.hxx"
#include "bout/sys/range.hxx"

class BoundaryRegionIter {
public:
  BoundaryRegionIter(int x, int y, int bx, int by, Mesh* mesh)
      : dir(bx + by), x(x), y(y), bx(bx), by(by), localmesh(mesh) {
    ASSERT3(bx * by == 0);
  }
  bool operator!=(const BoundaryRegionIter& rhs) { return ind() != rhs.ind(); }

  Ind3D ind() const { return xyz2ind(x, y, z); }
  BoundaryRegionIter& operator++() {
    ASSERT3(z < nz());
    z++;
    if (z == nz()) {
      z = 0;
      _next();
    }
    return *this;
  }
  virtual void _next() = 0;
  BoundaryRegionIter& operator*() { return *this; }

  void dirichlet_o2(Field3D& f, BoutReal value) const {
    ynext(f) = parallel_stencil::dirichlet_o2(1, f[ind()], 0.5, value);
  }

  BoutReal extrapolate_grad_o2(const Field3D& f) const { return f[ind()] - yprev(f); }

  BoutReal extrapolate_sheath_o2(const Field3D& f) const {
    return (f[ind()] * 3 - yprev(f)) * 0.5;
  }

  BoutReal extrapolate_next_o2(const Field3D& f) const { return 2 * f[ind()] - yprev(f); }

  BoutReal
  extrapolate_next_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& f) const {
    return 2 * f(0, ind()) - f(0, ind().yp(-by).xp(-bx));
  }

  BoutReal interpolate_sheath(const Field3D& f) const {
    return (f[ind()] + ynext(f)) * 0.5;
  }

  BoutReal& ynext(Field3D& f) const { return f[ind().yp(by).xp(bx)]; }
  const BoutReal& ynext(const Field3D& f) const { return f[ind().yp(by).xp(bx)]; }
  BoutReal& yprev(Field3D& f) const { return f[ind().yp(-by).xp(-bx)]; }
  const BoutReal& yprev(const Field3D& f) const { return f[ind().yp(-by).xp(-bx)]; }

  const int dir;

protected:
  int z{0};
  int x;
  int y;
  const int bx;
  const int by;

private:
  Mesh* localmesh;
  int nx() const { return localmesh->LocalNx; }
  int ny() const { return localmesh->LocalNy; }
  int nz() const { return localmesh->LocalNz; }

  Ind3D xyz2ind(int x, int y, int z) const {
    return Ind3D{(x * ny() + y) * nz() + z, ny(), nz()};
  }
};

class BoundaryRegionIterY : public BoundaryRegionIter {
public:
  BoundaryRegionIterY(RangeIterator r, int y, int dir, bool is_end, Mesh* mesh)
      : BoundaryRegionIter(r.ind, y, 0, dir, mesh), r(r), is_end(is_end) {}

  bool operator!=(const BoundaryRegionIterY& rhs) {
    ASSERT2(y == rhs.y);
    if (is_end) {
      if (rhs.is_end) {
        return false;
      }
      return !rhs.r.isDone();
    }
    if (rhs.is_end) {
      return !r.isDone();
    }
    return x != rhs.x;
  }

  virtual void _next() override {
    ++r;
    x = r.ind;
  }

private:
  RangeIterator r;
  bool is_end;
};

class NewBoundaryRegionY {
public:
  NewBoundaryRegionY(Mesh* mesh, bool lower, RangeIterator r)
      : mesh(mesh), lower(lower), r(std::move(r)) {}
  BoundaryRegionIterY begin(bool begin = true) {
    return BoundaryRegionIterY(r, lower ? mesh->ystart : mesh->yend, lower ? -1 : +1,
                               !begin, mesh);
  }
  BoundaryRegionIterY end() { return begin(false); }

private:
  Mesh* mesh;
  bool lower;
  RangeIterator r;
};
