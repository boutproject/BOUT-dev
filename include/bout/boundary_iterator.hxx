#pragma once

#include "bout/mesh.hxx"
#include "bout/parallel_boundary_region.hxx"
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

  BoutReal interpolate_sheath_o1(const Field3D& f) const {
    return (f[ind()] + ynext(f)) * 0.5;
  }
  BoutReal
  extrapolate_sheath_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& f) const {
    return 0.5 * (3 * f(0, ind()) - f(0, ind().yp(-by).xp(-bx)));
  }

  void limitFree(Field3D& f) const {
    const BoutReal fac =
        bout::parallel_boundary_region::limitFreeScale(yprev(f), ythis(f));
    BoutReal val = ythis(f);
    for (int i = 1; i <= localmesh->ystart; ++i) {
      val *= fac;
      f[ind().yp(by * i).xp(bx * i)] = val;
    }
  }

  void neumann_o1(Field3D& f, BoutReal grad) const {
    BoutReal val = ythis(f);
    for (int i = 1; i <= localmesh->ystart; ++i) {
      val += grad;
      f[ind().yp(by * i).xp(bx * i)] = val;
    }
  }

  void neumann_o2(Field3D& f, BoutReal grad) const {
    BoutReal val = yprev(f) + grad;
    for (int i = 1; i <= localmesh->ystart; ++i) {
      val += grad;
      f[ind().yp(by * i).xp(bx * i)] = val;
    }
  }

  BoutReal& ynext(Field3D& f) const { return f[ind().yp(by).xp(bx)]; }
  const BoutReal& ynext(const Field3D& f) const { return f[ind().yp(by).xp(bx)]; }
  BoutReal& yprev(Field3D& f) const { return f[ind().yp(-by).xp(-bx)]; }
  const BoutReal& yprev(const Field3D& f) const { return f[ind().yp(-by).xp(-bx)]; }
  BoutReal& ythis(Field3D& f) const { return f[ind()]; }
  const BoutReal& ythis(const Field3D& f) const { return f[ind()]; }

  void setYPrevIfValid(Field3D& f, BoutReal val) const { yprev(f) = val; }
  void setAll(Field3D& f, const BoutReal val) const {
    for (int i = -localmesh->ystart; i <= localmesh->ystart; ++i) {
      f[ind().yp(by * i).xp(bx * i)] = val;
    }
  }

  int abs_offset() const { return 1; }

#if BOUT_USE_METRIC_3D == 0
  BoutReal& ynext(Field2D& f) const { return f[ind().yp(by).xp(bx)]; }
  const BoutReal& ynext(const Field2D& f) const { return f[ind().yp(by).xp(bx)]; }
  BoutReal& yprev(Field2D& f) const { return f[ind().yp(-by).xp(-bx)]; }
  const BoutReal& yprev(const Field2D& f) const { return f[ind().yp(-by).xp(-bx)]; }
#endif

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
