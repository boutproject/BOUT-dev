
class BoundaryRegion;

#ifndef BOUT_BNDRY_REGION_H
#define BOUT_BNDRY_REGION_H

#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include "bout/sys/parallel_stencils.hxx"
#include <string>
#include <utility>

class Mesh;
namespace bout {
namespace globals {
extern Mesh* mesh; ///< Global mesh
} // namespace globals
} // namespace bout

/// Location of boundary
enum class BndryLoc {
  xin,
  xout,
  ydown,
  yup,
  all,
  par_fwd_xin, // Don't include parallel boundaries
  par_bkwd_xin,
  par_fwd_xout, // Don't include parallel boundaries
  par_bkwd_xout
};
constexpr BndryLoc BNDRY_XIN = BndryLoc::xin;
constexpr BndryLoc BNDRY_XOUT = BndryLoc::xout;
constexpr BndryLoc BNDRY_YDOWN = BndryLoc::ydown;
constexpr BndryLoc BNDRY_YUP = BndryLoc::yup;
constexpr BndryLoc BNDRY_ALL = BndryLoc::all;
constexpr BndryLoc BNDRY_PAR_FWD_XIN = BndryLoc::par_fwd_xin;
constexpr BndryLoc BNDRY_PAR_BKWD_XIN = BndryLoc::par_bkwd_xin;
constexpr BndryLoc BNDRY_PAR_FWD_XOUT = BndryLoc::par_fwd_xout;
constexpr BndryLoc BNDRY_PAR_BKWD_XOUT = BndryLoc::par_bkwd_xout;

class BoundaryRegionBase {
public:
  BoundaryRegionBase() = delete;
  BoundaryRegionBase(std::string name, Mesh* passmesh = nullptr)
      : localmesh(passmesh ? passmesh : bout::globals::mesh), label(std::move(name)) {}
  BoundaryRegionBase(std::string name, BndryLoc loc, Mesh* passmesh = nullptr)
      : localmesh(passmesh ? passmesh : bout::globals::mesh), label(std::move(name)),
        location(loc) {}

  virtual ~BoundaryRegionBase() = default;

  Mesh* localmesh; ///< Mesh does this boundary region belongs to

  std::string label; ///< Label for this boundary region

  BndryLoc location;       ///< Which side of the domain is it on?
  bool isParallel = false; ///< Is this a parallel boundary?

  virtual void first() = 0; ///< Move the region iterator to the start
  virtual void next() = 0;  ///< Get the next element in the loop
                            ///  over every element from inside out (in
                            ///  X or Y first)
  virtual bool
  isDone() = 0; ///< Returns true if outside domain. Can use this with nested nextX, nextY
};

class BoundaryRegionIter;
/// Describes a region of the boundary, and a means of iterating over it
class BoundaryRegion : public BoundaryRegionBase {
public:
  BoundaryRegion() = delete;
  BoundaryRegion(std::string name, BndryLoc loc, Mesh* passmesh = nullptr)
      : BoundaryRegionBase(name, loc, passmesh) {}
  BoundaryRegion(std::string name, int xd, int yd, Mesh* passmesh = nullptr)
      : BoundaryRegionBase(name, passmesh), bx(xd), by(yd), width(2) {}
  ~BoundaryRegion() override = default;

  int x, y;   ///< Indices of the point in the boundary
  int bx, by; ///< Direction of the boundary [x+bx][y+by] is going outwards

  int width; ///< Width of the boundary

  virtual void next1d() = 0; ///< Loop over the innermost elements
  virtual void nextX() = 0;  ///< Just loop over X
  virtual void nextY() = 0;  ///< Just loop over Y

  BoundaryRegionIter begin();
  BoundaryRegionIter end();
};

class BoundaryRegionIter {
public:
  BoundaryRegionIter(BoundaryRegion* rgn, bool is_end)
      : rgn(rgn), is_end(is_end), dir(rgn->bx + rgn->by) {
    //static_assert(std::is_base_of<BoundaryRegion, T>, "BoundaryRegionIter only works on BoundaryRegion");

    // Ensure only one is non-zero
    ASSERT3(rgn->bx * rgn->by == 0);
    if (!is_end) {
      rgn->first();
    }
  }
  bool operator!=(const BoundaryRegionIter& rhs) {
    if (is_end) {
      if (rhs.is_end || rhs.rgn->isDone()) {
        return false;
      } else {
        return true;
      }
    }
    if (rhs.is_end) {
      return !rgn->isDone();
    }
    return ind() != rhs.ind();
  }

  Ind3D ind() const { return xyz2ind(rgn->x - rgn->bx, rgn->y - rgn->by, z); }
  BoundaryRegionIter& operator++() {
    ASSERT3(z < nz());
    z++;
    if (z == nz()) {
      z = 0;
      rgn->next();
    }
    return *this;
  }
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
    return 2 * f(0, ind()) - f(0, ind().yp(-rgn->by).xp(-rgn->bx));
  }

  BoutReal interpolate_sheath(const Field3D& f) const {
    return (f[ind()] + ynext(f)) * 0.5;
  }

  BoutReal& ynext(Field3D& f) const { return f[ind().yp(rgn->by).xp(rgn->bx)]; }
  const BoutReal& ynext(const Field3D& f) const {
    return f[ind().yp(rgn->by).xp(rgn->bx)];
  }
  BoutReal& yprev(Field3D& f) const { return f[ind().yp(-rgn->by).xp(-rgn->bx)]; }
  const BoutReal& yprev(const Field3D& f) const {
    return f[ind().yp(-rgn->by).xp(-rgn->bx)];
  }

private:
  BoundaryRegion* rgn;
  const bool is_end;
  int z{0};

public:
  const int dir;

private:
  int nx() const { return rgn->localmesh->LocalNx; }
  int ny() const { return rgn->localmesh->LocalNy; }
  int nz() const { return rgn->localmesh->LocalNz; }

  Ind3D xyz2ind(int x, int y, int z) const {
    return Ind3D{(x * ny() + y) * nz() + z, ny(), nz()};
  }
};

class BoundaryRegionXIn : public BoundaryRegion {
public:
  BoundaryRegionXIn(std::string name, int ymin, int ymax, Mesh* passmesh = nullptr);

  void first() override;
  void next() override;
  void next1d() override;
  void nextX() override;
  void nextY() override;
  bool isDone() override;

private:
  int ys, ye;
};

class BoundaryRegionXOut : public BoundaryRegion {
public:
  BoundaryRegionXOut(std::string name, int ymin, int ymax, Mesh* passmesh = nullptr);

  void first() override;
  void next() override;
  void next1d() override;
  void nextX() override;
  void nextY() override;
  bool isDone() override;

private:
  int ys, ye;
};

class BoundaryRegionYDown : public BoundaryRegion {
public:
  BoundaryRegionYDown(std::string name, int xmin, int xmax, Mesh* passmesh = nullptr);

  void first() override;
  void next() override;
  void next1d() override;
  void nextX() override;
  void nextY() override;
  bool isDone() override;

private:
  int xs, xe;
};

class BoundaryRegionYUp : public BoundaryRegion {
public:
  BoundaryRegionYUp(std::string name, int xmin, int xmax, Mesh* passmesh = nullptr);

  void first() override;
  void next() override;
  void next1d() override;
  void nextX() override;
  void nextY() override;
  bool isDone() override;

private:
  int xs, xe;
};

#endif // BOUT_BNDRY_REGION_H
