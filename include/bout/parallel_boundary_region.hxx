#ifndef BOUT_PAR_BNDRY_H
#define BOUT_PAR_BNDRY_H

#include "bout/boundary_region.hxx"
#include "bout/bout_types.hxx"
#include <functional>
#include <vector>

#include "bout/sys/parallel_stencils.hxx"
#include <bout/field3d.hxx>
#include <bout/mesh.hxx>

/**
 * Boundary region for parallel direction. This contains a vector of points that are
 * inside the boundary.
 *
 */

namespace bout {
namespace parallel_boundary_region {

struct RealPoint {
  BoutReal s_x;
  BoutReal s_y;
  BoutReal s_z;
};

struct Indices {
  // Indices of the boundary point
  Ind3D index;
  // Intersection with boundary in index space
  RealPoint intersection;
  // Distance to intersection
  BoutReal length;
  // Angle between field line and boundary
  // BoutReal angle;
  // How many points we can go in the opposite direction
  signed char valid;
  signed char offset;
  unsigned char abs_offset;
};

using IndicesVec = std::vector<Indices>;
using IndicesIter = IndicesVec::iterator;
using IndicesIterConst = IndicesVec::const_iterator;

//}

template <class IndicesVec, class IndicesIter>
class BoundaryRegionParIterBase {
public:
  BoundaryRegionParIterBase(IndicesVec& bndry_points, IndicesIter bndry_position, int dir,
                            Mesh* localmesh)
      : bndry_points(bndry_points), bndry_position(bndry_position), dir(dir),
        localmesh(localmesh){};

  // getter
  Ind3D ind() const { return bndry_position->index; }
  BoutReal s_x() const { return bndry_position->intersection.s_x; }
  BoutReal s_y() const { return bndry_position->intersection.s_y; }
  BoutReal s_z() const { return bndry_position->intersection.s_z; }
  BoutReal length() const { return bndry_position->length; }
  signed char valid() const { return bndry_position->valid; }
  signed char offset() const { return bndry_position->offset; }
  unsigned char abs_offset() const { return bndry_position->abs_offset; }


  // extrapolate a given point to the boundary
  BoutReal extrapolate_sheath_o1(const Field3D& f) const { return f[ind()]; }
  BoutReal extrapolate_sheath_o2(const Field3D& f) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_sheath_o1(f);
    }
    return f[ind()] * (1 + length()) - f.ynext(-dir)[ind().yp(-dir)] * length();
  }
  inline BoutReal
  extrapolate_sheath_o1(const std::function<BoutReal(int yoffset, Ind3D ind)>& f) const {
    return f(0, ind());
  }
  inline BoutReal
  extrapolate_sheath_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& f) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_sheath_o1(f);
    }
    return f(0, ind()) * (1 + length()) - f(-dir, ind().yp(-dir)) * length();
  }

  inline BoutReal interpolate_sheath_o1(const Field3D& f) const {
    return f[ind()] * (1 - length()) + ynext(f) * length();
  }

  inline BoutReal extrapolate_next_o1(const Field3D& f) const { return f[ind()]; }
  inline BoutReal extrapolate_next_o2(const Field3D& f) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_next_o1(f);
    }
    return f[ind()] * 2 - f.ynext(-dir)[ind().yp(-dir)];
  }

  inline BoutReal
  extrapolate_next_o1(const std::function<BoutReal(int yoffset, Ind3D ind)>& f) const {
    return f(0, ind());
  }
  inline BoutReal
  extrapolate_next_o2(const std::function<BoutReal(int yoffset, Ind3D ind)>& f) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_sheath_o1(f);
    }
    return f(0, ind()) * 2 - f(-dir, ind().yp(-dir));
  }

  // extrapolate the gradient into the boundary
  inline BoutReal extrapolate_grad_o1(const Field3D& f) const { return 0; }
  inline BoutReal extrapolate_grad_o2(const Field3D& f) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return extrapolate_grad_o1(f);
    }
    return f[ind()] - f.ynext(-dir)[ind().yp(-dir)];
  }

  BoundaryRegionParIterBase& operator*() { return *this; }

  BoundaryRegionParIterBase& operator++() {
    ++bndry_position;
    return *this;
  }

  bool operator!=(const BoundaryRegionParIterBase& rhs) {
    return bndry_position != rhs.bndry_position;
  }

  // dirichlet boundary code
  void dirichlet_o1(Field3D& f, BoutReal value) const {
    f.ynext(dir)[ind().yp(dir)] = value;
  }

  void dirichlet_o2(Field3D& f, BoutReal value) const {
    if (length() < small_value) {
      return dirichlet_o1(f, value);
    }
    ynext(f) = parallel_stencil::dirichlet_o2(1, f[ind()], 1 - length(), value);
    // ynext(f) = f[ind()] * (1 + 1/length()) + value / length();
  }

  void dirichlet_o3(Field3D& f, BoutReal value) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return dirichlet_o2(f, value);
    }
    if (length() < small_value) {
      ynext(f) = parallel_stencil::dirichlet_o2(2, yprev(f), 1 - length(), value);
    } else {
      ynext(f) =
          parallel_stencil::dirichlet_o3(2, yprev(f), 1, f[ind()], 1 - length(), value);
    }
  }

  // NB: value needs to be scaled by dy
  // neumann_o1 is actually o2 if we would use an appropriate one-sided stencil.
  // But in general we do not, and thus for normal C2 stencils, this is 1st order.
  void neumann_o1(Field3D& f, BoutReal value) const { ynext(f) = f[ind()] + value; }

  // NB: value needs to be scaled by dy
  void neumann_o2(Field3D& f, BoutReal value) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return neumann_o1(f, value);
    }
    ynext(f) = yprev(f) + 2 * value;
  }

  // NB: value needs to be scaled by dy
  void neumann_o3(Field3D& f, BoutReal value) const {
    ASSERT3(valid() >= 0);
    if (valid() < 1) {
      return neumann_o1(f, value);
    }
    ynext(f) =
        parallel_stencil::neumann_o3(1 - length(), value, 1, f[ind()], 2, yprev(f));
  }

  // BoutReal get(const Field3D& f, int off)
  const BoutReal& ynext(const Field3D& f) const { return f.ynext(dir)[ind().yp(dir)]; }
  BoutReal& ynext(Field3D& f) const { return f.ynext(dir)[ind().yp(dir)]; }

  const BoutReal& yprev(const Field3D& f) const {
    ASSERT3(valid() > 0);
    return f.ynext(-dir)[ind().yp(-dir)];
  }
  BoutReal& yprev(Field3D& f) const {
    ASSERT3(valid() > 0);
    return f.ynext(-dir)[ind().yp(-dir)];
  }
  void setYPrevIfValid(Field3D& f, BoutReal val) const {
    if (valid() > 0) {
      yprev(f) = val;
    }
  }

#if BOUT_USE_METRIC_3D == 0
  const BoutReal& ynext(const Field2D& f) const { return f.ynext(dir)[ind().yp(dir)]; }
  BoutReal& ynext(Field2D& f) const { return f.ynext(dir)[ind().yp(dir)]; }

  const BoutReal& yprev(const Field2D& f) const {
    ASSERT3(valid() > 0);
    return f.ynext(-dir)[ind().yp(-dir)];
  }
  BoutReal& yprev(Field2D& f) const {
    ASSERT3(valid() > 0);
    return f.ynext(-dir)[ind().yp(-dir)];
  }
#endif

private:
  const IndicesVec& bndry_points;
  IndicesIter bndry_position;

  constexpr static BoutReal small_value = 1e-2;

public:
  const int dir;
  Mesh* localmesh;
};
} // namespace parallel_boundary_region
} // namespace bout
using BoundaryRegionParIter = bout::parallel_boundary_region::BoundaryRegionParIterBase<
    bout::parallel_boundary_region::IndicesVec,
    bout::parallel_boundary_region::IndicesIter>;
using BoundaryRegionParIterConst =
    bout::parallel_boundary_region::BoundaryRegionParIterBase<
        const bout::parallel_boundary_region::IndicesVec,
        bout::parallel_boundary_region::IndicesIterConst>;

class BoundaryRegionPar : public BoundaryRegionBase {
public:
  BoundaryRegionPar(const std::string& name, int dir, Mesh* passmesh)
      : BoundaryRegionBase(name, passmesh), dir(dir) {
    ASSERT0(std::abs(dir) == 1);
    BoundaryRegionBase::isParallel = true;
  }
  BoundaryRegionPar(const std::string& name, BndryLoc loc, int dir, Mesh* passmesh)
      : BoundaryRegionBase(name, loc, passmesh), dir(dir) {
    BoundaryRegionBase::isParallel = true;
    ASSERT0(std::abs(dir) == 1);
  }

  /// Add a point to the boundary
  void add_point(Ind3D ind, BoutReal x, BoutReal y, BoutReal z, BoutReal length,
                 char valid, signed char offset) {
    bndry_points.push_back({ind,
                            {x, y, z},
                            length,
                            valid,
                            offset,
                            static_cast<unsigned char>(std::abs(offset))});
  }
  void add_point(int ix, int iy, int iz, BoutReal x, BoutReal y, BoutReal z,
                 BoutReal length, char valid, signed char offset) {
    add_point(xyz2ind(ix, iy, iz, localmesh), x, y, z, length, valid, offset);
  }

  // final, so they can be inlined
  void first() final { bndry_position = std::begin(bndry_points); }
  void next() final { ++bndry_position; }
  bool isDone() final { return (bndry_position == std::end(bndry_points)); }

  bool contains(const BoundaryRegionPar& bndry) const {
    return std::binary_search(std::begin(bndry_points), std::end(bndry_points),
                              *bndry.bndry_position,
                              [](const bout::parallel_boundary_region::Indices& i1,
                                 const bout::parallel_boundary_region::Indices& i2) {
                                return i1.index < i2.index;
                              });
  }

  bool contains(const int ix, const int iy, const int iz) const {
    const auto i2 = xyz2ind(ix, iy, iz, localmesh);
    for (auto i1 : bndry_points) {
      if (i1.index == i2) {
        return true;
      }
    }
    return false;
  }
  // setter
  void setValid(char val) { bndry_position->valid = val; }

  // BoundaryRegionParIterConst begin() const {
  //   return BoundaryRegionParIterConst(bndry_points, bndry_points.begin(), dir);
  // }
  // BoundaryRegionParIterConst end() const {
  //   return BoundaryRegionParIterConst(bndry_points, bndry_points.begin(), dir);
  // }
  BoundaryRegionParIter begin() {
    return BoundaryRegionParIter(bndry_points, bndry_points.begin(), dir, localmesh);
  }
  BoundaryRegionParIter end() {
    return BoundaryRegionParIter(bndry_points, bndry_points.end(), dir, localmesh);
  }

  const int dir;

private:
  /// Vector of points in the boundary
  bout::parallel_boundary_region::IndicesVec bndry_points;
  /// Current position in the boundary points
  bout::parallel_boundary_region::IndicesIter bndry_position;

  static Ind3D xyz2ind(int x, int y, int z, Mesh* mesh) {
    const int ny = mesh->LocalNy;
    const int nz = mesh->LocalNz;
    return Ind3D{(x * ny + y) * nz + z, ny, nz};
  }
};

#endif //  BOUT_PAR_BNDRY_H
