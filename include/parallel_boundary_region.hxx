#ifndef __PAR_BNDRY_H__
#define __PAR_BNDRY_H__

#include "boundary_region.hxx"
#include "bout_types.hxx"
#include <vector>

/**
 * Boundary region for parallel direction. This contains a vector of points that are
 * inside the boundary.
 *
 */
class BoundaryRegionPar : public BoundaryRegionBase {

  struct IndexPoint {
    int jx;
    int jy;
    int jz;
  };

  struct RealPoint {
    BoutReal s_x;
    BoutReal s_y;
    BoutReal s_z;
  };

  struct Indices {
    // Indices of the boundary point
    IndexPoint index;
    // Intersection with boundary in index space
    RealPoint intersection;
    // Distance to intersection
    BoutReal length;
    // Angle between field line and boundary
    BoutReal angle;
  };

  typedef std::vector<Indices> IndicesVec;
  typedef IndicesVec::iterator IndicesIter;

  /// Vector of points in the boundary
  IndicesVec bndry_points;
  /// Current position in the boundary points
  IndicesIter bndry_position;

public:
  BoundaryRegionPar(const string &name,int dir) :
    BoundaryRegionBase(name), dir(dir) {
    BoundaryRegionBase::isParallel = true;}
  BoundaryRegionPar(const string &name, BndryLoc loc,int dir) :
    BoundaryRegionBase(name, loc), dir(dir) {
    BoundaryRegionBase::isParallel = true;}

  /// Add a point to the boundary
  void add_point(int jx,int jy,int jz,
                 BoutReal x,BoutReal y,BoutReal z,
                 BoutReal length,BoutReal angle);

  void first();
  void next();
  bool isDone();

  /// Index of the point in the boundary
  int x, y, z;
  BoutReal s_x, s_y, s_z;
  BoutReal length;
  BoutReal angle;

  const int dir;
};

#endif //  __PAR_BNDRY_H__
