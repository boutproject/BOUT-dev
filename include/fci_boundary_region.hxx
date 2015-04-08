#ifndef __FCI_BNDRY_H__
#define __FCI_BNDRY_H__

#include "boundary_region.hxx"
#include "bout_types.hxx"
#include <vector>

/**
 * Boundary region for FCI. This contains a vector of points that are
 * inside the boundary.
 *
 */
class BoundaryRegionFCI : public BoundaryRegion {

  struct Indices {
    int x;
    int y;
    int z;
    BoutReal length;
    BoutReal angle;
  };

  typedef std::vector<Indices> IndicesVec;
  typedef IndicesVec::iterator IndicesIter;

  /// Vector of points in the boundary
  IndicesVec bndry_points;
  /// Current position in the boundary points
  IndicesIter bndry_position;

public:
  BoundaryRegionFCI(const string &name, BndryLoc loc, const int dir) :
    BoundaryRegion(name, loc), dir(dir) {}

  /// Add a point to the boundary
  void add_point(const int x, const int y, const int z, const BoutReal length, const BoutReal angle);

  void first();
  void next();
  void next1d() {}
  void nextX() {}
  void nextY() {}
  bool isDone();

  /// Index of the point in the boundary
  int z;
  BoutReal length;
  BoutReal angle;

  const int dir;
};

#endif //  __FCI_BNDRY_H__
