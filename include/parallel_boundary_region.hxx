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
class BoundaryRegionPar {

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
  BoundaryRegionPar(const string &name, int dir, Mesh* passmesh) :
    localmesh(passmesh ? passmesh : mesh), label(std::move(name)), dir(dir) {}
  BoundaryRegionPar(const string &name, BndryLoc loc,int dir, Mesh* passmesh) :
    localmesh(passmesh ? passmesh : mesh), label(std::move(name)), location(loc), dir(dir) {}

  /// Add a point to the boundary
  void add_point(int jx,int jy,int jz,
                 const BoutReal x,BoutReal y,BoutReal z,
                 const BoutReal length,BoutReal angle);

  void first();  ///< Move the region iterator to the start
  void next();   ///< Get the next element in the loop
                 ///  over every element from inside out (in
                 ///  X or Y first)
  bool isDone(); ///< Returns true if outside domain. Can use this with nested nextX, nextY

  Mesh* localmesh; ///< Mesh does this boundary region belongs to

  string label; ///< Label for this boundary region

  BndryLoc location;         ///< Which side of the domain is it on?

  /// Index of the point in the boundary
  int x, y, z;
  BoutReal s_x, s_y, s_z;
  BoutReal length;
  BoutReal angle;

  const int dir;
};

#endif //  __PAR_BNDRY_H__
