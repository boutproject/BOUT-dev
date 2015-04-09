#include "parallel_boundary_region.hxx"

void BoundaryRegionPar::add_point(const int x, const int y, const int z, const BoutReal length, const BoutReal angle) {
  bndry_points.push_back({x, y, z, length, angle});
}

void BoundaryRegionPar::first() {
  bndry_position = begin(bndry_points);
  if (!isDone()) {
    x      = bndry_position->x;
    y      = bndry_position->y;
    z      = bndry_position->z;
    length = bndry_position->length;
    angle  = bndry_position->angle;
  }
}

void BoundaryRegionPar::next() {
  ++bndry_position;
  if (!isDone()) {
    x      = bndry_position->x;
    y      = bndry_position->y;
    z      = bndry_position->z;
    length = bndry_position->length;
    angle  = bndry_position->angle;
  }
}

bool BoundaryRegionPar::isDone() {
  return (bndry_position == end(bndry_points));
}

