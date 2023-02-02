#include "bout/parallel_boundary_region.hxx"

void BoundaryRegionPar::add_point(const int jx, const int jy, const int jz,
                                  const BoutReal x, const BoutReal y, const BoutReal z,
                                  const BoutReal length, const BoutReal angle) {
  bndry_points.push_back({{jx, jy, jz}, {x, y, z}, length, angle});
}

void BoundaryRegionPar::first() {
  bndry_position = begin(bndry_points);
  if (!isDone()) {
    x = bndry_position->index.jx;
    y = bndry_position->index.jy;
    z = bndry_position->index.jz;
    s_x = bndry_position->intersection.s_x;
    s_y = bndry_position->intersection.s_y;
    s_z = bndry_position->intersection.s_z;
    length = bndry_position->length;
    angle = bndry_position->angle;
  }
}

void BoundaryRegionPar::next() {
  ++bndry_position;
  if (!isDone()) {
    x = bndry_position->index.jx;
    y = bndry_position->index.jy;
    z = bndry_position->index.jz;
    s_x = bndry_position->intersection.s_x;
    s_y = bndry_position->intersection.s_y;
    s_z = bndry_position->intersection.s_z;
    length = bndry_position->length;
    angle = bndry_position->angle;
  }
}

bool BoundaryRegionPar::isDone() { return (bndry_position == end(bndry_points)); }
