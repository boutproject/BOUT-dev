#include <fci_boundary_region.hxx>

void BoundaryRegionFCI::add_point(const int x, const int y, const int z) {
  bndry_points.push_back({x, y, z});
}

void BoundaryRegionFCI::first() {
  bndry_position = bndry_points.begin();
  if (!isDone()) {
	x = bndry_position->x;
	y = bndry_position->y;
	z = bndry_position->z;
  }
}

void BoundaryRegionFCI::next() {
  ++bndry_position;
  if (!isDone()) {
	x = bndry_position->x;
	y = bndry_position->y;
	z = bndry_position->z;
  }
}

bool BoundaryRegionFCI::isDone() {
  return (bndry_position == bndry_points.end());
}
