#include "bout/yboundary_regions.hxx"

#include "bout/boundary_region_iter.hxx"
#include "bout/field_data.hxx"
#include "bout/mesh.hxx"

#include <memory>

namespace bout::boundary {
YBoundary::YBoundary(const Mesh& mesh, bool lower_y, bool upper_y)
    : _contains_low(&mesh, false), _contains_high(&mesh, false) {

  if (mesh.isFci()) {
    if (lower_y) {
      for (auto& bndry : mesh.getBoundariesPar(BoundaryParType::xout)) {
        boundary_regions_par.push_back(bndry);
      }
    }
    if (upper_y) {
      for (auto& bndry : mesh.getBoundariesPar(BoundaryParType::xin)) {
        boundary_regions_par.push_back(bndry);
      }
    }
  } else {
    for (auto& bndry : mesh.getBoundaries()) {
      if ((lower_y && bndry->location == BndryLoc::ydown)
          or (upper_y && bndry->location == BndryLoc::yup)) {
        boundary_regions.push_back(
            std::dynamic_pointer_cast<bout::boundary::BoundaryRegionY>(bndry));
      }
    }
  }

  // Cache boundary regions
  iter([&](const BoundaryIterator auto& point) {
    if (point.dir() == 1) {
      _contains_high[point.ind()] = true;
    } else if (point.dir() == -1) {
      _contains_low[point.ind()] = true;
    }
  });
}
} // namespace bout::boundary
