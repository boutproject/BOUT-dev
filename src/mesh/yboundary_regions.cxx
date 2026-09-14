#include "bout/yboundary_regions.hxx"

#include "bout/boundary_region_iter.hxx"
#include "bout/field_data.hxx"
#include "bout/mesh.hxx"

#include <memory>

namespace bout::boundary {
YBoundary::YBoundary(const Mesh& mesh, bool lower, bool upper)
    : lower(lower), upper(upper), _contains_lower(&mesh, false),
      _contains_upper(&mesh, false) {

  if (mesh.isFci()) {
    if (lower) {
      for (auto& bndry : mesh.getBoundariesPar(BoundaryParType::xout)) {
        boundary_regions_par.push_back(bndry);
      }
    }
    if (upper) {
      for (auto& bndry : mesh.getBoundariesPar(BoundaryParType::xin)) {
        boundary_regions_par.push_back(bndry);
      }
    }
  } else {
    for (auto& bndry : mesh.getBoundaries()) {
      if ((lower && bndry->location == BndryLoc::ydown)
          or (upper && bndry->location == BndryLoc::yup)) {
        boundary_regions.push_back(
            std::dynamic_pointer_cast<bout::boundary::BoundaryRegionY>(bndry));
      }
    }
  }

  // Cache boundary regions
  iter([&](const BoundaryIterator auto& point) {
    if (point.dir() == 1) {
      _contains_upper[point.ind()] = true;
    } else if (point.dir() == -1) {
      _contains_lower[point.ind()] = true;
    }
  });
}
} // namespace bout::boundary
