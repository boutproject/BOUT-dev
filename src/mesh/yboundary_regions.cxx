#include "bout/yboundary_regions.hxx"

#include "bout/boundary_region.hxx"
#include "bout/boundary_region_iter.hxx"
#include "bout/field_data.hxx"
#include "bout/options.hxx"

#include <memory>

namespace bout::boundary {
YBoundary::YBoundary(YBndryType type, Options* options_ptr, const Mesh& mesh)
    : _contains_low(&mesh, false), _contains_high(&mesh, false) {
  bool lower_y = true;
  bool upper_y = true;
  bool outer_x = true;
  bool inner_x = false;
  if (options_ptr != nullptr) {
    auto& options = *options_ptr;
    if (!mesh.isFci()) {
      lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(lower_y);
      upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(upper_y);
    } else {
      outer_x = options["outer_x"].doc("Boundary on outer x?").withDefault<bool>(outer_x);
      inner_x = options["inner_x"].doc("Boundary on inner x?").withDefault<bool>(inner_x);
    }
  }
  switch (type) {
  case YBndryType::sheath:
    break;
  case YBndryType::not_sheath:
    lower_y = !lower_y;
    upper_y = !upper_y;
    outer_x = !outer_x;
    inner_x = !inner_x;
    break;
  case YBndryType::all:
    lower_y = true;
    upper_y = true;
    outer_x = true;
    inner_x = true;
  }

  if (mesh.isFci()) {
    if (outer_x) {
      for (auto& bndry : mesh.getBoundariesPar(BoundaryParType::xout)) {
        boundary_regions_par.push_back(bndry);
      }
    }
    if (inner_x) {
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
