#pragma once

#include "./boundary_iterator.hxx"
#include "bout/assert.hxx"
#include "bout/boutexception.hxx"
#include "bout/field_data.hxx"
#include "bout/globals.hxx"
#include "bout/mask.hxx"
#include "bout/options.hxx"
#include "bout/parallel_boundary_region.hxx"
#include "bout/region.hxx"
#include <memory>
#include <vector>

/// This class allows to simplify iterating over y-boundaries.
///
/// It makes it easier to write code for FieldAligned boundaries, but if a bit
/// care is taken the code also works with FluxCoordinateIndependent code.
///
/// An example how to replace old code is given here:
/// ../../manual/sphinx/user_docs/boundary_options.rst

class YBoundary {
public:
  template <class Func>
  void iter_regions(const Func& func) {
    for (auto& region : boundary_regions) {
      func(*region);
    }
    for (auto& region : boundary_regions_par) {
      func(*region);
    }
  }
  template <class Func>
  void iter_points(const Func& func) {
    iter_regions([&](auto& region) {
      for (auto& point : region) {
        func(point);
      }
    });
  }

  template <class Func>
  void iter(const Func& func) {
    return iter_points(func);
  }

  YBoundary(YBndryType type, Options* options_ptr, const Mesh& mesh) {
    bool lower_y = true;
    bool upper_y = true;
    bool outer_x = true;
    bool inner_x = false;
    if (options_ptr != nullptr) {
      auto& options = *options_ptr;
      if (!mesh.isFci()) {
        lower_y =
            options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(lower_y);
        upper_y =
            options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(upper_y);
      } else {
        outer_x =
            options["outer_x"].doc("Boundary on outer x?").withDefault<bool>(outer_x);
        inner_x =
            options["inner_x"].doc("Boundary on inner x?").withDefault<bool>(inner_x);
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
      if (lower_y) {
        boundary_regions.push_back(
            std::make_shared<NewBoundaryRegionY>(mesh, true, mesh.iterateBndryLowerY()));
      }
      if (upper_y) {
        boundary_regions.push_back(
            std::make_shared<NewBoundaryRegionY>(mesh, false, mesh.iterateBndryUpperY()));
      }
    }
    // Cache boundary regions
    _contains.emplace_back(&mesh, false);
    _contains.emplace_back(&mesh, false);
    iter_points([&](const auto& point) {
      if (point.dir() == 1) {
        _contains[1][point.ind()] = true;
      } else if (point.dir() == -1) {
        _contains[0][point.ind()] = true;
      }
    });
  }

  bool contains_ylow(Ind3D ind) const { return _contains[0][ind]; }
  bool contains_yhigh(Ind3D ind) const { return _contains[1][ind]; }
  template <int dir>
  bool contains(Ind3D ind) const {
    static_assert(dir == 1 || dir == -1);
    if constexpr (dir == 1) {
      return _contains[1][ind];
    }
    if constexpr (dir == -1) {
      return _contains[0][ind];
    }
  }
  bool contains(int dir, Ind3D ind) const {
    if (dir == 1) {
      return contains<+1>(ind);
    }
    if (dir == -1) {
      return contains<-1>(ind);
    }
    throw BoutException("only dir == 1 and dir == -1 are implemented, not {}", dir);
  }

private:
  std::vector<std::shared_ptr<BoundaryRegionPar>> boundary_regions_par;
  std::vector<std::shared_ptr<NewBoundaryRegionY>> boundary_regions;

  std::vector<BoutMask> _contains;
};

inline std::shared_ptr<YBoundary> getYBoundary(Coordinates* coords,
                                               YBndryType type = YBndryType::sheath) {
  auto itype = static_cast<int>(type);
  if (coords->ybndrys[itype] == nullptr) {
    coords->ybndrys[itype] = coords->makeYBoundary(type);
  }
  return coords->ybndrys[itype];
}
