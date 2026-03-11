#pragma once

#include "./boundary_iterator.hxx"
#include "bout/assert.hxx"
#include "bout/boutexception.hxx"
#include "bout/globals.hxx"
#include "bout/mask.hxx"
#include "bout/options.hxx"
#include "bout/parallel_boundary_region.hxx"
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
  template <class F>
  void iter_regions(const F& f) {
    ASSERT1(is_init);
    for (auto& region : boundary_regions) {
      f(*region);
    }
    for (auto& region : boundary_regions_par) {
      f(*region);
    }
  }
  template <class F>
  void iter_pnts(const F& f) {
    iter_regions([&](auto& region) {
      for (auto& pnt : region) {
        f(pnt);
      }
    });
  }

  template <class F>
  void iter(const F& f) {
    return iter_regions(f);
  }

  void init(Options& options, Mesh* mesh = nullptr) {
    if (is_init) {
      if (optptr == &options) {
        return;
      }
      throw BoutException("YBoundary is already initialised!");
    }
    if (mesh == nullptr) {
      mesh = bout::globals::mesh;
    }

    const bool lower_y =
        options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
    const bool upper_y =
        options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);
    const bool outer_x =
        options["outer_x"].doc("Boundary on outer x?").withDefault<bool>(true);
    const bool inner_x =
        options["inner_x"].doc("Boundary on inner x?").withDefault<bool>(false);

    if (mesh->isFci()) {
      if (outer_x) {
        for (auto& bndry : mesh->getBoundariesPar(BoundaryParType::xout)) {
          boundary_regions_par.push_back(bndry);
        }
      }
      if (inner_x) {
        for (auto& bndry : mesh->getBoundariesPar(BoundaryParType::xin)) {
          boundary_regions_par.push_back(bndry);
        }
      }
    } else {
      if (lower_y) {
        boundary_regions.push_back(
            std::make_shared<NewBoundaryRegionY>(mesh, true, mesh->iterateBndryLowerY()));
      }
      if (upper_y) {
        boundary_regions.push_back(std::make_shared<NewBoundaryRegionY>(
            mesh, false, mesh->iterateBndryUpperY()));
      }
    }
    is_init = true;
    optptr = &options;
    // Cache boundary regions
    _contains.emplace_back(mesh, false);
    _contains.emplace_back(mesh, false);
    iter_pnts([&](const auto& pnt) {
      if (pnt.dir == 1) {
        _contains[1][pnt.ind()] = true;
      } else if (pnt.dir == -1) {
        _contains[0][pnt.ind()] = true;
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
  bool is_init{false};
  Options* optptr;
};
