/**************************************************************************
 * Implements the Flux Coordinate Independent scheme for parallel derivatives
 *
 * A method for field-aligned parallel derivatives which does not require
 * flux-coordinates in the perpendicular direction. Instead, parallel
 * derivatives are taken by following the field line to the adjacent
 * perpendicular planes and interpolating the function onto the field
 * line end-points. The interpolated function can then be used in a
 * finite differencing scheme.
 *
 * IMPORTANT NOTE: The FCI approach requires that the toroidal coordinate
 * be identified with the "parallel" direction. Due to the set-up of
 * BOUT++'s grids, this means that z is now the poloidal coordinate,
 * while y is taken to be toroidal. This means it is probably not
 * possible to just swap in the FCI approach for the standard BOUT++
 * Grad_par operator.
 **************************************************************************
 * Copyright 2014 - 2025 BOUT++ developers
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include "fci.hxx"
#include "bout/parallel_boundary_op.hxx"
#include "bout/parallel_boundary_region.hxx"
#include <bout/bout_types.hxx>
#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout/msg_stack.hxx>
#include <bout/utils.hxx>

#include <string>

namespace {
// Get a unique name for a field based on the sign/magnitude of the offset
std::string parallel_slice_field_name(std::string field, int offset) {
  const std::string direction = (offset > 0) ? "forward" : "backward";
  // We only have a suffix for parallel slices beyond the first
  // This is for backwards compatibility
  const std::string slice_suffix =
      (std::abs(offset) > 1) ? "_" + std::to_string(std::abs(offset)) : "";
  return direction + "_" + field + slice_suffix;
};

#if BOUT_USE_METRIC_3D
void set_parallel_metric_component(std::string name, Field3D& component, int offset,
                                   Field3D& data) {
  if (!component.hasParallelSlices()) {
    component.splitParallelSlices();
    component.allowCalcParallelSlices = false;
  }
  auto& pcom = component.ynext(offset);
  pcom.allocate();
  pcom.setRegion(fmt::format("RGN_YPAR_{:+d}", offset));
  pcom.name = name;
  BOUT_FOR(i, component.getRegion("RGN_NOBNDRY")) { pcom[i.yp(offset)] = data[i]; }
}

bool load_parallel_metric_component(std::string name, Field3D& component, int offset,
                                    bool doZero) {
  Mesh* mesh = component.getMesh();
  Field3D tmp{mesh};
  bool doload = mesh->sourceHasVar(name);
  bool isValid{false};
  if (doload) {
    const auto pname = parallel_slice_field_name(name, offset);
    isValid = mesh->get(tmp, pname, 0.0, false) == 0;
    if (not isValid) {
      throw BoutException("Could not read {:s} from grid file!\n"
                          "Regenerate the grid with a recent zoidberg!",
                          pname);
    }
  } else {
    auto lmin = min(component, true);
    auto lmax = max(component, true);
    if (lmin != lmax) {
      if (doZero) {
        lmin = lmax = 0.0;
      } else {
        throw BoutException("{:s} not in grid file but not constant!\n"
                            "  Cannot determine value for parallel slices.\n"
                            "  Regenerate the grid with a recent zoidberg!",
                            name);
      }
    } else {
      isValid = true;
    }
    tmp = lmin;
  }
  set_parallel_metric_component(name, component, offset, tmp);
  return isValid;
}
#endif

void load_parallel_metric_components([[maybe_unused]] Coordinates* coords,
                                     [[maybe_unused]] int offset) {
#if BOUT_USE_METRIC_3D
#define LOAD_PAR(var, doZero) \
  load_parallel_metric_component(#var, coords->var, offset, doZero)
  LOAD_PAR(g11, false);
  LOAD_PAR(g22, false);
  LOAD_PAR(g33, false);
  LOAD_PAR(g12, false);
  LOAD_PAR(g13, false);
  LOAD_PAR(g23, false);

  LOAD_PAR(g_11, false);
  LOAD_PAR(g_22, false);
  LOAD_PAR(g_33, false);
  LOAD_PAR(g_12, false);
  LOAD_PAR(g_13, false);
  LOAD_PAR(g_23, false);

  LOAD_PAR(dy, false);

  if (not LOAD_PAR(J, true)) {
    auto g =
        coords->g11.ynext(offset) * coords->g22.ynext(offset) * coords->g33.ynext(offset)
        + 2.0 * coords->g12.ynext(offset) * coords->g13.ynext(offset)
              * coords->g23.ynext(offset)
        - coords->g11.ynext(offset) * coords->g23.ynext(offset)
              * coords->g23.ynext(offset)
        - coords->g22.ynext(offset) * coords->g13.ynext(offset)
              * coords->g13.ynext(offset)
        - coords->g33.ynext(offset) * coords->g12.ynext(offset)
              * coords->g12.ynext(offset);

    const auto rgn = fmt::format("RGN_YPAR_{:+d}", offset);
    // Check that g is positive
    bout::checkPositive(g, "The determinant of g^ij", rgn);
    auto J = 1. / sqrt(g);
    auto& pcom = coords->J.ynext(offset);
    BOUT_FOR(i, J.getRegion(rgn)) { pcom[i] = J[i]; }
  }
  if (coords->Bxy.getMesh()->sourceHasVar(parallel_slice_field_name("Bxy", 1))) {
    LOAD_PAR(Bxy, true);
  } else {
    Field3D tmp{coords->Bxy.getMesh()};
    tmp.allocate();
    BOUT_FOR(iyp, coords->Bxy.getRegion("RGN_NOBNDRY")) {
      const auto i = iyp.ym(offset);
      tmp[i] = coords->Bxy[i] * coords->g_22[i] / coords->J[i]
               * coords->J.ynext(offset)[iyp] / coords->g_22.ynext(offset)[iyp];
    }
    set_parallel_metric_component("Bxy", coords->Bxy, offset, tmp);
  }
#undef LOAD_PAR
#endif
}

template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

} // namespace

FCIMap::FCIMap(Mesh& mesh, const Coordinates::FieldMetric& UNUSED(dy), Options& options,
               int offset_, const std::shared_ptr<BoundaryRegionPar>& inner_boundary,
               const std::shared_ptr<BoundaryRegionPar>& outer_boundary, bool zperiodic)
    : map_mesh(mesh), offset(offset_),
      region_no_boundary(map_mesh.getRegion("RGN_NOBNDRY")),
      corner_boundary_mask(map_mesh) {

  TRACE("Creating FCIMap for direction {:d}", offset);

  if (offset == 0) {
    throw BoutException(
        "FCIMap called with offset = 0; You probably didn't mean to do that");
  }

  auto& interpolation_options = options["xzinterpolation"];
  interp =
      XZInterpolationFactory::getInstance().create(&interpolation_options, &map_mesh);
  interp->setYOffset(offset);

  interp_corner =
      XZInterpolationFactory::getInstance().create(&interpolation_options, &map_mesh);
  interp_corner->setYOffset(offset);

  // Index-space coordinates of forward/backward points
  Field3D xt_prime{&map_mesh}, zt_prime{&map_mesh};

  // Real-space coordinates of grid points
  Field3D R{&map_mesh}, Z{&map_mesh};

  // Real-space coordinates of forward/backward points
  Field3D R_prime{&map_mesh}, Z_prime{&map_mesh};

  map_mesh.get(R, "R", 0.0, false);
  map_mesh.get(Z, "Z", 0.0, false);

  // If we can't read in any of these fields, things will silently not
  // work, so best throw
  if (map_mesh.get(xt_prime, parallel_slice_field_name("xt_prime", offset), 0.0, false)
      != 0) {
    throw BoutException("Could not read {:s} from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("xt_prime", offset));
  }
  if (map_mesh.get(zt_prime, parallel_slice_field_name("zt_prime", offset), 0.0, false)
      != 0) {
    throw BoutException("Could not read {:s} from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("zt_prime", offset));
  }
  if (map_mesh.get(R_prime, parallel_slice_field_name("R", offset), 0.0, false) != 0) {
    throw BoutException("Could not read {:s} from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("R", offset));
  }
  if (map_mesh.get(Z_prime, parallel_slice_field_name("Z", offset), 0.0, false) != 0) {
    throw BoutException("Could not read {:s} from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("Z", offset));
  }

  // Cell corners
  Field3D xt_prime_corner{emptyFrom(xt_prime)};
  Field3D zt_prime_corner{emptyFrom(zt_prime)};

  BOUT_FOR(i, xt_prime_corner.getRegion("RGN_NOBNDRY")) {
    // Point interpolated from (x+1/2, z+1/2)

    // Cache the offsets
    auto i_xplus = i.xp();
    auto i_zplus = i.zp();
    auto i_xzplus = i_zplus.xp();

    if ((xt_prime[i] < 0.0) || (xt_prime[i_xplus] < 0.0) || (xt_prime[i_xzplus] < 0.0)
        || (xt_prime[i_zplus] < 0.0)) {
      // Hit a boundary
      corner_boundary_mask(i.x(), i.y(), i.z()) = true;

      xt_prime_corner[i] = -1.0;
      zt_prime_corner[i] = -1.0;
    } else {
      xt_prime_corner[i] =
          0.25
          * (xt_prime[i] + xt_prime[i_xplus] + xt_prime[i_zplus] + xt_prime[i_xzplus]);

      zt_prime_corner[i] =
          0.25
          * (zt_prime[i] + zt_prime[i_xplus] + zt_prime[i_zplus] + zt_prime[i_xzplus]);
    }
  }

  interp_corner->setMask(corner_boundary_mask);

  {
    TRACE("FCImap: calculating corner weights");
    interp_corner->calcWeights(xt_prime_corner, zt_prime_corner);
  }

  {
    TRACE("FCImap: calculating weights");
    interp->calcWeights(xt_prime, zt_prime);
  }

  const int ncz = map_mesh.LocalNz;

  BoutMask to_remove(map_mesh);
  const int xend =
      map_mesh.xstart + (map_mesh.xend - map_mesh.xstart + 1) * map_mesh.getNXPE() - 1;
  // Default to the maximum number of points
  const int defValid{map_mesh.ystart - 1 + std::abs(offset)};
  // Serial loop because call to BoundaryRegionPar::addPoint
  // (probably?) can't be done in parallel
  BOUT_FOR_SERIAL(i, xt_prime.getRegion("RGN_NOBNDRY")) {
    // z is periodic, so make sure the z-index wraps around
    if (zperiodic) {
      zt_prime[i] = zt_prime[i]
                    - ncz * (static_cast<int>(zt_prime[i] / static_cast<BoutReal>(ncz)));

      if (zt_prime[i] < 0.0) {
        zt_prime[i] += ncz;
      }
    }

    if ((xt_prime[i] >= map_mesh.xstart) and (xt_prime[i] <= xend)) {
      // Not a boundary
      continue;
    }

    const auto x = i.x();
    const auto y = i.y();
    const auto z = i.z();

    //----------------------------------------
    // Boundary stuff
    //
    // If a field line leaves the domain, then the forward or backward
    // indices (forward/backward_xt_prime and forward/backward_zt_prime)
    // are set to -1

    to_remove(x, y, z) = true;

    // Need to specify the index of the boundary intersection, but
    // this may not be defined in general.
    // We do however have the real-space (R,Z) coordinates. Here we extrapolate,
    // using the change in R and Z to calculate the change in (x,z) indices
    //
    // ( dR ) = ( dR/dx  dR/dz ) ( dx )
    // ( dZ )   ( dZ/dx  dZ/dz ) ( dz )
    //
    // where (dR,dZ) is the change in (R,Z) along the field,
    // (dx,dz) is the change in (x,z) index along the field,
    // and the gradients dR/dx etc. are evaluated at (x,y,z)

    // Cache the offsets
    const auto i_zp = i.zp();
    const auto i_zm = i.zm();

    const BoutReal dR_dx = 0.5 * (R[i.xp()] - R[i.xm()]);
    const BoutReal dZ_dx = 0.5 * (Z[i.xp()] - Z[i.xm()]);

    const BoutReal dR_dz = 0.5 * (R[i_zp] - R[i_zm]);
    const BoutReal dZ_dz = 0.5 * (Z[i_zp] - Z[i_zm]);

    const BoutReal det = dR_dx * dZ_dz - dR_dz * dZ_dx; // Determinant of 2x2 matrix

    const BoutReal dR = R_prime[i] - R[i];
    const BoutReal dZ = Z_prime[i] - Z[i];

    // Invert 2x2 matrix to get change in index
    const BoutReal dx = (dZ_dz * dR - dR_dz * dZ) / det;
    const BoutReal dz = (dR_dx * dZ - dZ_dx * dR) / det;

    // Negative xt_prime means we've hit the inner boundary, otherwise the
    // outer boundary. However, if any of the surrounding points are negative,
    // that also means inner. So to differentiate between inner and outer we
    // need at least 2 points in the domain.
    ASSERT2(map_mesh.xend - map_mesh.xstart >= 2);
    auto boundary = (xt_prime[i] < map_mesh.xstart) ? inner_boundary : outer_boundary;
    if (!boundary->contains(x, y, z)) {
      boundary->add_point(x, y, z, x + dx, y + offset - sgn(offset) * 0.5,
                          z + dz, // Intersection point in local index space
                          std::abs(offset) - 0.5, // Distance to intersection
                          defValid, offset);
    }
  }
  region_no_boundary = region_no_boundary.mask(to_remove);

  interp->setRegion(region_no_boundary);

  const auto region = fmt::format("RGN_YPAR_{:+d}", offset);
  if (not map_mesh.hasRegion3D(region)) {
    // The valid region for this slice
    map_mesh.addRegion3D(
        region, Region<Ind3D>(map_mesh.xstart, map_mesh.xend, map_mesh.ystart + offset,
                              map_mesh.yend + offset, 0, map_mesh.LocalNz - 1,
                              map_mesh.LocalNy, map_mesh.LocalNz));
  }
}

Field3D FCIMap::integrate(Field3D& f) const {
  TRACE("FCIMap::integrate");

  ASSERT1(f.getDirectionY() == YDirectionType::Standard);
  ASSERT1(&map_mesh == f.getMesh());

  // Cell centre values
  Field3D centre = interp->interpolate(f);

  // Cell corner values (x+1/2, z+1/2)
  Field3D corner = interp_corner->interpolate(f);

  Field3D result{emptyFrom(f)};
#if CHECK > 2
  // The more general version of invalidate guards
  result = BoutNaN;
#endif

  BOUT_FOR(i, region_no_boundary) {
    const auto inext = i.yp(offset);
    const BoutReal f_c = centre[inext];
    const auto izm = i.zm();
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();
    const int zm = izm.z();
    if (corner_boundary_mask(x, y, z) || corner_boundary_mask(x - 1, y, z)
        || corner_boundary_mask(x, y, zm) || corner_boundary_mask(x - 1, y, zm)
        || (x == map_mesh.xstart)) {
      // One of the corners leaves the domain.
      // Use the cell centre value, since boundary conditions are not
      // currently applied to corners.
      result[inext] = f_c;
    } else {
      const BoutReal f_pp = corner[inext];           // (x+1/2, z+1/2)
      const BoutReal f_mp = corner[inext.xm()];      // (x-1/2, z+1/2)
      const BoutReal f_pm = corner[inext.zm()];      // (x+1/2, z-1/2)
      const BoutReal f_mm = corner[inext.xm().zm()]; // (x-1/2, z-1/2)

      // This uses a simple weighted average of centre and corners
      // A more sophisticated approach might be to use e.g. Gauss-Lobatto points
      // which would include cell edges and corners
      result[inext] = 0.5 * (f_c + 0.25 * (f_pp + f_mp + f_pm + f_mm));
    }
    ASSERT2(std::isfinite(result[inext]));
  }
  return result;
}

void FCITransform::checkInputGrid() {
  std::string parallel_transform;
  if (mesh.isDataSourceGridFile()
      && !mesh.get(parallel_transform, "parallel_transform")) {
    if (parallel_transform != "fci") {
      throw BoutException(
          "Incorrect parallel transform type '" + parallel_transform
          + "' used "
            "to generate metric components for FCITransform. Should be 'fci'.");
    }
  } // else: parallel_transform variable not found in grid input, indicates older input
    //       file or grid from options so must rely on the user having ensured the type is
    //       correct
}

void FCITransform::calcParallelSlices(Field3D& f) {
  TRACE("FCITransform::calcParallelSlices");

  ASSERT1(f.allowCalcParallelSlices);

  ASSERT1(f.getDirectionY() == YDirectionType::Standard);
  // Only have forward_map/backward_map for CELL_CENTRE, so can only deal with
  // CELL_CENTRE inputs
  ASSERT1(f.getLocation() == CELL_CENTRE);

  // Ensure that yup and ydown are different fields
  f.splitParallelSlices();

  // Interpolate f onto yup and ydown fields
  for (const auto& map : field_line_maps) {
    f.ynext(map.offset) = map.interpolate(f);
    f.ynext(map.offset).setRegion(fmt::format("RGN_YPAR_{:+d}", map.offset));
  }
}

void FCITransform::integrateParallelSlices(Field3D& f) {
  TRACE("FCITransform::integrateParallelSlices");

  ASSERT1(f.getDirectionY() == YDirectionType::Standard);
  // Only have forward_map/backward_map for CELL_CENTRE, so can only deal with
  // CELL_CENTRE inputs
  ASSERT1(f.getLocation() == CELL_CENTRE);

  // Ensure that yup and ydown are different fields
  f.splitParallelSlices();

  // Integrate f onto yup and ydown fields
  for (const auto& map : field_line_maps) {
    f.ynext(map.offset) = map.integrate(f);
  }
}

void FCITransform::loadParallelMetrics(Coordinates* coords) {
  for (int i = 1; i <= mesh.ystart; ++i) {
    load_parallel_metric_components(coords, -i);
    load_parallel_metric_components(coords, i);
  }
}

void FCITransform::outputVars(Options& output_options) {
  // Real-space coordinates of grid points
  output_options["R"].force(R, "FCI");
  output_options["Z"].force(Z, "FCI");
}
