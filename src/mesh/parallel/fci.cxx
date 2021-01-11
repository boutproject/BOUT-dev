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
 * Copyright 2014 B.D.Dudson, P. Hill
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
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
#include "parallel_boundary_op.hxx"
#include "parallel_boundary_region.hxx"
#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout_types.hxx>
#include <msg_stack.hxx>
#include <utils.hxx>

#include <string>

FCIMap::FCIMap(Mesh& mesh, const Coordinates::FieldMetric& dy, Options& options,
               int offset_, BoundaryRegionPar* inner_boundary,
               BoundaryRegionPar* outer_boundary, bool zperiodic)
    : map_mesh(mesh), offset(offset_), boundary_mask(map_mesh),
      corner_boundary_mask(map_mesh) {

  TRACE("Creating FCIMAP for direction {:d}", offset);

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

  // Get a unique name for a field based on the sign/magnitude of the offset
  const auto parallel_slice_field_name = [&](std::string field) -> std::string {
    const std::string direction = (offset > 0) ? "forward" : "backward";
    // We only have a suffix for parallel slices beyond the first
    // This is for backwards compatibility
    const std::string slice_suffix =
        (std::abs(offset) > 1) ? "_" + std::to_string(std::abs(offset)) : "";
    return direction + "_" + field + slice_suffix;
  };

  // If we can't read in any of these fields, things will silently not
  // work, so best throw
  if (map_mesh.get(xt_prime, parallel_slice_field_name("xt_prime"), 0.0, false) != 0) {
    throw BoutException("Could not read {:s} from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("xt_prime"));
  }
  if (map_mesh.get(zt_prime, parallel_slice_field_name("zt_prime"), 0.0, false) != 0) {
    throw BoutException("Could not read {:s} from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("zt_prime"));
  }
  if (map_mesh.get(R_prime, parallel_slice_field_name("R"), 0.0, false) != 0) {
    throw BoutException("Could not read {:s} from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("R"));
  }
  if (map_mesh.get(Z_prime, parallel_slice_field_name("Z"), 0.0, false) != 0) {
    throw BoutException("Could not read {:s} from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("Z"));
  }

  // Cell corners
  Field3D xt_prime_corner{emptyFrom(xt_prime)};
  Field3D zt_prime_corner{emptyFrom(xt_prime)};

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
      continue;
    }

    xt_prime_corner[i] =
        0.25 * (xt_prime[i] + xt_prime[i_xplus] + xt_prime[i_zplus] + xt_prime[i_xzplus]);

    zt_prime_corner[i] =
        0.25 * (zt_prime[i] + zt_prime[i_xplus] + zt_prime[i_zplus] + zt_prime[i_xzplus]);
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

    if ((xt_prime[i] >= 0.0) or (xt_prime[i] <= map_mesh.xend)) {
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

    boundary_mask(x, y, z) = true;

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
    const auto i_xp = i.xp();
    const auto i_xm = i.xm();
    const auto i_zp = i.zp();
    const auto i_zm = i.zm();

    const BoutReal dR_dx = 0.5 * (R[i_xp] - R[i_xm]);
    const BoutReal dZ_dx = 0.5 * (Z[i_xp] - Z[i_xm]);

    BoutReal dR_dz, dZ_dz;
    // Handle the edge cases in Z
    if (z == 0) {
      dR_dz = R[i_zp] - R[i];
      dZ_dz = Z[i_zp] - Z[i];

    } else if (z == map_mesh.LocalNz - 1) {
      dR_dz = R[i] - R[i_zm];
      dZ_dz = Z[i] - Z[i_zm];

    } else {
      dR_dz = 0.5 * (R[i_zp] - R[i_zm]);
      dZ_dz = 0.5 * (Z[i_zp] - Z[i_zm]);
    }

    const BoutReal det = dR_dx * dZ_dz - dR_dz * dZ_dx; // Determinant of 2x2 matrix

    const BoutReal dR = R_prime[i] - R[i];
    const BoutReal dZ = Z_prime[i] - Z[i];

    // Invert 2x2 matrix to get change in index
    const BoutReal dx = (dZ_dz * dR - dR_dz * dZ) / det;
    const BoutReal dz = (dR_dx * dZ - dZ_dx * dR) / det;

    // Negative xt_prime means we've hit the inner boundary, otherwise
    // the outer boundary
    auto* boundary = (xt_prime[i] < 0.0) ? inner_boundary : outer_boundary;
    boundary->add_point(x, y, z, x + dx, y + 0.5 * offset,
                        z + dz,      // Intersection point in local index space
                        0.5 * dy[i], // Distance to intersection
                        PI           // Right-angle intersection
    );
  }

  interp->setMask(boundary_mask);
}

Field3D FCIMap::integrate(Field3D &f) const {
  TRACE("FCIMap::integrate");

  ASSERT1(f.getDirectionY() == YDirectionType::Standard);
  ASSERT1(&map_mesh == f.getMesh());

  // Cell centre values
  Field3D centre = interp->interpolate(f);

  // Cell corner values (x+1/2, z+1/2)
  Field3D corner = interp_corner->interpolate(f);

  Field3D result{emptyFrom(f)};

  int nz = map_mesh.LocalNz;

  for(int x = map_mesh.xstart; x <= map_mesh.xend; x++) {
    for(int y = map_mesh.ystart; y <= map_mesh.yend; y++) {

      int ynext = y+offset;

      for(int z = 0; z < nz; z++) {
        if (boundary_mask(x,y,z))
          continue;

        int zm = z - 1;
        if (z == 0) {
          zm = nz-1;
        }

        BoutReal f_c  = centre(x,ynext,z);

        if (corner_boundary_mask(x, y, z) || corner_boundary_mask(x - 1, y, z) ||
            corner_boundary_mask(x, y, zm) || corner_boundary_mask(x - 1, y, zm) ||
            (x == map_mesh.xstart)) {
          // One of the corners leaves the domain.
          // Use the cell centre value, since boundary conditions are not
          // currently applied to corners.
          result(x, ynext, z) = f_c;

        } else {
          BoutReal f_pp = corner(x, ynext, z);      // (x+1/2, z+1/2)
          BoutReal f_mp = corner(x - 1, ynext, z);  // (x-1/2, z+1/2)
          BoutReal f_pm = corner(x, ynext, zm);     // (x+1/2, z-1/2)
          BoutReal f_mm = corner(x - 1, ynext, zm); // (x-1/2, z-1/2)

          // This uses a simple weighted average of centre and corners
          // A more sophisticated approach might be to use e.g. Gauss-Lobatto points
          // which would include cell edges and corners
          result(x, ynext, z) = 0.5 * (f_c + 0.25 * (f_pp + f_mp + f_pm + f_mm));

          ASSERT2(finite(result(x,ynext,z)));
        }
      }
    }
  }
  return result;
}

void FCITransform::checkInputGrid() {
  std::string parallel_transform;
  if (mesh.isDataSourceGridFile() && !mesh.get(parallel_transform, "parallel_transform")) {
    if (parallel_transform != "fci") {
      throw BoutException("Incorrect parallel transform type '"+parallel_transform+"' used "
          "to generate metric components for FCITransform. Should be 'fci'.");
    }
  } // else: parallel_transform variable not found in grid input, indicates older input
    //       file or grid from options so must rely on the user having ensured the type is
    //       correct
}

void FCITransform::calcParallelSlices(Field3D& f) {
  TRACE("FCITransform::calcParallelSlices");

  ASSERT1(f.getDirectionY() == YDirectionType::Standard);
  // Only have forward_map/backward_map for CELL_CENTRE, so can only deal with
  // CELL_CENTRE inputs
  ASSERT1(f.getLocation() == CELL_CENTRE);

  // Ensure that yup and ydown are different fields
  f.splitParallelSlices();

  // Interpolate f onto yup and ydown fields
  for (const auto& map : field_line_maps) {
    f.ynext(map.offset) = map.interpolate(f);
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
