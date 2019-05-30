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
#include "interpolation_factory.hxx"
#include "parallel_boundary_op.hxx"
#include "parallel_boundary_region.hxx"
#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout_types.hxx>
#include <msg_stack.hxx>
#include <utils.hxx>

#include <string>

FCIMap::FCIMap(Mesh& mesh, int offset_, BoundaryRegionPar* boundary, bool zperiodic)
    : map_mesh(mesh), offset(offset_), boundary_mask(map_mesh),
      corner_boundary_mask(map_mesh) {

  TRACE("Creating FCIMAP for direction %d", offset);

  if (offset == 0) {
    throw BoutException("FCIMap called with offset = 0; You probably didn't mean to do that");
  }

  interp =
      std::unique_ptr<Interpolation>(InterpolationFactory::getInstance()->create(&map_mesh));
  interp->setYOffset(offset);

  interp_corner =
      std::unique_ptr<Interpolation>(InterpolationFactory::getInstance()->create(&map_mesh));
  interp_corner->setYOffset(offset);

  // Index arrays contain guard cells in order to get subscripts right
  // x-index of bottom-left grid point
  auto i_corner = Tensor<int>(map_mesh.LocalNx, map_mesh.LocalNy, map_mesh.LocalNz);
  // z-index of bottom-left grid point
  auto k_corner = Tensor<int>(map_mesh.LocalNx, map_mesh.LocalNy, map_mesh.LocalNz);

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
    throw BoutException("Could not read %s from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("xt_prime").c_str());
  }
  if (map_mesh.get(zt_prime, parallel_slice_field_name("zt_prime"), 0.0, false) != 0) {
    throw BoutException("Could not read %s from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("zt_prime").c_str());
  }
  if (map_mesh.get(R_prime, parallel_slice_field_name("R"), 0.0, false) != 0) {
    throw BoutException("Could not read %s from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("R").c_str());
  }
  if (map_mesh.get(Z_prime, parallel_slice_field_name("Z"), 0.0, false) != 0) {
    throw BoutException("Could not read %s from grid file!\n"
                        "  Either add it to the grid file, or reduce MYG",
                        parallel_slice_field_name("Z").c_str());
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

    if ((xt_prime[i] < 0.0) || (xt_prime[i_xplus] < 0.0) || (xt_prime[i_xzplus] < 0.0) ||
        (xt_prime[i_zplus] < 0.0)) {
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
  
  int ncz = mesh.LocalNz;

#ifdef BOUT_HAS_Z_GUARD_CELLS_IMPLEMENTED
  throw BoutException("FCIMap/Zoidberg need updating to account for z-guard cells.");
#endif

  BoutReal t_x, t_z;

  Coordinates &coord = *(map_mesh.getCoordinates());

  BOUT_FOR(i, R.getRegion("RGN_NOBNDRY")) {

    const auto x = i.x(), y = i.y(), z = i.z();

    // The integer part of xt_prime, zt_prime are the indices of the cell
    // containing the field line end-point
    i_corner(x, y, z) = static_cast<int>(floor(xt_prime[i]));

    // z is periodic, so make sure the z-index wraps around
    if (zperiodic) {
      zt_prime[i] = zt_prime[i]
                    - ncz * (static_cast<int>(zt_prime[i] / static_cast<BoutReal>(ncz)));

      if (zt_prime[i] < 0.0)
        zt_prime[i] += ncz;
    }

    k_corner(x, y, z) = static_cast<int>(floor(zt_prime[i]));

    // t_x, t_z are the normalised coordinates \in [0,1) within the cell
    // calculated by taking the remainder of the floating point index
    t_x = xt_prime[i] - static_cast<BoutReal>(i_corner(x, y, z));
    t_z = zt_prime[i] - static_cast<BoutReal>(k_corner(x, y, z));

    //----------------------------------------
    // Boundary stuff
    //
    // If a field line leaves the domain, then the forward or backward
    // indices (forward/backward_xt_prime and forward/backward_zt_prime)
    // are set to -1

    if (xt_prime[i] < 0.0) {
      // Hit a boundary
      const auto xp = i.xp(), xm = i.xm();
      const auto zp = i.zp(), zm = i.zm();

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

      BoutReal dR_dx = 0.5 * (R[xp] - R[xm]);
      BoutReal dZ_dx = 0.5 * (Z[xp] - Z[xm]);

      BoutReal dR_dz, dZ_dz;
      // Handle the edge cases in Z
      if (z == 0) {
        dR_dz = R[zp] - R[i];
        dZ_dz = Z[zp] - Z[i];

      } else if (z == mesh.LocalNz - 1) {
        dR_dz = R[i] - R[zm];
        dZ_dz = Z[i] - Z[zm];

      } else {
        dR_dz = 0.5 * (R[zp] - R[zm]);
        dZ_dz = 0.5 * (Z[zp] - Z[zm]);
      }

      BoutReal det = dR_dx * dZ_dz - dR_dz * dZ_dx; // Determinant of 2x2 matrix

      BoutReal dR = R_prime[i] - R[i];
      BoutReal dZ = Z_prime[i] - Z[i];

      // Invert 2x2 matrix to get change in index
      BoutReal dx = (dZ_dz * dR - dR_dz * dZ) / det;
      BoutReal dz = (dR_dx * dZ - dZ_dx * dR) / det;
      boundary->add_point(
          x, y, z, x + dx, y + 0.5 * offset,
          z + dz,               // Intersection point in local index space
          0.5 * coord.dy(x, y), // sqrt( SQ(dR) + SQ(dZ) ),  // Distance to intersection
          PI                    // Right-angle intersection
          );
    }

    //----------------------------------------

    // Check that t_x and t_z are in range
    if ((t_x < 0.0) || (t_x > 1.0))
      throw BoutException("t_x=%e out of range at (%d,%d,%d)", t_x, x, y, z);

    if ((t_z < 0.0) || (t_z > 1.0))
      throw BoutException("t_z=%e out of range at (%d,%d,%d)", t_z, x, y, z);
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

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    const auto x = i.x(), y = i.y(), z = i.z();
    if (boundary_mask(x, y, z))
      continue;

    const auto ynext = i.y() + offset;
    const auto zm = i.zm();
    const auto xm = i.xm();

    BoutReal f_c = centre(x, ynext, z);

    if (corner_boundary_mask(x, y, z) || corner_boundary_mask(xm.ind, y, z)
        || corner_boundary_mask(x, y, zm.ind) || corner_boundary_mask(xm.ind, y, zm.ind)
        || (x == map_mesh.xstart)) {
      // One of the corners leaves the domain.
      // Use the cell centre value, since boundary conditions are not
      // currently applied to corners.
      result(x, ynext, z) = f_c;

    } else {
      BoutReal f_pp = corner(x, ynext, z);      // (x+1/2, z+1/2)
      BoutReal f_mp = corner(x - 1, ynext, z);  // (x-1/2, z+1/2)
      BoutReal f_pm = corner(x, ynext, zm.ind); // (x+1/2, z-1/2)
      BoutReal f_mm = corner(x - 1, ynext, zm.ind); // (x-1/2, z-1/2)

      // This uses a simple weighted average of centre and corners
      // A more sophisticated approach might be to use e.g. Gauss-Lobatto points
      // which would include cell edges and corners
      result(x, ynext, z) = 0.5 * (f_c + 0.25 * (f_pp + f_mp + f_pm + f_mm));

      ASSERT2(finite(result(x, ynext, z)));
    }
  }
  return result;
}

void FCITransform::checkInputGrid() {
  std::string coordinates_type = "";
  if (!mesh.get(coordinates_type, "coordinates_type")) {
    if (coordinates_type != "fci") {
      throw BoutException("Incorrect coordinate system type '"+coordinates_type+"' used "
          "to generate metric components for FCITransform. Should be 'fci'.");
    }
  } // else: coordinate_system variable not found in grid input, indicates older input
    //       file so must rely on the user having ensured the type is correct
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
