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
#include <bout_types.hxx> // See this for codes
#include <msg_stack.hxx>
#include <utils.hxx>

/**
 * Return the sign of val
 */
inline BoutReal sgn(BoutReal val) { return (BoutReal(0) < val) - (val < BoutReal(0)); }

// Calculate all the coefficients needed for the spline interpolation
// dir MUST be either +1 or -1
FCIMap::FCIMap(Mesh &mesh, int dir, bool zperiodic)
  : dir(dir), boundary_mask(mesh), corner_boundary_mask(mesh), y_prime(&mesh) {

  interp = InterpolationFactory::getInstance()->create(&mesh);
  interp->setYOffset(dir);

  interp_corner = InterpolationFactory::getInstance()->create(&mesh);
  interp_corner->setYOffset(dir);
  
  // Index arrays contain guard cells in order to get subscripts right
  // x-index of bottom-left grid point
  auto i_corner = Tensor<int>(mesh.LocalNx, mesh.LocalNy, mesh.LocalNz);
  // z-index of bottom-left grid point
  auto k_corner = Tensor<int>(mesh.LocalNx, mesh.LocalNy, mesh.LocalNz);

  Field3D xt_prime(&mesh), zt_prime(&mesh);
  Field3D R(&mesh), Z(&mesh); // Real-space coordinates of grid points
  Field3D R_prime(&mesh),
      Z_prime(&mesh); // Real-space coordinates of forward/backward points

  mesh.get(R, "R", 0.0, false);
  mesh.get(Z, "Z", 0.0, false);

  // Load the floating point indices from the grid file
  // Future, higher order parallel derivatives could require maps to +/-2 slices
  if (dir == +1) {
    mesh.get(xt_prime, "forward_xt_prime", 0.0, false);
    mesh.get(zt_prime, "forward_zt_prime", 0.0, false);
    mesh.get(R_prime, "forward_R", 0.0, false);
    mesh.get(Z_prime, "forward_Z", 0.0, false);
    boundary = new BoundaryRegionPar("FCI_forward", BNDRY_PAR_FWD, dir, &mesh);
  } else if (dir == -1) {
    mesh.get(xt_prime, "backward_xt_prime", 0.0, false);
    mesh.get(zt_prime, "backward_zt_prime", 0.0, false);
    mesh.get(R_prime, "backward_R", 0.0, false);
    mesh.get(Z_prime, "backward_Z", 0.0, false);
    boundary = new BoundaryRegionPar("FCI_backward", BNDRY_PAR_BKWD, dir, &mesh);
  } else {
    // Definitely shouldn't be called
    throw BoutException("FCIMap called with strange direction: %d. Only +/-1 currently supported.", dir);
  }

  // Add the boundary region to the mesh's vector of parallel boundaries
  mesh.addBoundaryPar(boundary);
  
  // Cell corners
  Field3D xt_prime_corner(&mesh), zt_prime_corner(&mesh);
  xt_prime_corner.allocate();
  zt_prime_corner.allocate();

  BOUT_FOR(i, R.getRegion("RGN_NOBNDRY")) {
    const auto xp = i.xp();
    const auto zp = i.zp();
    const auto xpzp = xp.zp();

    // Point interpolated from (x+1/2, z+1/2)
    if ((xt_prime[i] < 0.0) || (xt_prime[xp] < 0.0) || (xt_prime[xpzp] < 0.0)
        || (xt_prime[zp] < 0.0)) {
      // Hit a boundary
      corner_boundary_mask[i] = true;

      xt_prime_corner[i] = -1.0;
      zt_prime_corner[i] = -1.0;
      continue;
    }

    xt_prime_corner[i] =
        0.25 * (xt_prime[i] + xt_prime[xp] + xt_prime[zp] + xt_prime[xpzp]);

    zt_prime_corner[i] =
        0.25 * (zt_prime[i] + zt_prime[xp] + zt_prime[zp] + zt_prime[xpzp]);
  }

  interp_corner->setMask(corner_boundary_mask);
  interp_corner->calcWeights(xt_prime_corner, zt_prime_corner);
  
  interp->calcWeights(xt_prime, zt_prime);
  
  int ncz = mesh.LocalNz;

#ifdef BOUT_HAS_Z_GUARD_CELLS_IMPLEMENTED
  throw BoutException("FCIMap/Zoidberg need updating to account for z-guard cells.");
#endif

  BoutReal t_x, t_z;

  Coordinates &coord = *(mesh.getCoordinates());

  BOUT_FOR(i, R.getRegion("RGN_NOBNDRY")) {
    // The integer part of xt_prime, zt_prime are the indices of the cell
    // containing the field line end-point
    i_corner[i] = static_cast<int>(floor(xt_prime[i]));

    // z is periodic, so make sure the z-index wraps around
    if (zperiodic) {
      zt_prime[i] = zt_prime[i]
                    - ncz * (static_cast<int>(zt_prime[i] / static_cast<BoutReal>(ncz)));

      if (zt_prime[i] < 0.0)
        zt_prime[i] += ncz;
    }

    k_corner[i] = static_cast<int>(floor(zt_prime[i]));

    // t_x, t_z are the normalised coordinates \in [0,1) within the cell
    // calculated by taking the remainder of the floating point index
    t_x = xt_prime[i] - static_cast<BoutReal>(i_corner[i]);
    t_z = zt_prime[i] - static_cast<BoutReal>(k_corner[i]);

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
      const auto x = i.x(), y = i.y(), z = i.z();

      boundary_mask[i] = true;

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
          x, y, z, x + dx, y + 0.5 * dir,
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
    }
  }

  interp->setMask(boundary_mask);
}

const Field3D FCIMap::integrate(Field3D &f) const {
  TRACE("FCIMap::integrate");
  
  // Cell centre values
  Field3D centre = interp->interpolate(f);
  
  // Cell corner values (x+1/2, z+1/2)
  Field3D corner = interp_corner->interpolate(f);

  Field3D result(f.getMesh());
  result.allocate();
  result.setLocation(f.getLocation());

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    if (boundary_mask[i])
      continue;

    const auto ynext = i.y() + dir;
    const auto x = i.x(), y = i.y(), z = i.z();
    const auto zm = i.zm(), zp = i.zp();
    const auto xm = i.xm(), xp = i.xp(), xmzm = xm.zm();

    BoutReal f_c = centre(x, ynext, z);

    if (corner_boundary_mask[i] || corner_boundary_mask[xm] || corner_boundary_mask[zm]
        || corner_boundary_mask[xmzm] || (x == mesh->xstart)) {
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

      ASSERT2(finite(result(x, ynext, z)));
    }
  }
  return result;
}

void FCITransform::calcYUpDown(Field3D &f) {
  TRACE("FCITransform::calcYUpDown");

  // Ensure that yup and ydown are different fields
  f.splitYupYdown();

  // Interpolate f onto yup and ydown fields
  f.ynext(forward_map.dir) = forward_map.interpolate(f);
  f.ynext(backward_map.dir) = backward_map.interpolate(f);
}

void FCITransform::integrateYUpDown(Field3D &f) {
  TRACE("FCITransform::integrateYUpDown");
  
  // Ensure that yup and ydown are different fields
  f.splitYupYdown();

  // Integrate f onto yup and ydown fields
  f.ynext(forward_map.dir) = forward_map.integrate(f);
  f.ynext(backward_map.dir) = backward_map.integrate(f);
}
