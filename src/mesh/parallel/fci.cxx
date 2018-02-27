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
    : dir(dir), boundary_mask(mesh), y_prime(&mesh) {

  interp = InterpolationFactory::getInstance()->create();
  interp->setYOffset(dir);

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
    boundary = new BoundaryRegionPar("FCI_forward", BNDRY_PAR_FWD, dir);
  } else if (dir == -1) {
    mesh.get(xt_prime, "backward_xt_prime", 0.0, false);
    mesh.get(zt_prime, "backward_zt_prime", 0.0, false);
    mesh.get(R_prime, "backward_R", 0.0, false);
    mesh.get(Z_prime, "backward_Z", 0.0, false);
    boundary = new BoundaryRegionPar("FCI_backward", BNDRY_PAR_BKWD, dir);
  } else {
    // Definitely shouldn't be called
    throw BoutException("FCIMap called with strange direction: %d. Only +/-1 currently supported.", dir);
  }

  // Add the boundary region to the mesh's vector of parallel boundaries
  mesh.addBoundaryPar(boundary);

  interp->calcWeights(xt_prime, zt_prime);

  int ncz = mesh.LocalNz;
  BoutReal t_x, t_z;

  Coordinates &coord = *(mesh.coordinates());

  for (int x = mesh.xstart; x <= mesh.xend; x++) {
    for (int y = mesh.ystart; y <= mesh.yend; y++) {
      for (int z = 0; z < ncz; z++) {

        // The integer part of xt_prime, zt_prime are the indices of the cell
        // containing the field line end-point
        i_corner(x, y, z) = static_cast<int>(floor(xt_prime(x, y, z)));

        // z is periodic, so make sure the z-index wraps around
        if (zperiodic) {
          zt_prime(x, y, z) =
              zt_prime(x, y, z) -
              ncz * (static_cast<int>(zt_prime(x, y, z) / static_cast<BoutReal>(ncz)));

          if (zt_prime(x, y, z) < 0.0)
            zt_prime(x, y, z) += ncz;
        }

        k_corner(x, y, z) = static_cast<int>(floor(zt_prime(x, y, z)));

        // t_x, t_z are the normalised coordinates \in [0,1) within the cell
        // calculated by taking the remainder of the floating point index
        t_x = xt_prime(x, y, z) - static_cast<BoutReal>(i_corner(x, y, z));
        t_z = zt_prime(x, y, z) - static_cast<BoutReal>(k_corner(x, y, z));

        //----------------------------------------
        // Boundary stuff
        //
        // If a field line leaves the domain, then the forward or backward
        // indices (forward/backward_xt_prime and forward/backward_zt_prime)
        // are set to -1

        if (xt_prime(x, y, z) < 0.0) {
          // Hit a boundary

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

          BoutReal dR_dx = 0.5 * (R(x + 1, y, z) - R(x - 1, y, z));
          BoutReal dZ_dx = 0.5 * (Z(x + 1, y, z) - Z(x - 1, y, z));

          BoutReal dR_dz, dZ_dz;
          // Handle the edge cases in Z
          if (z == 0) {
            dR_dz = R(x, y, z + 1) - R(x, y, z);
            dZ_dz = Z(x, y, z + 1) - Z(x, y, z);

          } else if (z == mesh.LocalNz - 1) {
            dR_dz = R(x, y, z) - R(x, y, z - 1);
            dZ_dz = Z(x, y, z) - Z(x, y, z - 1);

          } else {
            dR_dz = 0.5 * (R(x, y, z + 1) - R(x, y, z - 1));
            dZ_dz = 0.5 * (Z(x, y, z + 1) - Z(x, y, z - 1));
          }

          BoutReal det = dR_dx * dZ_dz - dR_dz * dZ_dx; // Determinant of 2x2 matrix

          BoutReal dR = R_prime(x, y, z) - R(x, y, z);
          BoutReal dZ = Z_prime(x, y, z) - Z(x, y, z);

          // Invert 2x2 matrix to get change in index
          BoutReal dx = (dZ_dz * dR - dR_dz * dZ) / det;
          BoutReal dz = (dR_dx * dZ - dZ_dx * dR) / det;
          boundary->add_point(x, y, z, 
                              x + dx, y + 0.5*dir, z + dz,  // Intersection point in local index space
                              0.5*coord.dy(x,y), //sqrt( SQ(dR) + SQ(dZ) ),  // Distance to intersection
                              PI   // Right-angle intersection
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

void FCITransform::calcYUpDown(Field3D &f) {
  TRACE("FCITransform::calcYUpDown");

  // Ensure that yup and ydown are different fields
  f.splitYupYdown();

  // Interpolate f onto yup and ydown fields
  f.ynext(forward_map.dir) = forward_map.interpolate(f);
  f.ynext(backward_map.dir) = backward_map.interpolate(f);
}
