/**************************************************************************
 * Implements the Flux Coordinate Independent scheme for parallel derivatives

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
#include <bout/mesh.hxx>
#include <bout_types.hxx> // See this for codes
#include <msg_stack.hxx>
#include <utils.hxx>

/**
 * Return the sign of val
 */
inline BoutReal sgn(BoutReal val) {
    return (BoutReal(0) < val) - (val < BoutReal(0));
}

// Calculate all the coefficients needed for the spline interpolation
// dir MUST be either +1 or -1
FCIMap::FCIMap(Mesh& mesh, int dir, bool yperiodic, bool zperiodic) :
  dir(dir), boundary_mask(mesh) {

  interp = InterpolationFactory::getInstance()->create();
  interp->setYOffset(dir);

  // Index arrays contain guard cells in order to get subscripts right
  // x-index of bottom-left grid point
  int*** i_corner = i3tensor(mesh.LocalNx, mesh.LocalNy, mesh.LocalNz);
  // z-index of bottom-left grid point
  int*** k_corner = i3tensor(mesh.LocalNx, mesh.LocalNy, mesh.LocalNz);

  bool x_boundary;     // has the field line left the domain through the x-sides
  bool y_boundary;     // has the field line left the domain through the y-sides
  bool z_boundary;     // has the field line left the domain through the z-sides

  Field3D xt_prime, zt_prime;

  // Load the floating point indices from the grid file
  // Future, higher order parallel derivatives could require maps to +/-2 slices
  if (dir == +1) {
    mesh.get(xt_prime, "forward_xt_prime", 0.0, false);
    mesh.get(zt_prime, "forward_zt_prime", 0.0, false);
    boundary = new BoundaryRegionPar("FCI_forward", BNDRY_PAR_FWD, dir);
  } else if (dir == -1) {
    mesh.get(xt_prime, "backward_xt_prime", 0.0, false);
    mesh.get(zt_prime, "backward_zt_prime", 0.0, false);
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

  Coordinates& coord = *(mesh.coordinates());

  // Vector in real space
  struct RealVector {
    BoutReal x;
    BoutReal y;
    BoutReal z;
  };

  for(int x=mesh.xstart; x<=mesh.xend; x++) {
    for(int y=mesh.ystart;  y<=mesh.yend; y++) {
      for(int z=0; z<ncz; z++) {

        // Dot product of two vectors
        // Only needed in this function, so use a named lambda
        // Defined inside loop to capture x, y, z
        auto dot = [&](const RealVector &lhs, const RealVector &rhs) {
          BoutReal result;
          result = lhs.x*rhs.x*coord.g11(x, y)
          + lhs.y*rhs.y*coord.g22(x, y)
          + lhs.z*rhs.z*coord.g33(x, y);
          result += (lhs.x*rhs.y + lhs.y*rhs.x)*coord.g12(x, y)
          + (lhs.x*rhs.z + lhs.z*rhs.x)*coord.g13(x, y)
          + (lhs.y*rhs.z + lhs.z*rhs.y)*coord.g23(x, y);

          return result;
        };

        // The integer part of xt_prime, zt_prime are the indices of the cell
        // containing the field line end-point
        i_corner[x][y][z] = floor(xt_prime(x,y,z));

        // z is periodic, so make sure the z-index wraps around
        if (zperiodic) {
          zt_prime(x,y,z) = zt_prime(x,y,z) - ncz * ( (int)(zt_prime(x,y,z) / ((BoutReal) ncz)) );

          if (zt_prime(x,y,z) < 0.0)
            zt_prime(x,y,z) += ncz;
        }

        k_corner[x][y][z] = floor(zt_prime(x,y,z));

        // t_x, t_z are the normalised coordinates \in [0,1) within the cell
        // calculated by taking the remainder of the floating point index
        t_x = xt_prime(x,y,z) - (BoutReal)i_corner[x][y][z];
        t_z = zt_prime(x,y,z) - (BoutReal)k_corner[x][y][z];

        //----------------------------------------
        // Boundary stuff

        // Field line vector
        RealVector b_hat = {t_x*coord.dx(x,y), coord.dy(x,y), t_z*coord.dz};
        // Length of field line
        BoutReal length = sqrt(dot(b_hat, b_hat));

        // Parameterised distance to intersection point
        BoutReal s_intersect;

        // Parameterised distance to intersection with boundary
        // for the three different boundaries
        BoutReal s_intersect_x;
        BoutReal s_intersect_y;
        BoutReal s_intersect_z;

        // Total (signed) distance (in index space) from this point
        BoutReal p_x = xt_prime(x, y, z) - x;
        BoutReal p_z = zt_prime(x, y, z) - z;

        // Field line leaves through x boundary
        bool x_lower = (mesh.firstX() && (xt_prime(x,y,z) < mesh.xstart));
        bool x_upper = (mesh.lastX()  && (xt_prime(x,y,z) > mesh.xend));
        if (x_lower || x_upper) {
          x_boundary = true;
          // Total distance (in index space) from the boundary
          BoutReal boundary_dist;
          if (x_lower) {
            boundary_dist = abs((mesh.xstart - 0.5) - x);
          } else {
            boundary_dist = abs((mesh.xend + 0.5) - x);
          }
          s_intersect_x = boundary_dist / (abs(p_x));
       } else {
          x_boundary = false;
        }

        // Field line leaves through y boundary
        // Only add this point if the domain is NOT periodic in y
        bool y_lower = !yperiodic && (mesh.firstY() && (y + dir <= mesh.ystart - 0.5));
        bool y_upper = !yperiodic && (mesh.lastY()  && (y + dir >= mesh.yend + 0.5));
        if (y_lower || y_upper) {
          y_boundary = true;
          s_intersect_y = 0.5;
        } else {
          y_boundary = false;
        }

        // Field line leaves through z boundary
        // Only add this point if the domain is NOT periodic in Z
        bool z_lower = !zperiodic && (zt_prime(x,y,z) < 0.0);
        bool z_upper = !zperiodic && (zt_prime(x,y,z) > ncz - 1);
        if (z_lower || z_upper) {
          z_boundary = true;
          // Total distance (in index space) from the boundary
          BoutReal boundary_dist;
          if (z_lower) {
            boundary_dist = abs(-0.5 - z);
          } else {
            boundary_dist = abs((ncz - 0.5) - z);
          }
          s_intersect_z = boundary_dist / (abs(p_z));
        } else {
          z_boundary = false;
        }

        // Find the closest intersection with a boundary - seven
        // possible regions field line could end up in
        if (x_boundary && !y_boundary && !z_boundary) {
          // x
          s_intersect = s_intersect_x;
        } else if (!x_boundary && y_boundary && !z_boundary) {
          // y
          s_intersect = s_intersect_y;
        } else if (!x_boundary && !y_boundary && z_boundary) {
          // z
          s_intersect = s_intersect_z;
        } else if (!x_boundary && y_boundary && z_boundary) {
          // y & z
          s_intersect = std::min(s_intersect_y, s_intersect_z);
        } else if (x_boundary && !y_boundary && z_boundary) {
          // z & x
          s_intersect = std::min(s_intersect_x, s_intersect_z);
        } else if (x_boundary && y_boundary && !z_boundary) {
          // x & y
          s_intersect = std::min(s_intersect_x, s_intersect_y);
        } else if (x_boundary && y_boundary && z_boundary) {
          // x & y & z
          s_intersect = std::min(std::min(s_intersect_x, s_intersect_y), s_intersect_z);
        } else {
          // none
          s_intersect = 0;
        }

        // If field line leaves the domain at this point, then add it
        // to the boundary
        if (x_boundary || y_boundary || z_boundary) {
          // Normal to boundary
          RealVector norm;
          // s_intersect is set to that of the closest boundary, so set the normal
          // based on that
          if (s_intersect == s_intersect_x) {
            norm = {sgn(p_x), 0., 0.};
          } else if (s_intersect == s_intersect_y) {
            norm = {0., static_cast<BoutReal>(dir), 0.};
          } else if (s_intersect == s_intersect_z) {
            norm = {0., 0., sgn(p_z)};
          } else {
            // Shouldn't reach here - boundary set, but s_intersect not
            // equal to a boundary s_intersect_(x|y|z)
            throw BoutException("Something weird happened in FCIMap...");
          }

          // y-distance to boundary intersection
          BoutReal y_prime = coord.dy(x,y) * s_intersect;

          // Index-space coordinates of intersection
          BoutReal s_x = x + s_intersect*p_x;
          BoutReal s_y = y + s_intersect*dir;
          BoutReal s_z = z + s_intersect*p_z;

          // Angle between field line and boundary
          BoutReal angle = asin( dot(norm, b_hat) / length );

          // This would be correct, but need to work out correct modification to
          // the boundary conditions to use it
          // y_prime = s_intersect * length;

          boundary_mask(x, y, z) = true;
          boundary->add_point(x, y, z, s_x, s_y, s_z, y_prime, angle);
        }

        //----------------------------------------

        // Check that t_x and t_z are in range
        if( (t_x < 0.0) || (t_x > 1.0) )
          throw BoutException("t_x=%e out of range at (%d,%d,%d)", t_x, x,y,z);

        if( (t_z < 0.0) || (t_z > 1.0) )
          throw BoutException("t_z=%e out of range at (%d,%d,%d)", t_z, x,y,z);

      }
    }
  }

  interp->setMask(boundary_mask);

  free_i3tensor(i_corner);
  free_i3tensor(k_corner);
}


void FCITransform::calcYUpDown(Field3D &f) {
  TRACE("FCITransform::calcYUpDown");

  // Ensure that yup and ydown are different fields
  f.splitYupYdown();

  // Interpolate f onto yup and ydown fields
  f.ynext(forward_map.dir) = forward_map.interpolate(f);
  f.ynext(backward_map.dir) = backward_map.interpolate(f);
}
