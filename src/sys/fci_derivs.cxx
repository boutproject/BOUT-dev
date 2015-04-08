/**************************************************************************
 * Flux-coordinate Independent parallel derivatives
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
 *
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

#include <fci_boundary_region.hxx>
#include <fci_boundary_op.hxx>
#include <fci_derivs.hxx>
#include <derivs.hxx>
#include <msg_stack.hxx>
#include <bout/mesh.hxx>
#include <bout/assert.hxx>

#include <algorithm>

// Calculate all the coefficients needed for the spline interpolation
// dir MUST be either +1 or -1
FCIMap::FCIMap(Mesh& mesh, int dir, bool yperiodic, bool zperiodic) : dir(dir) {

  // Index arrays contain guard cells in order to get subscripts right
  i_corner = i3tensor(mesh.ngx, mesh.ngy, mesh.ngz-1);
  k_corner = i3tensor(mesh.ngx, mesh.ngy, mesh.ngz-1);

  bool x_boundary;     // has the field line left the domain through the x-sides
  bool y_boundary;     // has the field line left the domain through the y-sides
  bool z_boundary;     // has the field line left the domain through the z-sides

  // Make the boundary_mask the correct size
  // Ugly ugly code
  boundary_mask.resize(mesh.ngx, std::vector<std::vector<bool>>
                       (mesh.ngy, std::vector<bool>
                        (mesh.ngz-1)));

  Field3D xt_prime, zt_prime;

  // Load the floating point indices from the grid file
  // Future, higher order parallel derivatives could require maps to +/-2 slices
  if (dir == +1) {
    mesh.get(xt_prime, "forward_xt_prime");
    mesh.get(zt_prime, "forward_zt_prime");
    boundary = new BoundaryRegionFCI("FCI_forward", BNDRY_FCI_FWD);
  } else if (dir == -1) {
    mesh.get(xt_prime, "backward_xt_prime");
    mesh.get(zt_prime, "backward_zt_prime");
    boundary = new BoundaryRegionFCI("FCI_backward", BNDRY_FCI_BKWD);
  } else {
    // Definitely shouldn't be called
    throw BoutException("FCIMap called with strange direction: %d. Only +/-1 currently supported.", dir);
  }

  // Allocate Field3D members
  y_prime.allocate();
  h00_x.allocate();
  h01_x.allocate();
  h10_x.allocate();
  h11_x.allocate();
  h00_z.allocate();
  h01_z.allocate();
  h10_z.allocate();
  h11_z.allocate();

  int ncz = mesh.ngz-1;
  BoutReal t_x, t_z, temp;

  Coordinates& coord = *(mesh.coordinates());

  for(int x=mesh.xstart;x<=mesh.xend;x++) {
    for(int y=mesh.ystart; y<=mesh.yend;y++) {
      for(int z=0;z<ncz;z++) {

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

        // Distances to intersections with boundaries
        BoutReal y_prime_x;
        BoutReal y_prime_y;
        BoutReal y_prime_z;

        // Field line leaves through x boundary
        if (xt_prime(x,y,z) < 0 ||
            xt_prime(x,y,z) >= mesh.GlobalNx - 1) {
          x_boundary = true;

          BoutReal dx2 = coord.dx(x,y)/2.;
          BoutReal dy = coord.dy(x,y);
          y_prime_x =  dx2 * (dy / (t_x * coord.dx(x, y)));
        } else {
          x_boundary = false;
        }

        // Field line leaves through y boundary
        // Only add this point if the domain is NOT periodic in y
        if ((y + dir < 0 ||
             y + dir > mesh.GlobalNy - 1) && !yperiodic) {
          y_boundary = true;

          y_prime_y =  coord.dy(x,y) / 2.;
        } else {
          y_boundary = false;
        }

        // Field line leaves through z boundary
        // Only add this point if the domain is NOT periodic in Z
        if ((zt_prime(x,y,z) < 0 ||
             zt_prime(x,y,z) > ncz-1) && !zperiodic) {
          z_boundary = true;

          BoutReal dz2 = coord.dz/2.;
          BoutReal dy = coord.dy(x,y);
          y_prime_z =  dz2 * (dy / (t_z * coord.dz));
        } else {
          z_boundary = false;
        }

        // If field line leaves the domain at this point, then add it
        // to the boundary
        if (x_boundary || y_boundary || z_boundary) {
          boundary_mask[x][y][z] = true;
          boundary->add_point(x, y, z);
        }

        // Find the closest intersection with a boundary - seven
        // possible regions field line could end up in
        // Temp variables for convenience
        bool x_b = x_boundary;
        bool y_b = y_boundary;
        bool z_b = z_boundary;
        BoutReal temp;
        if (x_b && !y_b && !z_b) {
          // x
          temp = y_prime_x;
        } else if (!x_b && y_b && !z_b) {
          // y
          temp = y_prime_y;
        } else if (!x_b && !y_b && z_b) {
          // z
          temp = y_prime_z;
        } else if (!x_b && y_b && z_b) {
          // y & z
          temp = std::min(y_prime_y, y_prime_z);
        } else if (x_b && !y_b && z_b) {
          // z & x
          temp = std::min(y_prime_x, y_prime_z);
        } else if (x_b && y_b && !z_b) {
          // x & y
          temp = std::min(y_prime_x, y_prime_y);
        } else if (x_b && y_b && z_b) {
          // x & y & z
          temp = std::min(std::min(y_prime_x, y_prime_y), y_prime_z);
        } else {
          // none
          temp = 0;
        }
        y_prime(x, y, z) = temp;

        //----------------------------------------

        // Check that t_x and t_z are in range
        if( (t_x < 0.0) || (t_x > 1.0) )
          throw BoutException("t_x=%e out of range at (%d,%d,%d)", t_x, x,y,z);

        if( (t_z < 0.0) || (t_z > 1.0) )
          throw BoutException("t_z=%e out of range at (%d,%d,%d)", t_z, x,y,z);

        // NOTE: A (small) hack to avoid one-sided differences
        if( i_corner[x][y][z] == mesh.xend ) {
          i_corner[x][y][z] -= 1;
          t_x = 1.0;
        }

        h00_x(x, y, z) = 2.*t_x*t_x*t_x - 3.*t_x*t_x + 1.;
        h00_z(x, y, z) = 2.*t_z*t_z*t_z - 3.*t_z*t_z + 1.;

        h01_x(x, y, z) = -2.*t_x*t_x*t_x + 3.*t_x*t_x;
        h01_z(x, y, z) = -2.*t_z*t_z*t_z + 3.*t_z*t_z;

        h10_x(x, y, z) = t_x*(1.-t_x)*(1.-t_x);
        h10_z(x, y, z) = t_z*(1.-t_z)*(1.-t_z);

        h11_x(x, y, z) = t_x*t_x*t_x - t_x*t_x;
        h11_z(x, y, z) = t_z*t_z*t_z - t_z*t_z;
      }
    }
  }

}

/**
 * Use cubic Hermite splines to interpolate field f
 *
 * Use cubic Hermite splines to interpolate field f on the adjacent
 * toroidal slice. Spline coefficients and direction of slice are
 * stored in fcimap, and the interpolated field is stored in either
 * f.yup or f.ydown, according to the direction.
 *
 * @param f      The field to interpolate
 * @param fcimap Information on mapping field lines onto next slice
 */
void FCI::interpolate(Field3D &f, const FCIMap &fcimap) {

  if(!mesh.FCI)
    return; // Not using FCI method. Print error / warning?

  Field3D fx, fz, fxz;

  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  fx = mesh.indexDDX(f, CELL_DEFAULT, DIFF_DEFAULT);
  mesh.communicate(fx);
  fz = mesh.indexDDZ(f, CELL_DEFAULT, DIFF_DEFAULT, true);
  mesh.communicate(fz);
  fxz = mesh.indexDDX(fz, CELL_DEFAULT, DIFF_DEFAULT);
  mesh.communicate(fxz);

  Field3D& f_next = f.ynext(fcimap.dir);
  f_next = 0;

  for(int x=mesh.xstart;x<=mesh.xend;x++) {
    for(int y=mesh.ystart; y<=mesh.yend;y++) {
      for(int z=0;z<mesh.ngz-1;z++) {

        // If this field line leaves the domain through the
        // x-boundary, or through the z-boundary and the domain is not
        // periodic, skip it
        if (fcimap.boundary_mask[x][y][z]) continue;

        // Due to lack of guard cells in z-direction, we need to ensure z-index
        // wraps around
        int ncz = mesh.ngz-1;
        int z_mod = ((fcimap.k_corner[x][y][z] % ncz) + ncz) % ncz;
        int z_mod_p1 = (z_mod + 1) % ncz;

        // Interpolate f in X at Z
        BoutReal f_z = f(fcimap.i_corner[x][y][z], y + fcimap.dir, z_mod)*fcimap.h00_x(x,y,z)
          + f(fcimap.i_corner[x][y][z]+1, y + fcimap.dir, z_mod)*fcimap.h01_x(x,y,z)
          + fx( fcimap.i_corner[x][y][z], y + fcimap.dir, z_mod)*fcimap.h10_x(x,y,z)
          + fx( fcimap.i_corner[x][y][z]+1, y + fcimap.dir, z_mod)*fcimap.h11_x(x,y,z);

        // Interpolate f in X at Z+1
        BoutReal f_zp1 = f( fcimap.i_corner[x][y][z], y + fcimap.dir, z_mod_p1)*fcimap.h00_x(x,y,z)
          + f( fcimap.i_corner[x][y][z]+1, y + fcimap.dir, z_mod_p1)*fcimap.h01_x(x,y,z)
          + fx( fcimap.i_corner[x][y][z], y + fcimap.dir, z_mod_p1)*fcimap.h10_x(x,y,z)
          + fx( fcimap.i_corner[x][y][z]+1, y + fcimap.dir, z_mod_p1)*fcimap.h11_x(x,y,z);

        // Interpolate fz in X at Z
        BoutReal fz_z = fz(fcimap.i_corner[x][y][z], y + fcimap.dir, z_mod)*fcimap.h00_x(x,y,z)
          + fz( fcimap.i_corner[x][y][z]+1, y + fcimap.dir, z_mod)*fcimap.h01_x(x,y,z)
          + fxz(fcimap.i_corner[x][y][z], y + fcimap.dir, z_mod)*fcimap.h10_x(x,y,z)
          + fxz(fcimap.i_corner[x][y][z]+1, y + fcimap.dir, z_mod)*fcimap.h11_x(x,y,z);

        // Interpolate fz in X at Z+1
        BoutReal fz_zp1 = fz(fcimap.i_corner[x][y][z], y + fcimap.dir, z_mod_p1)*fcimap.h00_x(x,y,z)
          + fz( fcimap.i_corner[x][y][z]+1, y + fcimap.dir, z_mod_p1)*fcimap.h01_x(x,y,z)
          + fxz(fcimap.i_corner[x][y][z], y + fcimap.dir, z_mod_p1)*fcimap.h10_x(x,y,z)
          + fxz(fcimap.i_corner[x][y][z]+1, y + fcimap.dir, z_mod_p1)*fcimap.h11_x(x,y,z);

        // Interpolate in Z
        f_next(x,y + fcimap.dir,z) =
          + f_z    * fcimap.h00_z(x,y,z)
          + f_zp1  * fcimap.h01_z(x,y,z)
          + fz_z   * fcimap.h10_z(x,y,z)
          + fz_zp1 * fcimap.h11_z(x,y,z);
      }
    }
  }
}

/*******************************************************************************
 * Grad_par
 * The parallel derivative along unperturbed B-field
 *
 * If keep is true, then don't throw away the interpolated field
 *******************************************************************************/
const Field3D FCI::Grad_par(Field3D &f) {

#ifdef CHECK
  int msg_pos = msg_stack.push("FCI::Grad_par( Field3D )");
#endif

  Field3D result;

  result.allocate();

  Field3D &yup = f.yup();
  Field3D &ydown = f.ydown();

  Coordinates *coord = mesh.coordinates();

  for (int x=mesh.xstart;x<=mesh.xend;++x) {
    for (int y=mesh.ystart;y<=mesh.yend;++y) {
      for (int z=0;z<mesh.ngz-1;++z) {
	result(x,y,z) = (yup(x,y+1,z) - ydown(x,y-1,z))/(2*coord->dy(x,y)*sqrt(coord->g_22(x,y)));
      }
    }
  }

#ifdef TRACK
  result.name = "FCI::Grad_par("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}


/*******************************************************************************
 * Grad2_par2
 * second parallel derivative
 *
 * (b dot Grad)(b dot Grad)
 *
 * If keep is true, then don't throw away the interpolated field
 *******************************************************************************/
const Field3D FCI::Grad2_par2(Field3D &f) {

#ifdef CHECK
  int msg_pos = msg_stack.push("FCI::Grad2_par2( Field3D )");
#endif

  Field3D result;

  result.allocate();

  Coordinates *coord = mesh.coordinates();

  Field3D &yup = f.yup();
  Field3D &ydown = f.ydown();

  for (int x=mesh.xstart;x<=mesh.xend;++x) {
	for (int y=mesh.ystart;y<=mesh.yend;++y) {
	  for (int z=0;z<mesh.ngz-1;++z) {
		result(x,y,z) = (yup(x,y+1,z) - 2*f(x,y,z) + ydown(x,y-1,z))/(coord->dy(x,y) * coord->dy(x,y) * coord->g_22(x,y));
	  }
	}
  }
  
#ifdef TRACK
  result.name = "FCI::Grad2_par2("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

/*******************************************************************************
 * Div_par
 * parallel divergence operator B \partial_{||} (F/B)
 *
 * If keep is true, then don't throw away the interpolated field
 *******************************************************************************/
const Field3D FCI::Div_par(Field3D &f) {
#ifdef CHECK
  int msg_pos = msg_stack.push("FCI::Div_par( Field3D )");
#endif

  Coordinates *coord = mesh.coordinates();

  Field3D tmp = f/coord->Bxy;
  Field3D result = coord->Bxy*Grad_par(tmp);

#ifdef TRACK
  result.name = "FCI::Div_par("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

void FCI::applyBoundary(Field3D &f, BndryType bndry_type, FieldGenerator* upvalue, FieldGenerator* downvalue, BoutReal t) {

  BoundaryOpFCI* up_op;
  BoundaryOpFCI* down_op;

  switch(bndry_type) {
  case DIRICHLET:
    up_op = new BoundaryOpFCI_dirichlet(forward_map, upvalue);
    down_op = new BoundaryOpFCI_dirichlet(backward_map, downvalue);
    break;
  case NEUMANN:
    up_op = new BoundaryOpFCI_neumann(forward_map, upvalue);
    down_op = new BoundaryOpFCI_neumann(backward_map, downvalue);
    break;
  default:
    throw BoutException("Not a valid boundary type for FCI!");
  }

  up_op->apply(f, t);
  down_op->apply(f, t);

  delete up_op;
  delete down_op;

}

void FCI::applyBoundary(Field3D &f, BndryType bndry_type, FieldGenerator* upvalue, FieldGenerator* downvalue) {
  applyBoundary(f, bndry_type, upvalue, downvalue, 0);
}

void FCI::applyBoundary(Field3D &f, BndryType bndry_type, FieldGenerator* value, BoutReal t) {
  applyBoundary(f, bndry_type, value, value, t);
}

void FCI::applyBoundary(Field3D &f, BndryType bndry_type, FieldGenerator* value) {
  applyBoundary(f, bndry_type, value, value, 0);
}

void FCI::calcYUpDown(Field3D &f) {

  f.splitYupYdown();

  interpolate(f, forward_map);
  interpolate(f, backward_map);

}
