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

#include <fci_derivs.hxx>
#include <derivs.hxx>
#include <msg_stack.hxx>
#include <bout/mesh.hxx>
#include <bout/assert.hxx>

// Calculate all the coefficients needed for the spline interpolation
// dir MUST be either +1 or -1
FCIMap::FCIMap(Mesh& mesh, int dir) {

  // Index arrays contain guard cells in order to get subscripts right
  i_corner = i3tensor(mesh.ngx, mesh.ngy, mesh.ngz-1);
  k_corner = i3tensor(mesh.ngx, mesh.ngy, mesh.ngz-1);

  Field3D xt_prime, zt_prime;

  // Load the floating point indices from the grid file
  // Future, higher order parallel derivatives could require maps to +/-2 slices
  if (dir == +1) {
	mesh.get(xt_prime, "forward_xt_prime");
	mesh.get(zt_prime, "forward_zt_prime");
  } else if (dir == -1) {
	mesh.get(xt_prime, "backward_xt_prime");
	mesh.get(zt_prime, "backward_zt_prime");
  } else {
	// Definitely shouldn't be called
	throw BoutException("FCIMap called with strange direction: %d. Only +/-1 currently supported.", dir);
  }

  int ncz = mesh.ngz-1;
  BoutReal t_x, t_z, temp;

  for(int x=mesh.xstart;x<=mesh.xend;x++) {
	for(int y=mesh.ystart; y<=mesh.yend;y++) {
	  for(int z=0;z<ncz;z++) {
		// The integer part of xt_prime, zt_prime are the indices of the cell
		// containing the field line end-point
		i_corner[x][y][z] = (int)(xt_prime[x][y][z]);

		// z is periodic, so make sure the z-index wraps around
		zt_prime[x][y][z] = zt_prime[x][y][z] - ncz * ( (int) (zt_prime[x][y][z] / ((BoutReal) ncz)) );

		if(zt_prime[x][y][z] < 0.0)
		  zt_prime[x][y][z] += ncz;

		k_corner[x][y][z] = (int)(zt_prime[x][y][z]);

		// t_x, t_z are the normalised coordinates \in [0,1) within the cell
		// calculated by taking the remainder of the floating point index
		t_x = xt_prime[x][y][z] - (BoutReal)i_corner[x][y][z];
		t_z = zt_prime[x][y][z] - (BoutReal)k_corner[x][y][z];

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

		temp = 2.*t_x*t_x*t_x - 3.*t_x*t_x + 1.;
		h00_x.setData(x, y, z, &temp);
		temp = 2.*t_z*t_z*t_z - 3.*t_z*t_z + 1.;
		h00_z.setData(x, y, z, &temp);

		temp = -2.*t_x*t_x*t_x + 3.*t_x*t_x;
		h01_x.setData(x, y, z, &temp);
		temp = -2.*t_z*t_z*t_z + 3.*t_z*t_z;
		h01_z.setData(x, y, z, &temp);

		temp = t_x*(1.-t_x)*(1.-t_x);
		h10_x.setData(x, y, z, &temp);
		temp = t_z*(1.-t_z)*(1.-t_z);
		h10_z.setData(x, y, z, &temp);

		temp = t_x*t_x*t_x - t_x*t_x;
		h11_x.setData(x, y, z, &temp);
		temp = t_z*t_z*t_z - t_z*t_z;
		h11_z.setData(x, y, z, &temp);
	  }
	}
  }
}

// Use cubic Hermite splines to interpolate field f on the adjacent toroidal
// slice in direction dir. Spline coefficients are stored in fcimap and the
// interpolated field is stored in f_next
void FCI::interpolate(Field3D &f, Field3D &f_next, const FCIMap &fcimap, int dir) {

  if(!mesh.FCI)
	return; // Not using FCI method. Print error / warning?

  Field3D fx, fz, fxz;

  // If f_next has already been computed, don't bother doing it again
  if (f_next.isAllocated())
	return;

  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  fx = DDX(f) * mesh.dx;
  mesh.communicate(fx);
  fz = DDZ(f) * mesh.dz;
  mesh.communicate(fz);
  fxz = D2DXDZ(f) * mesh.dx * mesh.dz;
  mesh.communicate(fxz);

  f_next = 0;

  for(int x=mesh.xstart;x<=mesh.xend;x++) {
	for(int y=mesh.ystart; y<=mesh.yend;y++) {
	  for(int z=0;z<mesh.ngz-1;z++) {

		// Due to lack of guard cells in z-direction, we need to ensure z-index
		// wraps around
		int ncz = mesh.ngz-1;
		int z_mod = ((fcimap.k_corner[x][y][z] % ncz) + ncz) % ncz;
		int z_mod_p1 = (z_mod + 1) % ncz;

		// Interpolate f in X at Z
		BoutReal f_z = f(fcimap.i_corner[x][y][z], y + dir, z_mod)*fcimap.h00_x[x][y][z]
		  + f(fcimap.i_corner[x][y][z]+1, y + dir, z_mod)*fcimap.h01_x[x][y][z]
		  + fx( fcimap.i_corner[x][y][z], y + dir, z_mod)*fcimap.h10_x[x][y][z]
		  + fx( fcimap.i_corner[x][y][z]+1, y + dir, z_mod)*fcimap.h11_x[x][y][z];

		// Interpolate f in X at Z+1
		BoutReal f_zp1 = f( fcimap.i_corner[x][y][z], y + dir, z_mod_p1)*fcimap.h00_x[x][y][z]
		  + f( fcimap.i_corner[x][y][z]+1, y + dir, z_mod_p1)*fcimap.h01_x[x][y][z]
		  + fx( fcimap.i_corner[x][y][z], y + dir, z_mod_p1)*fcimap.h10_x[x][y][z]
		  + fx( fcimap.i_corner[x][y][z]+1, y + dir, z_mod_p1)*fcimap.h11_x[x][y][z];

		// Interpolate fz in X at Z
		BoutReal fz_z = fz(fcimap.i_corner[x][y][z], y + dir, z_mod)*fcimap.h00_x[x][y][z]
		  + fz( fcimap.i_corner[x][y][z]+1, y + dir, z_mod)*fcimap.h01_x[x][y][z]
		  + fxz(fcimap.i_corner[x][y][z], y + dir, z_mod)*fcimap.h10_x[x][y][z]
		  + fxz(fcimap.i_corner[x][y][z]+1, y + dir, z_mod)*fcimap.h11_x[x][y][z];

		// Interpolate fz in X at Z+1
		BoutReal fz_zp1 = fz(fcimap.i_corner[x][y][z], y + dir, z_mod_p1)*fcimap.h00_x[x][y][z]
		  + fz( fcimap.i_corner[x][y][z]+1, y + dir, z_mod_p1)*fcimap.h01_x[x][y][z]
		  + fxz(fcimap.i_corner[x][y][z], y + dir, z_mod_p1)*fcimap.h10_x[x][y][z]
		  + fxz(fcimap.i_corner[x][y][z]+1, y + dir, z_mod_p1)*fcimap.h11_x[x][y][z];

		// Interpolate in Z
		f_next(x,y + dir,z) =
		  + f_z    * fcimap.h00_z[x][y][z]
		  + f_zp1  * fcimap.h01_z[x][y][z]
		  + fz_z   * fcimap.h10_z[x][y][z]
		  + fz_zp1 * fcimap.h11_z[x][y][z];
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
const Field3D FCI::Grad_par(Field3D &f, bool keep) {

#ifdef CHECK
  int msg_pos = msg_stack.push("FCI::Grad_par( Field3D )");
#endif

  Field3D result;
  Field3D *yup, *ydown;

  result.allocate();

  yup = f.yup();
  ydown = f.ydown();

  // Should check if yup, ydown have already been calculated before calling interpolate
  interpolate(f, *yup, forward_map, +1);
  interpolate(f, *ydown, backward_map, -1);

  for (int x=mesh.xstart;x<=mesh.xend;++x) {
	for (int y=mesh.ystart;y<=mesh.yend;++y) {
	  for (int z=0;z<mesh.ngz-1;++z) {
		result(x,y,z) = ((*yup)(x,y+1,z) - (*ydown)(x,y-1,z))/(2*mesh.dy(x,y)*sqrt(mesh.g_22(x,y)));
	  }
	}
  }

  if (!keep) {
	f.resetFCI();
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
const Field3D FCI::Grad2_par2(Field3D &f, bool keep) {

#ifdef CHECK
  int msg_pos = msg_stack.push("FCI::Grad2_par2( Field3D )");
#endif

  Field3D result;
  Field3D *yup, *ydown;

  result.allocate();

  yup = f.yup();
  ydown = f.ydown();

  // Should check if yup, ydown have already been calculated before calling interpolate
  interpolate(f, *yup, forward_map, +1);
  interpolate(f, *ydown, backward_map, -1);

  for (int x=mesh.xstart;x<=mesh.xend;++x) {
	for (int y=mesh.ystart;y<=mesh.yend;++y) {
	  for (int z=0;z<mesh.ngz-1;++z) {
		result(x,y,z) = ((*yup)(x,y+1,z) - 2*f(x,y,z) + (*ydown)(x,y-1,z))/(mesh.dy(x,y) * mesh.dy(x,y) * mesh.g_22(x,y));
	  }
	}
  }

  if (!keep) {
	f.resetFCI();
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
const Field3D FCI::Div_par(Field3D &f, bool keep) {
#ifdef CHECK
  int msg_pos = msg_stack.push("FCI::Div_par( Field3D )");
#endif

  Field3D tmp = f/mesh.Bxy;
  Field3D result = mesh.Bxy*Grad_par(tmp, keep);

#ifdef TRACK
  result.name = "FCI::Div_par("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}
