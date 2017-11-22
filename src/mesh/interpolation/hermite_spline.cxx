/**************************************************************************
 * Copyright 2015 B.D.Dudson, P. Hill
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

#include "bout/mesh.hxx"
#include "globals.hxx"
#include "interpolation.hxx"

#include <vector>

HermiteSpline::HermiteSpline(int y_offset) :
  Interpolation(y_offset), msh(nullptr), 
  h00_x(msh), h01_x(msh),h10_x(msh),h11_x(msh),
  h00_z(msh), h01_z(msh),h10_z(msh),h11_z(msh)
{

  // Index arrays contain guard cells in order to get subscripts right
  i_corner = i3tensor(mesh->LocalNx, mesh->LocalNy, mesh->LocalNz);
  k_corner = i3tensor(mesh->LocalNx, mesh->LocalNy, mesh->LocalNz);

  // Allocate Field3D members
  h00_x.allocate();
  h01_x.allocate();
  h10_x.allocate();
  h11_x.allocate();
  h00_z.allocate();
  h01_z.allocate();
  h10_z.allocate();
  h11_z.allocate();
}

void HermiteSpline::calcWeights(const Field3D &delta_x, const Field3D &delta_z) {

  BoutReal t_x, t_z;

  for(int x=mesh->xstart;x<=mesh->xend;x++) {
    for(int y=mesh->ystart; y<=mesh->yend;y++) {
      for(int z=0;z<mesh->LocalNz;z++) {

        if (skip_mask(x, y, z)) continue;

        // The integer part of xt_prime, zt_prime are the indices of the cell
        // containing the field line end-point
        i_corner[x][y][z] = static_cast<int>(floor(delta_x(x,y,z)));
        k_corner[x][y][z] = static_cast<int>(floor(delta_z(x,y,z)));

        // t_x, t_z are the normalised coordinates \in [0,1) within the cell
        // calculated by taking the remainder of the floating point index
        t_x = delta_x(x,y,z) - static_cast<BoutReal>(i_corner[x][y][z]);
        t_z = delta_z(x,y,z) - static_cast<BoutReal>(k_corner[x][y][z]);

        // NOTE: A (small) hack to avoid one-sided differences
        if( i_corner[x][y][z] >= mesh->xend ) {
          i_corner[x][y][z] = mesh->xend-1;
          t_x = 1.0;
        }

        // Check that t_x and t_z are in range
        if( (t_x < 0.0) || (t_x > 1.0) )
          throw BoutException("t_x=%e out of range at (%d,%d,%d)", t_x, x,y,z);

        if( (t_z < 0.0) || (t_z > 1.0) )
          throw BoutException("t_z=%e out of range at (%d,%d,%d)", t_z, x,y,z);

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

void HermiteSpline::calcWeights(const Field3D &delta_x, const Field3D &delta_z, BoutMask mask) {
  skip_mask = mask;
  calcWeights(delta_x, delta_z);
}

Field3D HermiteSpline::interpolate(const Field3D& f) const {

  Field3D f_interp(f.getMesh());
  f_interp.allocate();

  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  Field3D fx = mesh->indexDDX(f, CELL_DEFAULT, DIFF_DEFAULT);
  mesh->communicateXZ(fx);
  Field3D fz = mesh->indexDDZ(f, CELL_DEFAULT, DIFF_DEFAULT, true);
  mesh->communicateXZ(fz);
  Field3D fxz = mesh->indexDDX(fz, CELL_DEFAULT, DIFF_DEFAULT);
  mesh->communicateXZ(fxz);

  for(int x=mesh->xstart;x<=mesh->xend;x++) {
    for(int y=mesh->ystart; y<=mesh->yend;y++) {
      for(int z=0;z<mesh->LocalNz;z++) {

        if (skip_mask(x, y, z)) continue;

        // Due to lack of guard cells in z-direction, we need to ensure z-index
        // wraps around
        int ncz = mesh->LocalNz;
        int z_mod = ((k_corner[x][y][z] % ncz) + ncz) % ncz;
        int z_mod_p1 = (z_mod + 1) % ncz;

        int y_next = y + y_offset;

        // Interpolate f in X at Z
        BoutReal f_z = f( i_corner[x][y][z],   y_next, z_mod)*h00_x(x,y,z)
          + f( i_corner[x][y][z]+1, y_next, z_mod)*h01_x(x,y,z)
          + fx(i_corner[x][y][z],   y_next, z_mod)*h10_x(x,y,z)
          + fx(i_corner[x][y][z]+1, y_next, z_mod)*h11_x(x,y,z);

        // Interpolate f in X at Z+1
        BoutReal f_zp1 = f( i_corner[x][y][z],   y_next, z_mod_p1)*h00_x(x,y,z)
          + f( i_corner[x][y][z]+1, y_next, z_mod_p1)*h01_x(x,y,z)
          + fx(i_corner[x][y][z],   y_next, z_mod_p1)*h10_x(x,y,z)
          + fx(i_corner[x][y][z]+1, y_next, z_mod_p1)*h11_x(x,y,z);

        // Interpolate fz in X at Z
        BoutReal fz_z = fz( i_corner[x][y][z],   y_next, z_mod)*h00_x(x,y,z)
          + fz( i_corner[x][y][z]+1, y_next, z_mod)*h01_x(x,y,z)
          + fxz(i_corner[x][y][z],   y_next, z_mod)*h10_x(x,y,z)
          + fxz(i_corner[x][y][z]+1, y_next, z_mod)*h11_x(x,y,z);

        // Interpolate fz in X at Z+1
        BoutReal fz_zp1 = fz( i_corner[x][y][z],   y_next, z_mod_p1)*h00_x(x,y,z)
          + fz( i_corner[x][y][z]+1, y_next, z_mod_p1)*h01_x(x,y,z)
          + fxz(i_corner[x][y][z],   y_next, z_mod_p1)*h10_x(x,y,z)
          + fxz(i_corner[x][y][z]+1, y_next, z_mod_p1)*h11_x(x,y,z);

        // Interpolate in Z
        f_interp(x,y_next,z) =
          + f_z    * h00_z(x,y,z)
          + f_zp1  * h01_z(x,y,z)
          + fz_z   * h10_z(x,y,z)
          + fz_zp1 * h11_z(x,y,z);
      }
    }
  }
  return f_interp;
}

Field3D HermiteSpline::interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z) {
  calcWeights(delta_x, delta_z);
  return interpolate(f);
}

Field3D HermiteSpline::interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z, BoutMask mask) {
  calcWeights(delta_x, delta_z, mask);
  return interpolate(f);
}
