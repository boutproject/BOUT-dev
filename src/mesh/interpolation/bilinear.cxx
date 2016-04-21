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

#include <string>
#include <vector>

Bilinear::Bilinear(int y_offset) :
  Interpolation(y_offset) {

  w0.allocate();
  w1.allocate();
  w2.allocate();
  w3.allocate();
}

void Bilinear::calcWeights(const Field3D &delta_x, const Field3D &delta_z) {

  for(int x=mesh->xstart;x<=mesh->xend;x++) {
    for(int y=mesh->ystart; y<=mesh->yend;y++) {
      for(int z=0;z<mesh->LocalNz;z++) {

        if (skip_mask(x, y, z)) continue;

        // The integer part of xt_prime, zt_prime are the indices of the cell
        // containing the field line end-point
        int i_corner = floor(delta_x(x,y,z));
        int k_corner = floor(delta_z(x,y,z));

        // t_x, t_z are the normalised coordinates \in [0,1) within the cell
        // calculated by taking the remainder of the floating point index
        BoutReal t_x = delta_x(x,y,z) - static_cast<BoutReal>(i_corner);
        BoutReal t_z = delta_z(x,y,z) - static_cast<BoutReal>(k_corner);
        BoutReal t_x1 = BoutReal(1.0) - t_x;
        BoutReal t_z1 = BoutReal(1.0) - t_z;

        // NOTE: A (small) hack to avoid one-sided differences
        if( i_corner == mesh->xend ) {
          i_corner -= 1;
          t_x = 1.0;
          t_x1 = 0.0;
        }

        // Check that t_x and t_z are in range
        if( (t_x < 0.0) || (t_x > 1.0) )
          throw BoutException("t_x=%e out of range at (%d,%d,%d)", t_x, x,y,z);

        if( (t_z < 0.0) || (t_z > 1.0) )
          throw BoutException("t_z=%e out of range at (%d,%d,%d)", t_z, x,y,z);

        w0(x,y,z) = t_x1 * t_z1;
        w1(x,y,z) = t_x  * t_z1;
        w2(x,y,z) = t_x1 * t_z;
        w3(x,y,z) = t_x  * t_z;

      }
    }
  }
}

void Bilinear::calcWeights(const Field3D &delta_x, const Field3D &delta_z, BoutMask mask) {
  skip_mask = mask;
  calcWeights(delta_x, delta_z);
}

const Field3D Bilinear::interpolate(const Field3D& f) const {

  Field3D f_interp;
  f_interp.allocate();

  for(int x=mesh->xstart;x<=mesh->xend;x++) {
    for(int y=mesh->ystart; y<=mesh->yend;y++) {
      for(int z=0;z<mesh->LocalNz;z++) {

        if (skip_mask(x, y, z)) continue;

        int y_next = y + y_offset;

        f_interp(x,y_next,z) = f(x,y_next,z) * w0(x,y,z) + f(x+1,y_next,z) * w1(x,y,z) + f(x,y_next,z+1) * w2(x,y,z) + f(x+1,y_next,z+1) * w3(x,y,z);

      }
    }
  }
  return f_interp;
}

const Field3D Bilinear::interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z) {
  calcWeights(delta_x, delta_z);
  return interpolate(f);
}

const Field3D Bilinear::interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z, BoutMask mask) {
  calcWeights(delta_x, delta_z, mask);
  return interpolate(f);
}
