/**************************************************************************
 * Copyright 2015-2018 B.D.Dudson, P. Hill
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

#include "globals.hxx"
#include "interpolation.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/mesh.hxx"

#include <vector>

HermiteSpline::HermiteSpline(int y_offset, Mesh *mesh)
    : Interpolation(y_offset, mesh),
      h00_x(localmesh), h01_x(localmesh), h10_x(localmesh), h11_x(localmesh),
      h00_z(localmesh), h01_z(localmesh), h10_z(localmesh), h11_z(localmesh) {

  // Index arrays contain guard cells in order to get subscripts right
  i_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);
  k_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);

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

  for (int x = localmesh->xstart; x <= localmesh->xend; x++) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      for (int z = 0; z < localmesh->LocalNz; z++) {

        if (skip_mask(x, y, z))
          continue;

        // The integer part of xt_prime, zt_prime are the indices of the cell
        // containing the field line end-point
        i_corner(x, y, z) = static_cast<int>(floor(delta_x(x, y, z)));
        k_corner(x, y, z) = static_cast<int>(floor(delta_z(x, y, z)));

        // t_x, t_z are the normalised coordinates \in [0,1) within the cell
        // calculated by taking the remainder of the floating point index
        t_x = delta_x(x, y, z) - static_cast<BoutReal>(i_corner(x, y, z));
        t_z = delta_z(x, y, z) - static_cast<BoutReal>(k_corner(x, y, z));

        // NOTE: A (small) hack to avoid one-sided differences
        if (i_corner(x, y, z) >= localmesh->xend) {
          i_corner(x, y, z) = localmesh->xend - 1;
          t_x = 1.0;
        }
        if (i_corner(x, y, z) < localmesh->xstart) {
          i_corner(x, y, z) = localmesh->xstart;
          t_x = 0.0;
        }

        // Check that t_x and t_z are in range
        if ((t_x < 0.0) || (t_x > 1.0)) {
          throw BoutException(
              "t_x={:e} out of range at ({:d},{:d},{:d}) (delta_x={:e}, i_corner={:d})", t_x, x, y,
              z, delta_x(x, y, z), i_corner(x, y, z));
        }

        if ((t_z < 0.0) || (t_z > 1.0)) {
          throw BoutException(
              "t_z={:e} out of range at ({:d},{:d},{:d}) (delta_z={:e}, k_corner={:d})", t_z, x, y,
              z, delta_z(x, y, z), k_corner(x, y, z));
        }

        h00_x(x, y, z) = (2. * t_x * t_x * t_x) - (3. * t_x * t_x) + 1.;
        h00_z(x, y, z) = (2. * t_z * t_z * t_z) - (3. * t_z * t_z) + 1.;

        h01_x(x, y, z) = (-2. * t_x * t_x * t_x) + (3. * t_x * t_x);
        h01_z(x, y, z) = (-2. * t_z * t_z * t_z) + (3. * t_z * t_z);

        h10_x(x, y, z) = t_x * (1. - t_x) * (1. - t_x);
        h10_z(x, y, z) = t_z * (1. - t_z) * (1. - t_z);

        h11_x(x, y, z) = (t_x * t_x * t_x) - (t_x * t_x);
        h11_z(x, y, z) = (t_z * t_z * t_z) - (t_z * t_z);
      }
    }
  }
}

void HermiteSpline::calcWeights(const Field3D &delta_x, const Field3D &delta_z, const BoutMask &mask) {
  skip_mask = mask;
  calcWeights(delta_x, delta_z);
}

/*!
 * Return position and weight of points needed to approximate the function value at the
 * point that the field line through (i,j,k) meets the (j+1)-plane. For the case where
 * only the z-direction is not aligned to grid points, the approximation is: f(i,j+1,k*) =
 * h00_z * f(i,j+1,k) + h01_z * f(i,j+1,k+1)
 *                 + h10_z * dfdz(i,j+1,k) + h11_z * dfdz(i,j+1,k+1)
 *               = h00_z * f(i,j+1,k) + h01_z * f(i,j+1,k+1)
 *                 + h10_z * ( f(i,j+1,k+1) - f(i,j+1,k-1) ) / 2
 *                 + h11_z * ( f(i,j+1,k+2) - f(i,j+1,k) ) / 2
 * for k* a point between k and k+1. Therefore, this function returns
 *   position 		weight
 *   (i, j+1, k-1)	- h10_z / 2
 *   (i, j+1, k)	h00_z - h11_z / 2
 *   (i, j+1, k+1)	h01_z + h10_z / 2
 *   (i, j+1, k+2)	h11_z / 2
 */
std::vector<ParallelTransform::PositionsAndWeights>
HermiteSpline::getWeightsForYApproximation(int i, int j, int k, int yoffset) {
  std::vector<ParallelTransform::PositionsAndWeights> pw;
  ParallelTransform::PositionsAndWeights p;

  int ncz = localmesh->LocalNz;
  int k_mod = ((k_corner(i, j, k) % ncz) + ncz) % ncz;
  int k_mod_m1 = (k_mod > 0) ? (k_mod - 1) : (ncz - 1);
  int k_mod_p1 = (k_mod + 1) % ncz;
  int k_mod_p2 = (k_mod + 2) % ncz;

  // Same x, y for all:
  p.i = i;
  p.j = j + yoffset;

  p.k = k_mod_m1;
  p.weight = -0.5 * h10_z(i, j, k);
  pw.push_back(p);

  p.k = k_mod;
  p.weight = h00_z(i, j, k) - 0.5 * h11_z(i, j, k);
  pw.push_back(p);

  p.k = k_mod_p1;
  p.weight = h01_z(i, j, k) + 0.5 * h10_z(i, j, k);
  pw.push_back(p);

  p.k = k_mod_p2;
  p.weight = 0.5 * h11_z(i, j, k);
  pw.push_back(p);

  return pw;
}

Field3D HermiteSpline::interpolate(const Field3D &f) const {

  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  Field3D fx = bout::derivatives::index::DDX(f, CELL_DEFAULT, "DEFAULT");
  localmesh->communicateXZ(fx);
  Field3D fz = bout::derivatives::index::DDZ(f, CELL_DEFAULT, "DEFAULT", "RGN_ALL");
  localmesh->communicateXZ(fz);
  Field3D fxz = bout::derivatives::index::DDX(fz, CELL_DEFAULT, "DEFAULT");
  localmesh->communicateXZ(fxz);

  for (int x = localmesh->xstart; x <= localmesh->xend; x++) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      for (int z = 0; z < localmesh->LocalNz; z++) {

        if (skip_mask(x, y, z))
          continue;

        // Due to lack of guard cells in z-direction, we need to ensure z-index
        // wraps around
        int ncz = localmesh->LocalNz;
        int z_mod = ((k_corner(x, y, z) % ncz) + ncz) % ncz;
        int z_mod_p1 = (z_mod + 1) % ncz;

        int y_next = y + y_offset;

        // Interpolate f in X at Z
        BoutReal f_z = f(i_corner(x, y, z), y_next, z_mod) * h00_x(x, y, z) +
                       f(i_corner(x, y, z) + 1, y_next, z_mod) * h01_x(x, y, z) +
                       fx(i_corner(x, y, z), y_next, z_mod) * h10_x(x, y, z) +
                       fx(i_corner(x, y, z) + 1, y_next, z_mod) * h11_x(x, y, z);

        // Interpolate f in X at Z+1
        BoutReal f_zp1 = f(i_corner(x, y, z), y_next, z_mod_p1) * h00_x(x, y, z) +
                         f(i_corner(x, y, z) + 1, y_next, z_mod_p1) * h01_x(x, y, z) +
                         fx(i_corner(x, y, z), y_next, z_mod_p1) * h10_x(x, y, z) +
                         fx(i_corner(x, y, z) + 1, y_next, z_mod_p1) * h11_x(x, y, z);

        // Interpolate fz in X at Z
        BoutReal fz_z = fz(i_corner(x, y, z), y_next, z_mod) * h00_x(x, y, z) +
                        fz(i_corner(x, y, z) + 1, y_next, z_mod) * h01_x(x, y, z) +
                        fxz(i_corner(x, y, z), y_next, z_mod) * h10_x(x, y, z) +
                        fxz(i_corner(x, y, z) + 1, y_next, z_mod) * h11_x(x, y, z);

        // Interpolate fz in X at Z+1
        BoutReal fz_zp1 = fz(i_corner(x, y, z), y_next, z_mod_p1) * h00_x(x, y, z) +
                          fz(i_corner(x, y, z) + 1, y_next, z_mod_p1) * h01_x(x, y, z) +
                          fxz(i_corner(x, y, z), y_next, z_mod_p1) * h10_x(x, y, z) +
                          fxz(i_corner(x, y, z) + 1, y_next, z_mod_p1) * h11_x(x, y, z);

        // Interpolate in Z
        f_interp(x, y_next, z) = +f_z * h00_z(x, y, z) + f_zp1 * h01_z(x, y, z) +
                                 fz_z * h10_z(x, y, z) + fz_zp1 * h11_z(x, y, z);

        ASSERT2(finite(f_interp(x, y_next, z)));
      }
    }
  }
  return f_interp;
}

Field3D HermiteSpline::interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z) {
  calcWeights(delta_x, delta_z);
  return interpolate(f);
}

Field3D HermiteSpline::interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z, const BoutMask &mask) {
  calcWeights(delta_x, delta_z, mask);
  return interpolate(f);
}
