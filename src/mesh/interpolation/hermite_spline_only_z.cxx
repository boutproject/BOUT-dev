/**************************************************************************
 * Copyright 2015-2019 B.D.Dudson, P. Hill, J. Omotani, J. Parker
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

HermiteSplineOnlyZ::HermiteSplineOnlyZ(int y_offset, Mesh *mesh)
    : Interpolation(y_offset, mesh),
      h00_z(localmesh), h01_z(localmesh), h10_z(localmesh), h11_z(localmesh) {

  // Index arrays contain guard cells in order to get subscripts right
  k_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);

  // Initialise in order to avoid 'uninitialized value' errors from Valgrind when using
  // guard-cell values
  k_corner = -1;

  // Allocate Field3D members
  h00_z.allocate();
  h01_z.allocate();
  h10_z.allocate();
  h11_z.allocate();
}

void HermiteSplineOnlyZ::calcWeights(const Field3D& delta_x, const Field3D& delta_z) {

  BoutReal t_z;

  for (int x = localmesh->xstart; x <= localmesh->xend; x++) {
    for (int y = localmesh->ystart; y <= localmesh->yend; y++) {
      for (int z = 0; z < localmesh->LocalNz; z++) {

        if (skip_mask(x, y, z))
          continue;

        // The integer part of zt_prime are the indices of the cell
        // containing the field line end-point
        k_corner(x, y, z) = static_cast<int>(floor(delta_z(x, y, z)));

        // t_z is the normalised coordinate \in [0,1) within the cell
        // calculated by taking the remainder of the floating point index
        t_z = delta_z(x, y, z) - static_cast<BoutReal>(k_corner(x, y, z));

        // Check that t_z is in range
        if ((t_z < 0.0) || (t_z > 1.0)) {
          throw BoutException(
              "t_z=%e out of range at (%d,%d,%d) (delta_z=%e, k_corner=%d)", t_z, x, y,
              z, delta_z(x, y, z), k_corner(x, y, z));
        }

        h00_z(x, y, z) = (2. * t_z * t_z * t_z) - (3. * t_z * t_z) + 1.;

        h01_z(x, y, z) = (-2. * t_z * t_z * t_z) + (3. * t_z * t_z);

        h10_z(x, y, z) = t_z * (1. - t_z) * (1. - t_z);

        h11_z(x, y, z) = (t_z * t_z * t_z) - (t_z * t_z);
      }
    }
  }
}

void HermiteSplineOnlyZ::calcWeights(const Field3D &delta_x, const Field3D &delta_z, const BoutMask &mask) {
  skip_mask = mask;
  calcWeights(delta_x, delta_z);
}

/*!
 * Return position and weight of points needed to approximate the function value at the point that the
 * field line through (i,j,k) meets the (j+1)-plane.
 * For the case where only the z-direction is not aligned to grid points, the approximation is:
 *   f(i,j+1,k*) = h00_z * f(i,j+1,k) + h01_z * f(i,j+1,k+1)
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
std::vector<ParallelTransform::positionsAndWeights> HermiteSplineOnlyZ::getWeightsForYApproximation(int i, int j, int k, int yoffset) {
  std::vector<ParallelTransform::positionsAndWeights> pw;
  ParallelTransform::positionsAndWeights p;

  int ncz = localmesh->LocalNz;
  int k_mod = ((k_corner(i, j, k) % ncz) + ncz) % ncz;
  int k_mod_m1 = (k_mod > 0) ? (k_mod-1) : (ncz-1);
  int k_mod_p1 = (k_mod + 1) % ncz;
  int k_mod_p2 = (k_mod + 2) % ncz;

  // Same x, y for all:
  p.i = i;
  p.j = j + yoffset;

  p.k = k_mod_m1;
  p.weight = -0.5*h10_z(i,j,k);
  pw.push_back(p);

  p.k = k_mod;
  p.weight = h00_z(i,j,k) - 0.5*h11_z(i,j,k);
  pw.push_back(p);

  p.k = k_mod_p1;
  p.weight = h01_z(i,j,k) + 0.5*h10_z(i,j,k);
  pw.push_back(p);

  p.k = k_mod_p2;
  p.weight = 0.5*h11_z(i,j,k);
  pw.push_back(p);

  return pw;

}

Field3D HermiteSplineOnlyZ::interpolate(const Field3D &f) const {

  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  Field3D fz = bout::derivatives::index::DDZ(f, CELL_DEFAULT, "DEFAULT", "RGN_ALL");
  localmesh->communicateXZ(fz);

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

        // Interpolate in Z
        f_interp(x, y_next, z) = +f(x, y_next, z_mod) * h00_z(x, y, z) +
                                 f(x, y_next, z_mod_p1) * h01_z(x, y, z) +
                                 fz(x, y_next, z_mod) * h10_z(x, y, z) +
                                 fz(x, y_next, z_mod_p1) * h11_z(x, y, z);

        ASSERT2(finite(f_interp(x, y_next, z)));
      }
    }
  }
  return f_interp;
}

Field3D HermiteSplineOnlyZ::interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z) {
  calcWeights(delta_x, delta_z);
  return interpolate(f);
}

Field3D HermiteSplineOnlyZ::interpolate(const Field3D& f, const Field3D &delta_x, const Field3D &delta_z, const BoutMask &mask) {
  calcWeights(delta_x, delta_z, mask);
  return interpolate(f);
}
