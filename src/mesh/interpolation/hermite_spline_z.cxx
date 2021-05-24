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
#include "interpolation_z.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/mesh.hxx"

#include <vector>

ZHermiteSpline::ZHermiteSpline(int y_offset, Mesh* mesh, Region<Ind3D> region_in)
    : ZInterpolation(y_offset, mesh, region_in),
      fz_region(region_in.size() != 0 ? "RGN_ALL"
                                      : y_offset == 0 ? "RGN_NOBNDRY" : "RGN_NOX"),
      h00(localmesh), h01(localmesh), h10(localmesh), h11(localmesh) {

  if (region_in.size() != 0) {
    output_warn << "Custom region passed to ZInterpolation. ZHermiteSpline requires an "
                << "offset region for calculating DDZ(f): using RGN_ALL, which may not "
                << "be optimal." << endl;
  }

  // Index arrays contain guard cells in order to get subscripts right
  const int n_total = localmesh->LocalNx*localmesh->LocalNy*localmesh->LocalNz;
  k_corner.reallocate(n_total);

  // Initialise in order to avoid 'uninitialized value' errors from Valgrind when using
  // guard-cell values
  std::fill(std::begin(k_corner), std::end(k_corner),
            Ind3D(-1, localmesh->LocalNy, localmesh->LocalNz));

  // Allocate Field3D members
  h00.allocate();
  h01.allocate();
  h10.allocate();
  h11.allocate();
}

void ZHermiteSpline::calcWeights(const Field3D& delta_z) {

  const int ncy = localmesh->LocalNy;
  const int ncz = localmesh->LocalNz;

  // Calculate weights for all points if y_offset==0 in case they are needed, otherwise
  // only calculate weights for RGN_NOY, which should be a superset of 'region'
  const auto& local_region = (y_offset == 0) ? delta_z.getRegion("RGN_ALL") : delta_z.getRegion("RGN_NOY");

  BOUT_FOR(i, local_region) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    // The integer part of zt_prime are the indices of the cell
    // containing the field line end-point
    int corner_zind = floor(delta_z(x, y, z));

    // t_z is the normalised coordinate \in [0,1) within the cell
    // calculated by taking the remainder of the floating point index
    const BoutReal t_z = delta_z(x, y, z) - static_cast<BoutReal>(corner_zind);

    // Make corner_zind be in the range 0<=corner_zind<nz
    // This needs to be done after calculating t_z (the coordinate within the
    // cell) because delta_z is allowed to be less than 0 or greater than ncz.
    corner_zind = ((corner_zind % ncz) + ncz) % ncz;

    // Convert z-index to Ind3D
    k_corner[i.ind] = Ind3D((x*ncy + y)*ncz + corner_zind, ncy, ncz);

    // Check that t_z is in range
    if ((t_z < 0.0) || (t_z > 1.0)) {
      throw BoutException(
          "t_z={:e} out of range at ({:d},{:d},{:d}) (delta_z={:e}, k_corner={:d})", t_z,
          x, y, z, delta_z(x, y, z), k_corner[i.ind].ind);
    }

    h00(x, y, z) = (2. * t_z * t_z * t_z) - (3. * t_z * t_z) + 1.;
    h01(x, y, z) = (-2. * t_z * t_z * t_z) + (3. * t_z * t_z);
    h10(x, y, z) = t_z * (1. - t_z) * (1. - t_z);
    h11(x, y, z) = (t_z * t_z * t_z) - (t_z * t_z);
  }

#if CHECK > 2
  bout::checkFinite(h00, "h00", "RGN_NOY");
  bout::checkFinite(h01, "h01", "RGN_NOY");
  bout::checkFinite(h10, "h10", "RGN_NOY");
  bout::checkFinite(h11, "h11", "RGN_NOY");
#endif
}

/*!
 * Return position and weight of points needed to approximate the function value at the
 * point that the field line through (i,j,k) meets the (j+1)-plane. For the case where
 * only the z-direction is not aligned to grid points, the approximation is: f(i,j+1,k*) =
 * h00 * f(i,j+1,k) + h01 * f(i,j+1,k+1)
 *                 + h10 * dfdz(i,j+1,k) + h11 * dfdz(i,j+1,k+1)
 *               = h00 * f(i,j+1,k) + h01 * f(i,j+1,k+1)
 *                 + h10 * ( f(i,j+1,k+1) - f(i,j+1,k-1) ) / 2
 *                 + h11 * ( f(i,j+1,k+2) - f(i,j+1,k) ) / 2
 * for k* a point between k and k+1. Therefore, this function returns
 *   position 		weight
 *   (i, j+1, k-1)	- h10 / 2
 *   (i, j+1, k)	h00 - h11 / 2
 *   (i, j+1, k+1)	h01 + h10 / 2
 *   (i, j+1, k+2)	h11 / 2
 */
std::vector<ParallelTransform::PositionsAndWeights>
ZHermiteSpline::getWeightsForYApproximation(int i, int j, int k, int yoffset) const {
  ASSERT3(i >= 0);
  ASSERT3(i <= localmesh->LocalNx);
  ASSERT3(j >= localmesh->ystart);
  ASSERT3(j <= localmesh->yend);
  ASSERT3(k >= 0);
  ASSERT3(k <= localmesh->LocalNz);

  const int ncz = localmesh->LocalNz;
  const auto corner = k_corner[(i*localmesh->LocalNy + j)*ncz + k];
  const int k_mod = corner.z();
  const int k_mod_m1 = corner.zm().z();
  const int k_mod_p1 = corner.zp().z();
  const int k_mod_p2 = corner.zpp().z();

  return {{i, j + yoffset, k_mod_m1, -0.5 * h10(i, j, k)},
          {i, j + yoffset, k_mod,    h00(i, j, k) - 0.5 * h11(i, j, k)},
          {i, j + yoffset, k_mod_p1, h01(i, j, k) + 0.5 * h10(i, j, k)},
          {i, j + yoffset, k_mod_p2, 0.5 * h11(i, j, k)}};
}

Field3D ZHermiteSpline::interpolate(const Field3D& f, const std::string& region_str) const {

  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  std::string local_fz_region;
  if (region_str == "DEFAULT") {
    local_fz_region = fz_region;
  } else {
    // Cannot use non-default regions when y_offset!=0 as calcWeights() used the default
    // region
    ASSERT1(y_offset == 0);
    local_fz_region = region_str;
  }
  const auto& local_region = (region_str == "DEFAULT") ? region : f.getRegion(region_str);

  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  Field3D fz = bout::derivatives::index::DDZ(f, CELL_DEFAULT, "DEFAULT", local_fz_region);

  BOUT_FOR(i, local_region) {
    const auto corner = k_corner[i.ind].yp(y_offset);
    const auto corner_zp1 = corner.zp();

    // Interpolate in Z
    f_interp[i.yp(y_offset)] = f[corner] * h00[i] + f[corner_zp1] * h01[i]
                               + fz[corner] * h10[i] + fz[corner_zp1] * h11[i];

    ASSERT2(finite(f_interp[i.yp(y_offset)]) || i.x() < localmesh->xstart
            || i.x() > localmesh->xend);
  }
  return f_interp;
}

Field3D ZHermiteSpline::interpolate(const Field3D& f, const Field3D& delta_z,
                                    const std::string& region_str) {
  calcWeights(delta_z);
  return interpolate(f, region_str);
}

void ZInterpolationFactory::ensureRegistered() {}

namespace {
RegisterZInterpolation<ZHermiteSpline> registerzinterphermitespline{"hermitespline"};
} // namespace

constexpr decltype(ZInterpolationFactory::type_name) ZInterpolationFactory::type_name;
constexpr decltype(
    ZInterpolationFactory::section_name) ZInterpolationFactory::section_name;
constexpr decltype(ZInterpolationFactory::option_name) ZInterpolationFactory::option_name;
constexpr decltype(
    ZInterpolationFactory::default_type) ZInterpolationFactory::default_type;
