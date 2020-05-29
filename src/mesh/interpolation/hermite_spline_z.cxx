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

ZHermiteSpline::ZHermiteSpline(BoutMask mask, int y_offset, Mesh* mesh)
    : ZInterpolation(mask, y_offset, mesh), h00(localmesh), h01(localmesh),
      h10(localmesh), h11(localmesh) {

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

void ZHermiteSpline::calcWeights(const Field3D& delta_z, const std::string& region) {

  const int ncy = localmesh->LocalNy;
  const int ncz = localmesh->LocalNz;

  BOUT_FOR(i, delta_z.getRegion(region)) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    if (skip_mask(x, y, z))
      continue;

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
}

void ZHermiteSpline::calcWeights(const Field3D& delta_z, const BoutMask& mask,
                                 const std::string& region) {
  skip_mask = mask;
  has_mask = true;
  calcWeights(delta_z, region);
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

Field3D ZHermiteSpline::interpolate(const Field3D& f, const std::string& region) const {
  // Template with two branches for the body of this method so that if
  // has_mask=false we can optimize out the conditional with skip_mask from the
  // tight loop.
  if (has_mask) {
    return interpolate_internal<true>(f, region);
  }
  return interpolate_internal<false>(f, region);
}

template<bool with_mask>
Field3D ZHermiteSpline::interpolate_internal(const Field3D& f, const std::string& region) const {

  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  const std::string fz_region = (region == "RGN_NOBNDRY" or region == "RGN_NOX")
                                    ? y_offset == 0 ? "RGN_NOBNDRY" : "RGN_NOX"
                                    : "RGN_ALL";
  Field3D fz = bout::derivatives::index::DDZ(f, CELL_DEFAULT, "DEFAULT", fz_region);
  if (region != "RGN_NOBNDRY" and region != "RGN_NOX") {
    localmesh->communicateXZ(fz);
  }
  if (y_offset != 0 or (region != "RGN_NOBNDRY" and region != "RGN_NOY")) {
    // communicate in y, but do not calculate parallel slices
    auto h = localmesh->sendY(fz);
    localmesh->wait(h);
  }

  BOUT_FOR(i, f.getRegion(region)) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    if (with_mask) {
      if (skip_mask(x, y, z))
        continue;
    }

    // Due to lack of guard cells in z-direction, we need to ensure z-index
    // wraps around
    const auto corner = k_corner[i.ind].yp(y_offset);
    const auto corner_zp1 = corner.zp();

    // Interpolate in Z
    f_interp[i.yp(y_offset)] = f[corner] * h00[i] + f[corner_zp1] * h01[i]
                               + fz[corner] * h10[i] + fz[corner_zp1] * h11[i];

    ASSERT2(finite(f_interp[i.yp(y_offset)]) || x < localmesh->xstart
            || x > localmesh->xend);
  }
  return f_interp;
}

Field3D ZHermiteSpline::interpolate(const Field3D& f, const Field3D& delta_z,
                                    const std::string& region) {
  calcWeights(delta_z, region);
  return interpolate(f, region);
}

Field3D ZHermiteSpline::interpolate(const Field3D& f, const Field3D& delta_z,
                                    const BoutMask& mask, const std::string& region) {
  calcWeights(delta_z, mask, region);
  return interpolate(f, region);
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
