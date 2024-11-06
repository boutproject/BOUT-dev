/**************************************************************************
 * Copyright 2024 The BOUT++ Team
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

#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/globalindexer.hxx"
#include "bout/globals.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/operatorstencil.hxx"
#include "bout/paralleltransform.hxx"
#include "bout/petsc_interface.hxx"
#include "bout/petsclib.hxx"
#include "bout/region.hxx"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

PetscXZHermiteSpline::PetscXZHermiteSpline(int y_offset, Mesh* mesh)
    : XZInterpolation(y_offset, mesh),
      petsclib(&Options::root()["mesh:paralleltransform:xzinterpolation:hermitespline"]),
      indexer(std::make_shared<bout::SimpleGlobalIndexer>(localmesh)),
      weights(indexer, false), h00_x(localmesh), h01_x(localmesh), h10_x(localmesh),
      h11_x(localmesh), h00_z(localmesh), h01_z(localmesh), h10_z(localmesh),
      h11_z(localmesh) {

  // Index arrays contain guard cells in order to get subscripts right
  i_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);
  k_corner.reallocate(localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz);

  // Initialise in order to avoid 'uninitialized value' errors from Valgrind when using
  // guard-cell values
  k_corner = -1;

  // Allocate Field3D members
  h00_x.allocate();
  h01_x.allocate();
  h10_x.allocate();
  h11_x.allocate();
  h00_z.allocate();
  h01_z.allocate();
  h10_z.allocate();
  h11_z.allocate();

  // The stencil has 16 elements, so maximum number of columns in any
  // given row (whether on local or remote process) is 16
  // NOLINTNEXTLINE(misc-include-cleaner)
  MatMPIAIJSetPreallocation(*weights.get(), 16, nullptr, 16, nullptr);
}

void PetscXZHermiteSpline::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                                       const std::string& region) {

  const int ny = localmesh->LocalNy;
  const int nz = localmesh->LocalNz;
  const int xend = (localmesh->xend - localmesh->xstart + 1) * localmesh->getNXPE()
                   + localmesh->xstart - 1;

  // TODO: work out why using `getRegion(region)` directly causes
  // issues on more than one core
  const auto actual_region = getRegion(region);
  BOUT_FOR(i, actual_region) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    // The integer part of xt_prime, zt_prime are the indices of the cell
    // containing the field line end-point
    int i_corn = static_cast<int>(std::floor(delta_x(x, y, z)));
    k_corner(x, y, z) = static_cast<int>(std::floor(delta_z(x, y, z)));

    // t_x, t_z are the normalised coordinates \in [0,1) within the cell
    // calculated by taking the remainder of the floating point index
    BoutReal t_x = delta_x(x, y, z) - static_cast<BoutReal>(i_corn);
    const BoutReal t_z = delta_z(x, y, z) - static_cast<BoutReal>(k_corner(x, y, z));

    // NOTE: A (small) hack to avoid one-sided differences
    if (i_corn >= xend) {
      i_corn = xend - 1;
      t_x = 1.0;
    }
    if (i_corn < localmesh->xstart) {
      i_corn = localmesh->xstart;
      t_x = 0.0;
    }

    k_corner(x, y, z) = ((k_corner(x, y, z) % nz) + nz) % nz;

    // Check that t_x and t_z are in range
    if ((t_x < 0.0) || (t_x > 1.0)) {
      throw BoutException(
          "t_x={:e} out of range at ({:d},{:d},{:d}) (delta_x={:e}, i_corn={:d})", t_x, x,
          y, z, delta_x(x, y, z), i_corn);
    }

    if ((t_z < 0.0) || (t_z > 1.0)) {
      throw BoutException(
          "t_z={:e} out of range at ({:d},{:d},{:d}) (delta_z={:e}, k_corner={:d})", t_z,
          x, y, z, delta_z(x, y, z), k_corner(x, y, z));
    }

    h00_x[i] = (2. * t_x * t_x * t_x) - (3. * t_x * t_x) + 1.;
    h00_z[i] = (2. * t_z * t_z * t_z) - (3. * t_z * t_z) + 1.;

    h01_x[i] = (-2. * t_x * t_x * t_x) + (3. * t_x * t_x);
    h01_z[i] = (-2. * t_z * t_z * t_z) + (3. * t_z * t_z);

    h10_x[i] = t_x * (1. - t_x) * (1. - t_x);
    h10_z[i] = t_z * (1. - t_z) * (1. - t_z);

    h11_x[i] = (t_x * t_x * t_x) - (t_x * t_x);
    h11_z[i] = (t_z * t_z * t_z) - (t_z * t_z);

    // Need to convert from global indices to local
    const auto i_c = Field3D::ind_type(
        (((localmesh->getLocalXIndex(i_corn) * ny) + (y + y_offset)) * nz
         + localmesh->getLocalZIndex(k_corner(x, y, z))),
        ny, nz);

    // f[ic] * h00_x[i] + f[icxp] * h01_x[i] + fx[ic] * h10_x[i] + fx[icxp] * h11_x[i];
    // f[iczp] * h00_x[i] + f[icxpzp] * h01_x[i] + fx[iczp] * h10_x[i] + fx[icxpzp] * h11_x[i];
    // fz[ic] * h00_x[i] + fz[icxp] * h01_x[i] + fxz[ic] * h10_x[i]+ fxz[icxp] * h11_x[i];
    // fz[iczp] * h00_x[i] + fz[icxpzp] * h01_x[i] + fxz[iczp] * h10_x[i] + fxz[icxpzp] * h11_x[i];

    // Quite a few duplicated terms, could pre-calculate
    weights(i, i_c.xm()) = (h10_x[i] * h11_z[i] / 4) - (h10_x[i] * h00_z[i] / 2);
    weights(i, i_c.xm().zm()) = (h10_x[i] * h10_z[i] / 4);
    weights(i, i_c.xm().zp()) = -(h10_x[i] * h01_z[i] / 2) - (h10_x[i] * h10_z[i] / 4);
    weights(i, i_c.xm().zp(2)) = -(h10_x[i] * h11_z[i] / 4);
    weights(i, i_c.xp()) = (h01_x[i] * h00_z[i]) + (h10_x[i] * h00_z[i] / 2)
                           - (h01_x[i] * h11_z[i] / 2) - (h10_x[i] * h11_z[i] / 4);
    weights(i, i_c.xp().zm()) = -(h01_x[i] * h10_z[i] / 2) - (h10_x[i] * h10_z[i] / 4);
    weights(i, i_c.xp().zp()) = (h01_x[i] * h01_z[i]) + (h01_x[i] * h10_z[i] / 2)
                                + (h10_x[i] * h01_z[i] / 2) + (h10_x[i] * h10_z[i] / 4);
    weights(i, i_c.xp().zp(2)) = (h01_x[i] * h11_z[i] / 2) + (h10_x[i] * h11_z[i] / 4);
    weights(i, i_c.xp(2)) = (h11_x[i] * h00_z[i] / 2) - (h11_x[i] * h11_z[i] / 4);
    weights(i, i_c.xp(2).zm()) = -(h11_x[i] * h10_z[i] / 4);
    weights(i, i_c.xp(2).zp()) = (h11_x[i] * h01_z[i] / 2) + (h11_x[i] * h10_z[i] / 4);
    weights(i, i_c.xp(2).zp(2)) = (h11_x[i] * h11_z[i] / 4);
    weights(i, i_c) = (h00_x[i] * h00_z[i]) + (h11_x[i] * h11_z[i] / 4)
                      - (h00_x[i] * h11_z[i] / 2) - (h11_x[i] * h00_z[i] / 2);
    weights(i, i_c.zm()) = (h11_x[i] * h10_z[i] / 4) - (h00_x[i] * h10_z[i] / 2);
    weights(i, i_c.zp()) = (h00_x[i] * h01_z[i]) + (h00_x[i] * h10_z[i] / 2)
                           - (h11_x[i] * h01_z[i] / 2) - (h11_x[i] * h10_z[i] / 4);
    weights(i, i_c.zp(2)) = (h00_x[i] * h11_z[i] / 2) - (h11_x[i] * h11_z[i] / 4);
  }

  weights.assemble();
}

void PetscXZHermiteSpline::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                                       const BoutMask& mask, const std::string& region) {
  setMask(mask);
  calcWeights(delta_x, delta_z, region);
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
PetscXZHermiteSpline::getWeightsForYApproximation(int i, int j, int k, int yoffset) {
  const int nz = localmesh->LocalNz;
  const int k_mod = k_corner(i, j, k);
  const int k_mod_m1 = (k_mod > 0) ? (k_mod - 1) : (nz - 1);
  const int k_mod_p1 = (k_mod == nz) ? 0 : k_mod + 1;
  const int k_mod_p2 = (k_mod_p1 == nz) ? 0 : k_mod_p1 + 1;

  return {{i, j + yoffset, k_mod_m1, -0.5 * h10_z(i, j, k)},
          {i, j + yoffset, k_mod, h00_z(i, j, k) - 0.5 * h11_z(i, j, k)},
          {i, j + yoffset, k_mod_p1, h01_z(i, j, k) + 0.5 * h10_z(i, j, k)},
          {i, j + yoffset, k_mod_p2, 0.5 * h11_z(i, j, k)}};
}

Field3D
PetscXZHermiteSpline::interpolate(const Field3D& f,
                                  [[maybe_unused]] const std::string& region) const {

  ASSERT1(f.getMesh() == localmesh);
  const PetscVector<Field3D> f_petsc(f, indexer);
  const PetscVector<Field3D> f_result = weights * f_petsc;
  // TODO: we should only consider the passed-in region
  // const auto region2 = y_offset == 0 ? region : fmt::format("RGN_YPAR_{:+d}", y_offset);
  return f_result.toField();
}

Field3D PetscXZHermiteSpline::interpolate(const Field3D& f, const Field3D& delta_x,
                                          const Field3D& delta_z,
                                          const std::string& region) {
  calcWeights(delta_x, delta_z, region);
  return interpolate(f, region);
}

Field3D PetscXZHermiteSpline::interpolate(const Field3D& f, const Field3D& delta_x,
                                          const Field3D& delta_z, const BoutMask& mask,
                                          const std::string& region) {
  calcWeights(delta_x, delta_z, mask, region);
  return interpolate(f, region);
}

#endif
