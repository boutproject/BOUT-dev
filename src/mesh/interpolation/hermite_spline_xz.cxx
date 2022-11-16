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
#include "interpolation_xz.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/mesh.hxx"

#include <vector>

XZHermiteSpline::XZHermiteSpline(int y_offset, Mesh *mesh)
    : XZInterpolation(y_offset, mesh),
      h00_x(localmesh), h01_x(localmesh), h10_x(localmesh), h11_x(localmesh),
      h00_z(localmesh), h01_z(localmesh), h10_z(localmesh), h11_z(localmesh) {

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


  newWeights.reserve(16);
  for (int w=0; w<16;++w){
    newWeights.emplace_back(localmesh);
    newWeights[w].allocate();
  }
}

void XZHermiteSpline::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                                  const std::string& region) {

  const int ny = localmesh->LocalNy;
  const int nz = localmesh->LocalNz;
  BOUT_FOR(i, getRegion(region)) {
    const int x = i.x();
    const int y = i.y();
    const int z = i.z();

    // The integer part of xt_prime, zt_prime are the indices of the cell
    // containing the field line end-point
    int i_corn = static_cast<int>(floor(delta_x(x, y, z)));
    k_corner(x, y, z) = static_cast<int>(floor(delta_z(x, y, z)));

    // t_x, t_z are the normalised coordinates \in [0,1) within the cell
    // calculated by taking the remainder of the floating point index
    BoutReal t_x = delta_x(x, y, z) - static_cast<BoutReal>(i_corn);
    BoutReal t_z = delta_z(x, y, z) - static_cast<BoutReal>(k_corner(x, y, z));

    // NOTE: A (small) hack to avoid one-sided differences
    if (i_corn >= localmesh->xend) {
      i_corn = localmesh->xend - 1;
      t_x = 1.0;
    }
    if (i_corn < localmesh->xstart) {
      i_corn = localmesh->xstart;
      t_x = 0.0;
    }

    k_corner(x, y, z) = ((k_corner(x,y,z) % nz) + nz) % nz;

    // Check that t_x and t_z are in range
    if ((t_x < 0.0) || (t_x > 1.0)) {
      throw BoutException(
          "t_x={:e} out of range at ({:d},{:d},{:d}) (delta_x={:e}, i_corn={:d})", t_x,
          x, y, z, delta_x(x, y, z), i_corn);
    }

    if ((t_z < 0.0) || (t_z > 1.0)) {
      throw BoutException(
          "t_z={:e} out of range at ({:d},{:d},{:d}) (delta_z={:e}, k_corner={:d})", t_z,
          x, y, z, delta_z(x, y, z), k_corner(x, y, z));
    }

    i_corner[i] = SpecificInd<IND_TYPE::IND_3D>(
        (((i_corn * ny) + (y + y_offset)) * nz + k_corner(x, y, z)), ny, nz);

    h00_x[i] = (2. * t_x * t_x * t_x) - (3. * t_x * t_x) + 1.;
    h00_z[i] = (2. * t_z * t_z * t_z) - (3. * t_z * t_z) + 1.;

    h01_x[i] = (-2. * t_x * t_x * t_x) + (3. * t_x * t_x);
    h01_z[i] = (-2. * t_z * t_z * t_z) + (3. * t_z * t_z);

    h10_x[i] = t_x * (1. - t_x) * (1. - t_x);
    h10_z[i] = t_z * (1. - t_z) * (1. - t_z);

    h11_x[i] = (t_x * t_x * t_x) - (t_x * t_x);
    h11_z[i] = (t_z * t_z * t_z) - (t_z * t_z);

#define USE_NEW_WEIGHTS 1
#if USE_NEW_WEIGHTS

    for (int w =0; w<16;++w){
      newWeights[w][i]=0;
    }
    // The distribution of our weights:
    //  0   4   8    12
    //  1   5   9    13
    //  2   6   10   14
    //  3   7   11   15
    // e.g. 1 == ic.xm(); 4 == ic.zm(); 5 == ic;  7 == ic.zp(2);

    // f[ic] * h00_x[i] + f[icxp] * h01_x[i] + fx[ic] * h10_x[i] + fx[icxp] * h11_x[i];
    newWeights[5][i]  += h00_x[i] * h00_z[i];
    newWeights[9][i]  += h01_x[i] * h00_z[i];
    newWeights[9][i]  += h10_x[i] * h00_z[i] / 2;
    newWeights[1][i]  -= h10_x[i] * h00_z[i] / 2;
    newWeights[13][i] += h11_x[i] * h00_z[i] / 2;
    newWeights[5][i]  -= h11_x[i] * h00_z[i] / 2;

    // f[iczp] * h00_x[i] + f[icxpzp] * h01_x[i] +
    // fx[iczp] * h10_x[i] + fx[icxpzp] * h11_x[i];
    newWeights[6][i]  += h00_x[i] * h01_z[i];
    newWeights[10][i] += h01_x[i] * h01_z[i];
    newWeights[10][i] += h10_x[i] * h01_z[i] / 2;
    newWeights[2][i]  -= h10_x[i] * h01_z[i] / 2;
    newWeights[14][i] += h11_x[i] * h01_z[i] / 2;
    newWeights[6][i]  -= h11_x[i] * h01_z[i] / 2;

    // fz[ic] * h00_x[i] + fz[icxp] * h01_x[i] +
    // fxz[ic] * h10_x[i]+ fxz[icxp] * h11_x[i];
    newWeights[6][i]  += h00_x[i] * h10_z[i] / 2;
    newWeights[4][i]  -= h00_x[i] * h10_z[i] / 2;
    newWeights[10][i] += h01_x[i] * h10_z[i] / 2;
    newWeights[8][i]  -= h01_x[i] * h10_z[i] / 2;
    newWeights[10][i] += h10_x[i] * h10_z[i] / 4;
    newWeights[8][i]  -= h10_x[i] * h10_z[i] / 4;
    newWeights[2][i]  -= h10_x[i] * h10_z[i] / 4;
    newWeights[0][i]  += h10_x[i] * h10_z[i] / 4;
    newWeights[14][i] += h11_x[i] * h10_z[i] / 4;
    newWeights[12][i] -= h11_x[i] * h10_z[i] / 4;
    newWeights[6][i] -= h11_x[i] * h10_z[i] / 4;
    newWeights[4][i]  += h11_x[i] * h10_z[i] / 4;

    // fz[iczp] * h00_x[i] + fz[icxpzp] * h01_x[i] +
    // fxz[iczp] * h10_x[i] + fxz[icxpzp] * h11_x[i];
    newWeights[7][i]  += h00_x[i] * h11_z[i] / 2;
    newWeights[5][i]  -= h00_x[i] * h11_z[i] / 2;
    newWeights[11][i] += h01_x[i] * h11_z[i] / 2;
    newWeights[9][i]  -= h01_x[i] * h11_z[i] / 2;
    newWeights[11][i] += h10_x[i] * h11_z[i] / 4;
    newWeights[9][i]  -= h10_x[i] * h11_z[i] / 4;
    newWeights[3][i]  -= h10_x[i] * h11_z[i] / 4;
    newWeights[1][i]  += h10_x[i] * h11_z[i] / 4;
    newWeights[15][i] += h11_x[i] * h11_z[i] / 4;
    newWeights[13][i] -= h11_x[i] * h11_z[i] / 4;
    newWeights[7][i]  -= h11_x[i] * h11_z[i] / 4;
    newWeights[5][i]  += h11_x[i] * h11_z[i] / 4;


    //   // f[ic] * h00_x[i] + f[icxp] * h01_x[i] + fx[ic] * h10_x[i] + fx[icxp] * h11_x[i];
    // newWeights[5][i]  += h00_x[i] * h00_z[i];
    // newWeights[9][i]  += h01_x[i] * h00_z[i];
    // newWeights[9][i]  += h10_x[i] * h00_z[i] / 2 / localmesh->dx[ic];
    // newWeights[1][i]  -= h10_x[i] * h00_z[i] / 2 / localmesh->dx[ic];
    // newWeights[13][i] += h11_x[i] * h00_z[i] / 2 / localmesh->dx[ic.xp()];
    // newWeights[5][i]  -= h11_x[i] * h00_z[i] / 2 / localmesh->dx[ic.xp()];

    // // f[iczp] * h00_x[i] + f[icxpzp] * h01_x[i] +
    // // fx[iczp] * h10_x[i] + fx[icxpzp] * h11_x[i];
    // newWeights[6][i]  += h00_x[i] * h01_z[i];
    // newWeights[10][i]  += h01_x[i] * h01_z[i];
    // newWeights[10][i]  += h10_x[i] * h01_z[i] / 2/ localmesh->dx[ic.zp()];
    // newWeights[2][i]  -= h10_x[i] * h01_z[i] / 2/ localmesh->dx[ic.zp()];
    // newWeights[14][i] += h11_x[i] * h01_z[i] / 2/ localmesh->dx[ic.zp().xp()];
    // newWeights[6][i]  -= h11_x[i] * h01_z[i] / 2/ localmesh->dx[ic.zp().xp()];

    // // fz[ic] * h00_x[i] + fz[icxp] * h01_x[i] +
    // // fxz[ic] * h10_x[i]+ fxz[icxp] * h11_x[i];
    // newWeights[6][i]  += h00_x[i] * h10_z[i] / 2 / localmesh->dz[ic];
    // newWeights[4][i]  -= h00_x[i] * h10_z[i] / 2 / localmesh->dz[ic];
    // newWeights[10][i] += h01_x[i] * h10_z[i] / 2 / localmesh->dz[ic.xp()];
    // newWeights[8][i]  -= h01_x[i] * h10_z[i] / 2 / localmesh->dz[ic.xp()];
    // newWeights[10][i] += h10_x[i] * h10_z[i] / 4 / localmesh->dz[ic] / localmesh->dx[ic];
    // newWeights[8][i]  -= h10_x[i] * h10_z[i] / 4 / localmesh->dz[ic] / localmesh->dx[ic];
    // newWeights[2][i]  -= h10_x[i] * h10_z[i] / 4 / localmesh->dz[ic] / localmesh->dx[ic];
    // newWeights[0][i]  += h10_x[i] * h10_z[i] / 4 / localmesh->dz[ic] / localmesh->dx[ic];
    // newWeights[14][i] += h11_x[i] * h10_z[i] / 4 / localmesh->dz[ic.xp()] / localmesh->dx[ic.xp()];
    // newWeights[12][i] -= h11_x[i] * h10_z[i] / 4 / localmesh->dz[ic.xp()] / localmesh->dx[ic.xp()];
    // newWeights[6][i] -= h11_x[i] * h10_z[i] / 4 / localmesh->dz[ic.xp()] / localmesh->dx[ic.xp()];
    // newWeights[4][i]  += h11_x[i] * h10_z[i] / 4 / localmesh->dz[ic.xp()] / localmesh->dx[ic.xp()];

    // // fz[iczp] * h00_x[i] + fz[icxpzp] * h01_x[i] +
    // // fxz[iczp] * h10_x[i] + fxz[icxpzp] * h11_x[i];
    // newWeights[7][i]  += h00_x[i] * h11_z[i] / 2 / localmesh->dz[ic.zp()];
    // newWeights[5][i]  -= h00_x[i] * h11_z[i] / 2 / localmesh->dz[ic.zp()];
    // newWeights[11][i] += h01_x[i] * h11_z[i] / 2 / localmesh->dz[ic.zp().xp()];
    // newWeights[9][i]  -= h01_x[i] * h11_z[i] / 2 / localmesh->dz[ic.zp().xp()];
    // newWeights[11][i] += h10_x[i] * h11_z[i] / 4 / localmesh->dz[ic.zp()] / localmesh->dx[ic.zp()];
    // newWeights[9][i]  -= h10_x[i] * h11_z[i] / 4 / localmesh->dz[ic.zp()] / localmesh->dx[ic.zp()];
    // newWeights[3][i]  -= h10_x[i] * h11_z[i] / 4 / localmesh->dz[ic.zp()] / localmesh->dx[ic.zp()];
    // newWeights[1][i]  += h10_x[i] * h11_z[i] / 4 / localmesh->dz[ic.zp()] / localmesh->dx[ic.zp()];
    // newWeights[15][i] += h11_x[i] * h11_z[i] / 4 / localmesh->dz[ic.zp().xp()] / localmesh->dx[ic.zp().xp()];
    // newWeights[13][i] -= h11_x[i] * h11_z[i] / 4 / localmesh->dz[ic.zp().xp()] / localmesh->dx[ic.zp().xp()];
    // newWeights[7][i]  -= h11_x[i] * h11_z[i] / 4 / localmesh->dz[ic.zp().xp()] / localmesh->dx[ic.zp().xp()];
    // newWeights[5][i]  += h11_x[i] * h11_z[i] / 4 / localmesh->dz[ic.zp().xp()] / localmesh->dx[ic.zp().xp()];
#endif
  }
}

void XZHermiteSpline::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
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
XZHermiteSpline::getWeightsForYApproximation(int i, int j, int k, int yoffset) {
  const int ncz = localmesh->LocalNz;
  const int k_mod = k_corner(i, j, k);
  const int k_mod_m1 = (k_mod > 0) ? (k_mod - 1) : (ncz - 1);
  const int k_mod_p1 = (k_mod == ncz) ? 0 : k_mod + 1;
  const int k_mod_p2 = (k_mod_p1 == ncz) ? 0 : k_mod_p1 + 1;

  return {{i, j + yoffset, k_mod_m1, -0.5 * h10_z(i, j, k)},
          {i, j + yoffset, k_mod,    h00_z(i, j, k) - 0.5 * h11_z(i, j, k)},
          {i, j + yoffset, k_mod_p1, h01_z(i, j, k) + 0.5 * h10_z(i, j, k)},
          {i, j + yoffset, k_mod_p2, 0.5 * h11_z(i, j, k)}};
}

Field3D XZHermiteSpline::interpolate(const Field3D& f, const std::string& region) const {

  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};


#if USE_NEW_WEIGHTS
  BOUT_FOR(i, getRegion(region)) {
    auto ic =  i_corner[i];
    auto iyp = i.yp(y_offset);

    f_interp[iyp]=0;
    for (int w = 0; w < 4; ++w){
      f_interp[iyp] += newWeights[w*4+0][i] * f[ic.zm().xp(w-1)];
      f_interp[iyp] += newWeights[w*4+1][i] * f[ic.xp(w-1)];
      f_interp[iyp] += newWeights[w*4+2][i] * f[ic.zp().xp(w-1)];
      f_interp[iyp] += newWeights[w*4+3][i] * f[ic.zp(2).xp(w-1)];
    }
  }
  return f_interp;
#else
  // Derivatives are used for tension and need to be on dimensionless
  // coordinates
  const auto region2 = fmt::format("RGN_YPAR_{:+d}", y_offset);
  // f has been communcated, and thus we can assume that the x-boundaries are
  // also valid in the y-boundary.  Thus the differentiated field needs no
  // extra comms.
  Field3D fx = bout::derivatives::index::DDX(f, CELL_DEFAULT, "DEFAULT", region2);
  Field3D fz = bout::derivatives::index::DDZ(f, CELL_DEFAULT, "DEFAULT", region2);
  Field3D fxz = bout::derivatives::index::DDZ(fx, CELL_DEFAULT, "DEFAULT", region2);

  BOUT_FOR(i, getRegion(region)) {
    const auto iyp = i.yp(y_offset);

    const auto ic = i_corner[i];
    const auto iczp = ic.zp();
    const auto icxp = ic.xp();
    const auto icxpzp = iczp.xp();

    // Interpolate f in X at Z
    const BoutReal f_z =
        f[ic] * h00_x[i] + f[icxp] * h01_x[i] + fx[ic] * h10_x[i] + fx[icxp] * h11_x[i];

    // Interpolate f in X at Z+1
    const BoutReal f_zp1 = f[iczp] * h00_x[i] + f[icxpzp] * h01_x[i] + fx[iczp] * h10_x[i]
                           + fx[icxpzp] * h11_x[i];

    // Interpolate fz in X at Z
    const BoutReal fz_z = fz[ic] * h00_x[i] + fz[icxp] * h01_x[i] + fxz[ic] * h10_x[i]
                          + fxz[icxp] * h11_x[i];

    // Interpolate fz in X at Z+1
    const BoutReal fz_zp1 = fz[iczp] * h00_x[i] + fz[icxpzp] * h01_x[i]
                            + fxz[iczp] * h10_x[i] + fxz[icxpzp] * h11_x[i];

    // Interpolate in Z
    f_interp[iyp] =
        +f_z * h00_z[i] + f_zp1 * h01_z[i] + fz_z * h10_z[i] + fz_zp1 * h11_z[i];

    ASSERT2(std::isfinite(f_interp[iyp]) || i.x() < localmesh->xstart
            || i.x() > localmesh->xend);
  }
  return f_interp;
# endif
}

Field3D XZHermiteSpline::interpolate(const Field3D& f, const Field3D& delta_x,
                                     const Field3D& delta_z, const std::string& region) {
  calcWeights(delta_x, delta_z, region);
  return interpolate(f, region);
}

Field3D XZHermiteSpline::interpolate(const Field3D& f, const Field3D& delta_x,
                                     const Field3D& delta_z, const BoutMask& mask,
                                     const std::string& region) {
  calcWeights(delta_x, delta_z, mask, region);
  return interpolate(f, region);
}
