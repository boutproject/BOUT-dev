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

#include <bout/build_defines.hxx>

#if BOUT_HAS_PETSC

#include "petsc_hermite_spline_xz.hxx"

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/globals.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/mesh.hxx"
#include "bout/openmpwrap.hxx"
#include "bout/paralleltransform.hxx"
#include "bout/region.hxx"

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace {
class IndConverter {
public:
  IndConverter(Mesh* mesh) : mesh(mesh), nxpe(mesh->getNXPE()), nype(mesh->getNYPE()) {}
  // ix and iy are global indices
  // iy is local
  int fromMeshToGlobal(int ix, int iy, int iz) {
    const int xstart = mesh->xstart;
    const int lnx = mesh->LocalNx - (xstart * 2);
    // x-proc-id
    const int pex = std::clamp(divToNeg(ix - xstart, lnx), 0, nxpe - 1);
    const int zstart = 0;
    const int lnz = mesh->LocalNz - (zstart * 2);
    // z-proc-id
    // pez only for wrapping around ; later needs similar treatment than pey
    const int pez = divToNeg(iz - zstart, lnz);
    // y proc-id - y is already local
    const int ystart = mesh->ystart;
    const int lny = mesh->LocalNy - (ystart * 2);
    const int pey_offset = divToNeg(iy - ystart, lny);
    int pey = pey_offset + mesh->getYProcIndex();
    while (pey < 0) {
      pey += nype;
    }
    while (pey >= nype) {
      pey -= nype;
    }
    ASSERT2(pex >= 0);
    ASSERT2(pex < nxpe);
    ASSERT2(pey >= 0);
    ASSERT2(pey < nype);
    return fromLocalToGlobal(ix - (pex * lnx), iy - (pey_offset * lny), iz - (pez * lnz),
                             pex, pey, 0);
  }
  int fromLocalToGlobal(const int ilocalx, const int ilocaly, const int ilocalz) {
    return fromLocalToGlobal(ilocalx, ilocaly, ilocalz, mesh->getXProcIndex(),
                             mesh->getYProcIndex(), 0);
  }
  int fromLocalToGlobal(const int ilocalx, const int ilocaly, const int ilocalz,
                        const int pex, const int pey, const int pez) {
    ASSERT3(ilocalx >= 0);
    ASSERT3(ilocaly >= 0);
    ASSERT3(ilocalz >= 0);
    const int ilocal = (((ilocalx * mesh->LocalNy) + ilocaly) * mesh->LocalNz) + ilocalz;
    const int ret = ilocal
                    + (mesh->LocalNx * mesh->LocalNy * mesh->LocalNz
                       * ((pey * nxpe + pex) * nzpe + pez));
    ASSERT3(ret >= 0);
    ASSERT3(ret < nxpe * nype * mesh->LocalNx * mesh->LocalNy * mesh->LocalNz);
    return ret;
  }

private:
  // number of procs
  Mesh* mesh;
  int nxpe;
  int nype;
  int nzpe{1};
  static int divToNeg(const int n, const int d) {
    return (n < 0) ? ((n - d + 1) / d) : (n / d);
  }
};
} // namespace

XZPetscHermiteSpline::XZPetscHermiteSpline(int y_offset, Mesh* meshin)
    : XZInterpolation(y_offset, meshin), h00_x(localmesh), h01_x(localmesh),
      h10_x(localmesh), h11_x(localmesh), h00_z(localmesh), h01_z(localmesh),
      h10_z(localmesh), h11_z(localmesh),
      petsclib(
          &Options::root()["mesh:paralleltransform:xzinterpolation:petschermitespline"]) {

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

  const int m = localmesh->LocalNx * localmesh->LocalNy * localmesh->LocalNz;
  const int M = m * localmesh->getNXPE() * localmesh->getNYPE();
  MatCreateAIJ(BoutComm::get(), m, m, M, M, 16, nullptr, 16, nullptr, &petscWeights);
}

void XZPetscHermiteSpline::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                                       const std::string& region) {

  const int ny = localmesh->LocalNy;
  const int nz = localmesh->LocalNz;
  const int xend = ((localmesh->xend - localmesh->xstart + 1) * localmesh->getNXPE())
                   + localmesh->xstart - 1;

  IndConverter conv{localmesh};

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
    const BoutReal t_z = delta_z(x, y, z) - static_cast<BoutReal>(k_corner(x, y, z));

    // NOTE: A (small) hack to avoid one-sided differences. We need at
    // least 2 interior points due to an awkwardness with the
    // boundaries. The splines need derivatives in x, but we don't
    // know the value in the boundaries, so _any_ interpolation in the
    // last interior cell can't be done. Instead, we fudge the
    // interpolation in the last cell to be at the extreme right-hand
    // edge of the previous cell (that is, exactly on the last
    // interior point). However, this doesn't work with only one
    // interior point, because we have to do something similar to the
    // _first_ cell, and these two fudges cancel out and we end up
    // indexing into the boundary anyway.
    // TODO(peter): Can we remove this if we apply (dirichlet?) BCs to
    // the X derivatives? Note that we need at least _2_
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

    i_corner[i] =
        Ind3D(((((i_corn * ny) + (y + y_offset)) * nz) + k_corner(x, y, z)), ny, nz);

    h00_x[i] = (2. * t_x * t_x * t_x) - (3. * t_x * t_x) + 1.;
    h00_z[i] = (2. * t_z * t_z * t_z) - (3. * t_z * t_z) + 1.;

    h01_x[i] = (-2. * t_x * t_x * t_x) + (3. * t_x * t_x);
    h01_z[i] = (-2. * t_z * t_z * t_z) + (3. * t_z * t_z);

    h10_x[i] = t_x * (1. - t_x) * (1. - t_x);
    h10_z[i] = t_z * (1. - t_z) * (1. - t_z);

    h11_x[i] = (t_x * t_x * t_x) - (t_x * t_x);
    h11_z[i] = (t_z * t_z * t_z) - (t_z * t_z);

    std::array<BoutReal, 16> weights{};
    weights.fill(0.0);

    // The distribution of our weights:
    //  0   4   8    12
    //  1   5   9    13
    //  2   6   10   14
    //  3   7   11   15
    // e.g. 1 == ic.xm(); 4 == ic.zm(); 5 == ic;  7 == ic.zp(2);

    // f[ic] * h00_x[i] + f[icxp] * h01_x[i] + fx[ic] * h10_x[i] + fx[icxp] * h11_x[i];
    weights[5] += h00_x[i] * h00_z[i];
    weights[9] += h01_x[i] * h00_z[i];
    weights[9] += h10_x[i] * h00_z[i] / 2;
    weights[1] -= h10_x[i] * h00_z[i] / 2;
    weights[13] += h11_x[i] * h00_z[i] / 2;
    weights[5] -= h11_x[i] * h00_z[i] / 2;

    // f[iczp] * h00_x[i] + f[icxpzp] * h01_x[i] +
    // fx[iczp] * h10_x[i] + fx[icxpzp] * h11_x[i];
    weights[6] += h00_x[i] * h01_z[i];
    weights[10] += h01_x[i] * h01_z[i];
    weights[10] += h10_x[i] * h01_z[i] / 2;
    weights[2] -= h10_x[i] * h01_z[i] / 2;
    weights[14] += h11_x[i] * h01_z[i] / 2;
    weights[6] -= h11_x[i] * h01_z[i] / 2;

    // fz[ic] * h00_x[i] + fz[icxp] * h01_x[i] +
    // fxz[ic] * h10_x[i]+ fxz[icxp] * h11_x[i];
    weights[6] += h00_x[i] * h10_z[i] / 2;
    weights[4] -= h00_x[i] * h10_z[i] / 2;
    weights[10] += h01_x[i] * h10_z[i] / 2;
    weights[8] -= h01_x[i] * h10_z[i] / 2;
    weights[10] += h10_x[i] * h10_z[i] / 4;
    weights[8] -= h10_x[i] * h10_z[i] / 4;
    weights[2] -= h10_x[i] * h10_z[i] / 4;
    weights[0] += h10_x[i] * h10_z[i] / 4;
    weights[14] += h11_x[i] * h10_z[i] / 4;
    weights[12] -= h11_x[i] * h10_z[i] / 4;
    weights[6] -= h11_x[i] * h10_z[i] / 4;
    weights[4] += h11_x[i] * h10_z[i] / 4;

    // fz[iczp] * h00_x[i] + fz[icxpzp] * h01_x[i] +
    // fxz[iczp] * h10_x[i] + fxz[icxpzp] * h11_x[i];
    weights[7] += h00_x[i] * h11_z[i] / 2;
    weights[5] -= h00_x[i] * h11_z[i] / 2;
    weights[11] += h01_x[i] * h11_z[i] / 2;
    weights[9] -= h01_x[i] * h11_z[i] / 2;
    weights[11] += h10_x[i] * h11_z[i] / 4;
    weights[9] -= h10_x[i] * h11_z[i] / 4;
    weights[3] -= h10_x[i] * h11_z[i] / 4;
    weights[1] += h10_x[i] * h11_z[i] / 4;
    weights[15] += h11_x[i] * h11_z[i] / 4;
    weights[13] -= h11_x[i] * h11_z[i] / 4;
    weights[7] -= h11_x[i] * h11_z[i] / 4;
    weights[5] += h11_x[i] * h11_z[i] / 4;

    const PetscInt idxn = {conv.fromLocalToGlobal(x, y + y_offset, z)};
    constexpr std::size_t stencil_size = 4;
    for (std::size_t j = 0; j < stencil_size; ++j) {
      std::array<PetscInt, stencil_size> idxm{};
      std::array<PetscScalar, stencil_size> vals{};
      for (std::size_t k = 0; k < stencil_size; ++k) {
        idxm[k] = conv.fromMeshToGlobal(i_corn - 1 + static_cast<int>(j), y + y_offset,
                                        k_corner(x, y, z) - 1 + static_cast<int>(k));
        vals[k] = weights[(j * stencil_size) + k];
      }
      BOUT_OMP(critical)
      MatSetValues(petscWeights, 1, &idxn, stencil_size, idxm.data(), vals.data(),
                   INSERT_VALUES);
    }
  }

  MatAssemblyBegin(petscWeights, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(petscWeights, MAT_FINAL_ASSEMBLY);
  if (!isInit) {
    MatCreateVecs(petscWeights, &rhs, &result);
  }
  isInit = true;
}

void XZPetscHermiteSpline::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
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
XZPetscHermiteSpline::getWeightsForYApproximation(int i, int j, int k, int yoffset) {
  const int nz = localmesh->LocalNz;
  const int k_mod = k_corner(i, j, k);
  const int k_mod_m1 = (k_mod > 0) ? (k_mod - 1) : (nz - 1);
  const int k_mod_p1 = (k_mod == nz) ? 0 : k_mod + 1;
  const int k_mod_p2 = (k_mod_p1 == nz) ? 0 : k_mod_p1 + 1;

  return {{i, j + yoffset, k_mod_m1, -0.5 * h10_z(i, j, k)},
          {i, j + yoffset, k_mod, h00_z(i, j, k) - (0.5 * h11_z(i, j, k))},
          {i, j + yoffset, k_mod_p1, h01_z(i, j, k) + (0.5 * h10_z(i, j, k))},
          {i, j + yoffset, k_mod_p2, 0.5 * h11_z(i, j, k)}};
}

Field3D XZPetscHermiteSpline::interpolate(const Field3D& f,
                                          const std::string& region) const {

  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  const auto region2 =
      y_offset == 0 ? "RGN_NOY" : fmt::format("RGN_YPAR_{:+d}", y_offset);

  BoutReal* ptr = nullptr;
  const BoutReal* cptr = nullptr;
  VecGetArray(rhs, &ptr);
  BOUT_FOR(i, f.getRegion("RGN_NOY")) { ptr[int(i)] = f[i]; }
  VecRestoreArray(rhs, &ptr);
  MatMult(petscWeights, rhs, result);
  VecGetArrayRead(result, &cptr);
  BOUT_FOR(i, f.getRegion(region2)) {
    f_interp[i] = cptr[int(i)];
    ASSERT2(std::isfinite(cptr[int(i)]));
  }
  VecRestoreArrayRead(result, &cptr);

  f_interp.setRegion(region2);
  ASSERT2(f_interp.getRegionID());
  return f_interp;
}

Field3D XZPetscHermiteSpline::interpolate(const Field3D& f, const Field3D& delta_x,
                                          const Field3D& delta_z,
                                          const std::string& region) {
  calcWeights(delta_x, delta_z, region);
  return interpolate(f, region);
}

Field3D XZPetscHermiteSpline::interpolate(const Field3D& f, const Field3D& delta_x,
                                          const Field3D& delta_z, const BoutMask& mask,
                                          const std::string& region) {
  calcWeights(delta_x, delta_z, mask, region);
  return interpolate(f, region);
}

#endif // BOUT_HAS_PETSC
