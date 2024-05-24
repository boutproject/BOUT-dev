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

#include "../impls/bout/boutmesh.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/globalindexer.hxx"
#include "bout/globals.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/operatorstencil.hxx"
#include "bout/petsc_interface.hxx"

#include <array>
#include <cstddef>
#include <vector>

class IndConverter {
public:
  IndConverter(Mesh* mesh)
      : mesh(dynamic_cast<BoutMesh*>(mesh)), nxpe(mesh->getNXPE()),
        nype(mesh->getNYPE()) {}
  // ix and iy are global indices
  // iy is local
  int fromMeshToGlobal(int ix, int iy, int iz) {
    const int xstart = mesh->xstart;
    const int lnx = mesh->LocalNx - xstart * 2;
    // x-proc-id
    int pex = divToNeg(ix - xstart, lnx);
    if (pex < 0) {
      pex = 0;
    }
    if (pex >= nxpe) {
      pex = nxpe - 1;
    }
    const int zstart = 0;
    const int lnz = mesh->LocalNz - zstart * 2;
    // z-proc-id
    // pez only for wrapping around ; later needs similar treatment than pey
    const int pez = divToNeg(iz - zstart, lnz);
    // y proc-id - y is already local
    const int ystart = mesh->ystart;
    const int lny = mesh->LocalNy - ystart * 2;
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
    return fromLocalToGlobal(ix - pex * lnx, iy - pey_offset * lny, iz - pez * lnz, pex,
                             pey, 0);
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
    const int ilocal = ((ilocalx * mesh->LocalNy) + ilocaly) * mesh->LocalNz + ilocalz;
    const int ret = ilocal
                    + mesh->LocalNx * mesh->LocalNy * mesh->LocalNz
                          * ((pey * nxpe + pex) * nzpe + pez);
    ASSERT3(ret >= 0);
    ASSERT3(ret < nxpe * nype * mesh->LocalNx * mesh->LocalNy * mesh->LocalNz);
    return ret;
  }

private:
  // number of procs
  BoutMesh* mesh;
  int nxpe;
  int nype;
  int nzpe{1};
  static int divToNeg(const int n, const int d) {
    return (n < 0) ? ((n - d + 1) / d) : (n / d);
  }
};

namespace {
// A 4x4 stencil in x-z, covering (x-1, z-1) to (x+2, z+2)
auto XZstencil(const Mesh& localmesh) {
  using ind_type = Field3D::ind_type;
  OperatorStencil<ind_type> stencil;
  IndexOffset<ind_type> zero;

  // Stencil pattern
  //        x -->
  //    0   4   8    12
  // z  1   5   9    13
  // |  2   6   10   14
  // V  3   7   11   15

  const std::vector<IndexOffset<ind_type>> offsets = {
      zero.xm().zm(),      zero.xm(),      zero.xm().zp(),      zero.xm().zp().zp(),
      zero.zm(),           zero,           zero.zp(),           zero.zp().zp(),
      zero.xp().zm(),      zero.xp(),      zero.xp().zp(),      zero.xp().zp().zp(),
      zero.xp().xp().zm(), zero.xp().xp(), zero.xp().xp().zp(), zero.xp().xp().zp().zp(),
  };

  stencil.add(
      [&localmesh](ind_type ind) -> bool {
        return (localmesh.xstart <= ind.x() and ind.x() <= localmesh.xend)
               and (localmesh.ystart <= ind.y() and ind.y() <= localmesh.yend)
               and (localmesh.zstart <= ind.z() and ind.z() <= localmesh.zend);
      },
      offsets);
  stencil.add([]([[maybe_unused]] auto ind) -> bool { return true; }, {zero});

  return stencil;
}
} // namespace

PetscXZHermiteSpline::PetscXZHermiteSpline(int y_offset, Mesh* mesh)
    : XZInterpolation(y_offset, mesh),
      petsclib(&Options::root()["mesh:paralleltransform:xzinterpolation:hermitespline"]),
      indexer(std::make_shared<GlobalIndexer<Field3D>>(localmesh, XZstencil(*localmesh))),
      weights(indexer), h00_x(localmesh), h01_x(localmesh), h10_x(localmesh),
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

  newWeights.reserve(16);
  for (int w = 0; w < 16; ++w) {
    newWeights.emplace_back(localmesh);
    newWeights[w].allocate();
  }

  const int grid_size = localmesh->LocalNx * localmesh->LocalNy * localmesh->LocalNz;
  const int proc_size = grid_size * localmesh->getNXPE() * localmesh->getNYPE();
  MatCreateAIJ(MPI_COMM_WORLD, grid_size, grid_size, proc_size, proc_size, 16, nullptr,
               16, nullptr, &petscWeights);
}

PetscXZHermiteSpline::~PetscXZHermiteSpline() {
  if (!isInit) {
    return;
  }
  MatDestroy(&petscWeights);
  VecDestroy(&rhs);
  VecDestroy(&result);
}

void PetscXZHermiteSpline::calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                                       const std::string& region) {

  const int ny = localmesh->LocalNy;
  const int nz = localmesh->LocalNz;
  const int xend = (localmesh->xend - localmesh->xstart + 1) * localmesh->getNXPE()
                   + localmesh->xstart - 1;

  IndConverter conv{localmesh};

  const auto actual_region = getRegion(region);

  // PETSc doesn't like mixing `INSERT` and `ADD`, so first make sure matrix is zeroed
  // TODO: PetscMatrix could buffer and insert all at once
  weights.partialAssemble();
  BOUT_FOR(i, actual_region) { weights(i, i) = 0.0; }
  weights.partialAssemble();

  BOUT_FOR(i, actual_region) {
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

    for (int w = 0; w < 16; ++w) {
      newWeights[w][i] = 0;
    }
    // The distribution of our weights:
    //  0   4   8    12
    //  1   5   9    13
    //  2   6   10   14
    //  3   7   11   15
    // e.g. 1 == ic.xm(); 4 == ic.zm(); 5 == ic;  7 == ic.zp(2);

    // f[ic] * h00_x[i] + f[icxp] * h01_x[i] + fx[ic] * h10_x[i] + fx[icxp] * h11_x[i];
    weights(i, i) += (h00_x[i] * h00_z[i]) - (h11_x[i] * h00_z[i] / 2);
    weights(i, i.xp()) += (h01_x[i] * h00_z[i]) + (h10_x[i] * h00_z[i] / 2);
    weights(i, i.xm()) += -(h10_x[i] * h00_z[i] / 2);
    weights(i, i.xp(2)) += (h11_x[i] * h00_z[i] / 2);

    newWeights[5][i] += h00_x[i] * h00_z[i];
    newWeights[9][i] += h01_x[i] * h00_z[i];
    newWeights[9][i] += h10_x[i] * h00_z[i] / 2;
    newWeights[1][i] -= h10_x[i] * h00_z[i] / 2;
    newWeights[13][i] += h11_x[i] * h00_z[i] / 2;
    newWeights[5][i] -= h11_x[i] * h00_z[i] / 2;

    // f[iczp] * h00_x[i] + f[icxpzp] * h01_x[i] +
    // fx[iczp] * h10_x[i] + fx[icxpzp] * h11_x[i];
    weights(i, i.zp()) += (h00_x[i] * h01_z[i]) - (h11_x[i] * h01_z[i] / 2);
    weights(i, i.zp().xp()) += (h01_x[i] * h01_z[i]) + (h10_x[i] * h01_z[i] / 2);
    weights(i, i.zp().xm()) += -(h10_x[i] * h01_z[i] / 2);
    weights(i, i.zp().xp(2)) += h11_x[i] * h01_z[i] / 2;

    newWeights[6][i] += h00_x[i] * h01_z[i];
    newWeights[10][i] += h01_x[i] * h01_z[i];
    newWeights[10][i] += h10_x[i] * h01_z[i] / 2;
    newWeights[2][i] -= h10_x[i] * h01_z[i] / 2;
    newWeights[14][i] += h11_x[i] * h01_z[i] / 2;
    newWeights[6][i] -= h11_x[i] * h01_z[i] / 2;

    // fz[ic] * h00_x[i] + fz[icxp] * h01_x[i] +
    // fxz[ic] * h10_x[i]+ fxz[icxp] * h11_x[i];
    weights(i, i.zp()) += (h00_x[i] * h10_z[i] / 2) - (h11_x[i] * h10_z[i] / 4);
    weights(i, i.zm()) += -(h00_x[i] * h10_z[i] / 2) + (h11_x[i] * h10_z[i] / 4);
    weights(i, i.xp().zp()) += (h01_x[i] * h10_z[i] / 2) + (h10_x[i] * h10_z[i] / 4);
    weights(i, i.xp().zm()) += -(h01_x[i] * h10_z[i] / 2) - (h10_x[i] * h10_z[i] / 4);
    weights(i, i.xm().zp()) += -(h10_x[i] * h10_z[i] / 4);
    weights(i, i.xm().zm()) += (h10_x[i] * h10_z[i] / 4);
    weights(i, i.xp(2).zp()) += h11_x[i] * h10_z[i] / 4;
    weights(i, i.xp(2).zm()) += -(h11_x[i] * h10_z[i] / 4);

    newWeights[6][i] += h00_x[i] * h10_z[i] / 2;
    newWeights[4][i] -= h00_x[i] * h10_z[i] / 2;
    newWeights[10][i] += h01_x[i] * h10_z[i] / 2;
    newWeights[8][i] -= h01_x[i] * h10_z[i] / 2;
    newWeights[10][i] += h10_x[i] * h10_z[i] / 4;
    newWeights[8][i] -= h10_x[i] * h10_z[i] / 4;
    newWeights[2][i] -= h10_x[i] * h10_z[i] / 4;
    newWeights[0][i] += h10_x[i] * h10_z[i] / 4;
    newWeights[14][i] += h11_x[i] * h10_z[i] / 4;
    newWeights[12][i] -= h11_x[i] * h10_z[i] / 4;
    newWeights[6][i] -= h11_x[i] * h10_z[i] / 4;
    newWeights[4][i] += h11_x[i] * h10_z[i] / 4;

    // fz[iczp] * h00_x[i] + fz[icxpzp] * h01_x[i] +
    // fxz[iczp] * h10_x[i] + fxz[icxpzp] * h11_x[i];
    weights(i, i.zp(2)) += (h00_x[i] * h11_z[i] / 2) - (h11_x[i] * h11_z[i] / 4);
    weights(i, i) += -(h00_x[i] * h11_z[i] / 2) + (h11_x[i] * h11_z[i] / 4);
    weights(i, i.xp().zp(2)) += (h01_x[i] * h11_z[i] / 2) + (h10_x[i] * h11_z[i] / 4);
    weights(i, i.xp()) += -(h01_x[i] * h11_z[i] / 2) - (h10_x[i] * h11_z[i] / 4);
    weights(i, i.xm().zp(2)) += -(h10_x[i] * h11_z[i] / 4);
    weights(i, i.xm()) += (h10_x[i] * h11_z[i] / 4);
    weights(i, i.xp(2).zp(2)) += h11_x[i] * h11_z[i] / 4;
    weights(i, i.xp(2)) += -(h11_x[i] * h11_z[i] / 4);
    weights(i, i.zp(2)) += -(h11_x[i] * h11_z[i] / 4);

    newWeights[7][i] += h00_x[i] * h11_z[i] / 2;
    newWeights[5][i] -= h00_x[i] * h11_z[i] / 2;
    newWeights[11][i] += h01_x[i] * h11_z[i] / 2;
    newWeights[9][i] -= h01_x[i] * h11_z[i] / 2;
    newWeights[11][i] += h10_x[i] * h11_z[i] / 4;
    newWeights[9][i] -= h10_x[i] * h11_z[i] / 4;
    newWeights[3][i] -= h10_x[i] * h11_z[i] / 4;
    newWeights[1][i] += h10_x[i] * h11_z[i] / 4;
    newWeights[15][i] += h11_x[i] * h11_z[i] / 4;
    newWeights[13][i] -= h11_x[i] * h11_z[i] / 4;
    newWeights[7][i] -= h11_x[i] * h11_z[i] / 4;
    newWeights[5][i] += h11_x[i] * h11_z[i] / 4;

    std::array<PetscInt, 1> idxn = {conv.fromLocalToGlobal(x, y + y_offset, z)};
    for (int j = 0; j < 4; ++j) {
      std::array<PetscInt, 4> idxm{};
      std::array<PetscScalar, 4> vals{};
      for (std::size_t k = 0; k < 4; ++k) {
        idxm[k] = conv.fromMeshToGlobal(i_corn - 1 + j, y + y_offset,
                                        k_corner(x, y, z) - 1 + k);
        vals[k] = newWeights[j * 4 + k][i];
      }
      MatSetValues(petscWeights, 1, idxn.data(), 4, idxm.data(), vals.data(),
                   INSERT_VALUES);
    }
  }

  weights.assemble();

  MatAssemblyBegin(petscWeights, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(petscWeights, MAT_FINAL_ASSEMBLY);
  if (!isInit) {
    MatCreateVecs(petscWeights, &rhs, &result);
  }
  isInit = true;
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

Field3D PetscXZHermiteSpline::interpolate(const Field3D& f,
                                          const std::string& region) const {

  ASSERT1(f.getMesh() == localmesh);
  Field3D f_interp{emptyFrom(f)};

  BoutReal* ptr = nullptr;
  const BoutReal* cptr = nullptr;
  VecGetArray(rhs, &ptr);
  BOUT_FOR(i, f.getRegion("RGN_NOY")) { ptr[int(i)] = f[i]; }
  VecRestoreArray(rhs, &ptr);
  MatMult(petscWeights, rhs, result);
  VecGetArrayRead(result, &cptr);
  const auto region2 = y_offset == 0 ? region : fmt::format("RGN_YPAR_{:+d}", y_offset);
  BOUT_FOR(i, f.getRegion(region2)) {
    f_interp[i] = cptr[int(i)];
    ASSERT2(std::isfinite(cptr[int(i)]));
  }
  VecRestoreArrayRead(result, &cptr);

  PetscVector<Field3D> f_petsc(f, indexer);
  PetscVector<Field3D> f_result = weights * f_petsc;
  auto new_result = f_result.toField();
  // Using the Petsc interface classes gives pretty shoddy results
  // return new_result;
  return f_interp;
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
