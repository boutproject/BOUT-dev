/**************************************************************************
 * Various differential operators defined on BOUT grid
 *
 **************************************************************************
 * Copyright 2010 - 2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#include "bout/build_config.hxx"
#include "bout/build_defines.hxx"
#include "bout/dcomplex.hxx"
#include "bout/metric_tensor.hxx"
#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/coordinates.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/fft.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/globals.hxx>
#include <bout/interpolation.hxx>
#include <bout/invert_laplace.hxx> // Delp2 uses same coefficients as inversion code
#include <bout/mesh.hxx>
#include <bout/region.hxx>
#include <bout/solver.hxx>
#include <bout/unused.hxx>
#include <bout/utils.hxx>

#include <cmath>
#include <limits>
#include <string>

/*******************************************************************************
* Grad_par
* The parallel derivative along unperturbed B-field
*******************************************************************************/

bout::FieldMetric Grad_par(const Field2D& var, CELL_LOC outloc,
                           const std::string& method) {
  return DDY(var, outloc, method) * var.getCoordinates(outloc)->invSg();
}

Field3D Grad_par(const Field3DParallel& var, CELL_LOC outloc, const std::string& method) {
  return DDY(var, outloc, method) * var.getCoordinates(outloc)->invSg();
}

/*******************************************************************************
* Grad_parP
*
* Derivative along perturbed field-line
*
* b0 dot Grad  -  (1/B)b0 x Grad(apar) dot Grad
*
* Combines the parallel and perpendicular calculation to include
* grid-points at the corners.
*******************************************************************************/

Field3D Grad_parP(const Field3D& apar, const Field3D& f) {
  ASSERT1_FIELDS_COMPATIBLE(apar, f);
  ASSERT1(f.hasParallelSlices());

  const Mesh* mesh = apar.getMesh();

  Field3D result{emptyFrom(f)};

  const int ncz = mesh->LocalNz;

  const Coordinates* metric = apar.getCoordinates();

  Field3D gys{emptyFrom(f)};

  // Need Y derivative everywhere
  for (int x = 1; x <= mesh->LocalNx - 2; x++) {
    for (int y = 1; y <= mesh->LocalNy - 2; y++) {
      for (int z = 0; z < ncz; z++) {
        gys(x, y, z) = (f.yup()(x, y + 1, z) - f.ydown()(x, y - 1, z))
                       / ((0.5 * metric->dy(x, y + 1, z)) + metric->dy(x, y, z)
                          + (0.5 * metric->dy(x, y - 1, z)));
      }
    }
  }

  for (int x = 1; x <= mesh->LocalNx - 2; x++) {
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      for (int z = 0; z < ncz; z++) {
        const BoutReal by = 1. / sqrt(metric->g_22(x, y, z));
        // Z indices zm and zp
        const int zm = (z - 1 + ncz) % ncz;
        const int zp = (z + 1) % ncz;

        // bx = -DDZ(apar)
        const BoutReal bx = (apar(x, y, zm) - apar(x, y, zp))
                            / ((0.5 * metric->dz(x, y, zm)) + metric->dz(x, y, z)
                               + (0.5 * metric->dz(x, y, zp)));
        // bz = DDX(f)
        const BoutReal bz = (apar(x + 1, y, z) - apar(x - 1, y, z))
                            / ((0.5 * metric->dx(x - 1, y, z)) + metric->dx(x, y, z)
                               + (0.5 * metric->dx(x + 1, y, z)));

        // Now calculate (bx*d/dx + by*d/dy + bz*d/dz) f

        // Length dl for predictor
        BoutReal dl = fabs(metric->dx(x, y, z)) / (fabs(bx) + 1e-16);
        dl = BOUTMIN(dl, fabs(metric->dy(x, y, z)) / (fabs(by) + 1e-16));
        dl = BOUTMIN(dl, fabs(metric->dz(x, y, z)) / (fabs(bz) + 1e-16));

        BoutReal fp, fm;

        // X differencing
        fp =
            f(x + 1, y, z)
            + (0.25 * dl / metric->dz(x, y, z)) * bz * (f(x + 1, y, zm) - f(x + 1, y, zp))
            - 0.5 * dl * by * gys(x + 1, y, z);

        fm =
            f(x - 1, y, z)
            + (0.25 * dl / metric->dz(x, y, z)) * bz * (f(x - 1, y, zm) - f(x - 1, y, zp))
            - 0.5 * dl * by * gys(x - 1, y, z);

        result(x, y, z) = bx * (fp - fm)
                          / (0.5 * metric->dx(x - 1, y, z) + metric->dx(x, y, z)
                             + 0.5 * metric->dx(x + 1, y, z));

        // Z differencing

        fp =
            f(x, y, zp)
            + (0.25 * dl / metric->dx(x, y, z)) * bx * (f(x - 1, y, zp) - f(x + 1, y, zp))
            - 0.5 * dl * by * gys(x, y, zp);

        fm =
            f(x, y, zm)
            + (0.25 * dl / metric->dx(x, y, z)) * bx * (f(x - 1, y, zm) - f(x + 1, y, zm))
            - 0.5 * dl * by * gys(x, y, zm);

        result(x, y, z) += bz * (fp - fm)
                           / (0.5 * metric->dz(x, y, zm) + metric->dz(x, y, z)
                              + 0.5 * metric->dz(x, y, zp));

        // Y differencing

        fp = f.yup()(x, y + 1, z)
             - 0.5 * dl * bx * (f.yup()(x + 1, y + 1, z) - f.yup()(x - 1, y + 1, z))
                   / (0.5 * metric->dx(x - 1, y, z) + metric->dx(x, y, z)
                      + 0.5 * metric->dx(x + 1, y, z))

             + (0.25 * dl / metric->dz(x, y, z)) * bz
                   * (f.yup()(x, y + 1, zm) - f.yup()(x, y + 1, zp));

        fm = f.ydown()(x, y - 1, z)
             - 0.5 * dl * bx * (f.ydown()(x + 1, y - 1, z) - f.ydown()(x - 1, y - 1, z))
                   / (0.5 * metric->dx(x - 1, y, z) + metric->dx(x, y, z)
                      + 0.5 * metric->dx(x + 1, y, z))
             + (0.25 * dl / metric->dz(x, y, z)) * bz
                   * (f.ydown()(x, y - 1, zm) - f.ydown()(x, y - 1, zp));

        result(x, y, z) += by * (fp - fm)
                           / (0.5 * metric->dy(x, y - 1, z) + metric->dy(x, y, z)
                              + 0.5 * metric->dy(x, y + 1, z));
      }
    }
  }

  ASSERT2(result.getLocation() == f.getLocation());

  return result;
}

/*******************************************************************************
* Vpar_Grad_par
* vparallel times the parallel derivative along unperturbed B-field
*******************************************************************************/

bout::FieldMetric Vpar_Grad_par(const Field2D& v, const Field2D& f, CELL_LOC outloc,
                                const std::string& method) {
  return VDDY(v, f, outloc, method) * f.getCoordinates(outloc)->invSg();
}

Field3D Vpar_Grad_par(const Field3DParallel& v, const Field3DParallel& f, CELL_LOC outloc,
                      const std::string& method) {
  return VDDY(v, f, outloc, method) * f.getCoordinates(outloc)->invSg();
}

/*******************************************************************************
* Div_par
* parallel divergence operator B \partial_{||} (F/B)
*******************************************************************************/
bout::FieldMetric Div_par(const Field2D& f, CELL_LOC outloc, const std::string& method) {
  const auto& Bxy_outloc = f.getCoordinates(outloc)->Bxy();
  // Need Bxy at location of f, which might be different from outloc
  const auto& Bxy_floc = f.getCoordinates()->Bxy();

  return Bxy_outloc * Grad_par(bout::FieldMetric{f / Bxy_floc}, outloc, method);
}

Field3D Div_par(const Field3DParallel& f, CELL_LOC outloc, const std::string& method) {
  const auto& Bxy_outloc = f.getCoordinates(outloc)->Bxy();
  // Need Bxy at location of f, which might be different from outloc
  const auto& Bxy_floc = f.getCoordinates()->Bxy();

  return Bxy_outloc * Grad_par(f / Bxy_floc, outloc, method);
}

Field3D Div_par(const Field3DParallel& f, const Field3DParallel& v) {
  ASSERT1_FIELDS_COMPATIBLE(f, v);
  ASSERT1(f.hasParallelSlices());
  ASSERT1(v.hasParallelSlices());

  // Parallel divergence, using velocities at cell boundaries
  // Note: Not guaranteed to be flux conservative
  const Mesh* mesh = f.getMesh();

  Field3D result{emptyFrom(f)};

  const Coordinates* coord = f.getCoordinates();

  for (int i = mesh->xstart; i <= mesh->xend; i++) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = mesh->zstart; k <= mesh->zend; k++) {
        // Value of f and v at left cell face
        const BoutReal fL = 0.5 * (f(i, j, k) + f.ydown()(i, j - 1, k));
        const BoutReal vL = 0.5 * (v(i, j, k) + v.ydown()(i, j - 1, k));

        const BoutReal fR = 0.5 * (f(i, j, k) + f.yup()(i, j + 1, k));
        const BoutReal vR = 0.5 * (v(i, j, k) + v.yup()(i, j + 1, k));

        // Calculate flux at right boundary (y+1/2)
        const BoutReal fluxRight =
            fR * vR * (coord->J(i, j, k) + coord->J(i, j + 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j + 1, k)));

        // Calculate at left boundary (y-1/2)
        const BoutReal fluxLeft =
            fL * vL * (coord->J(i, j, k) + coord->J(i, j - 1, k))
            / (sqrt(coord->g_22(i, j, k)) + sqrt(coord->g_22(i, j - 1, k)));

        result(i, j, k) =
            (fluxRight - fluxLeft) / (coord->dy(i, j, k) * coord->J(i, j, k));
      }
    }
  }

  return result;
}

//////// Flux methods

Field3D Div_par_flux(const Field3DParallel& v, const Field3DParallel& f, CELL_LOC outloc,
                     const std::string& method) {
  const Coordinates* metric = f.getCoordinates(outloc);

  auto Bxy_floc = f.getCoordinates()->Bxy();

  if (!f.hasParallelSlices()) {
    const Field3D f_B = f / Bxy_floc;
    return metric->Bxy() * FDDY(v, f_B, outloc, method) / sqrt(metric->g_22());
  }

  // Need to modify yup and ydown fields
  Field3D f_B = f / Bxy_floc;
  // Distinct fields
  f_B.splitParallelSlices();
  f_B.yup() = f.yup() / Bxy_floc;
  f_B.ydown() = f.ydown() / Bxy_floc;
  return metric->Bxy() * FDDY(v, f_B, outloc, method) / sqrt(metric->g_22());
}

/*******************************************************************************
* Grad2_par2
* second parallel derivative
*
* (b dot Grad)(b dot Grad)
*
* Note: For parallel Laplacian use LaplacePar
*******************************************************************************/

bout::FieldMetric Grad2_par2(const Field2D& f, CELL_LOC outloc,
                             const std::string& method) {
  const auto& coords = *f.getCoordinates(outloc);

  return coords.Grad2_par2_DDY_invSg(outloc, method) * DDY(f, outloc, method)
         + D2DY2(f, outloc, method) / coords.g_22();
}

Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc, const std::string& method) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  const auto& coords = *f.getCoordinates(outloc);

  return coords.Grad2_par2_DDY_invSg(outloc, method) * DDY(f, outloc, method)
         + D2DY2(f, outloc, method) / coords.g_22();
}

/*******************************************************************************
* Div_par_K_Grad_par
* Parallel divergence of diffusive flux, K*Grad_par
*******************************************************************************/

bout::FieldMetric Div_par_K_Grad_par(BoutReal kY, const Field2D& f, CELL_LOC outloc) {
  return kY * Grad2_par2(f, outloc);
}

Field3D Div_par_K_Grad_par(BoutReal kY, const Field3DParallel& f, CELL_LOC outloc) {
  return kY * Grad2_par2(f, outloc);
}

bout::FieldMetric Div_par_K_Grad_par(const Field2D& kY, const Field2D& f,
                                     CELL_LOC outloc) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return interp_to(kY, outloc) * Grad2_par2(f, outloc)
         + Div_par(kY, outloc) * Grad_par(f, outloc);
}

Field3D Div_par_K_Grad_par(const Field2D& kY, const Field3DParallel& f, CELL_LOC outloc) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return interp_to(kY, outloc) * Grad2_par2(f, outloc)
         + Div_par(kY, outloc) * Grad_par(f, outloc);
}

Field3D Div_par_K_Grad_par(const Field3DParallel& kY, const Field2D& f, CELL_LOC outloc) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return interp_to(kY, outloc) * Grad2_par2(f, outloc)
         + Div_par(kY, outloc) * Grad_par(f, outloc);
}

Field3D Div_par_K_Grad_par(const Field3DParallel& kY, const Field3DParallel& f,
                           CELL_LOC outloc) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return interp_to(kY, outloc) * Grad2_par2(f, outloc)
         + Div_par(kY, outloc) * Grad_par(f, outloc);
}

namespace {
template <bout::ConductionMethod conduction_method>
Field3D Div_par_K_Grad_par_mod_impl(const Field3DParallel& Kin,
                                    const Field3DParallel& fin, Field3D& flow_ylow,
                                    bool bndry_flux) {
  ASSERT2(Kin.getLocation() == fin.getLocation());

  const Mesh* mesh = Kin.getMesh();
  const Coordinates* coord = fin.getCoordinates();

  if (Kin.isFci()) {
    ASSERT1(Kin.hasParallelSlices());
    ASSERT1(fin.hasParallelSlices());
    // Using parallel slices.
    // Note: Y slices may use different coordinate systems
    //       -> Only B, dy and g_22 can be used in yup/ydown
    //          Others (e.g J) may not be averaged between y planes.

    const auto& K_up = Kin.yup();
    const auto& K_down = Kin.ydown();

    const auto& f_up = fin.yup();
    const auto& f_down = fin.ydown();

    Field3D result{zeroFrom(fin)};
    flow_ylow = zeroFrom(fin);

    BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
      const auto iyp = i.yp();
      const auto iym = i.ym();

      // Upper cell edge
      const BoutReal c_up = 0.5 * (Kin[i] + K_up[iyp]); // K at the upper boundary
      const BoutReal J_up =
          0.5 * (coord->J()[i] + coord->J().yup()[iyp]); // Jacobian at boundary
      const BoutReal g_22_up = 0.5 * (coord->g_22()[i] + coord->g_22().yup()[iyp]);
      const BoutReal gradient_up =
          2. * (f_up[iyp] - fin[i]) / (coord->dy()[i] + coord->dy().yup()[iyp]);

      const BoutReal flux_up = c_up * J_up * gradient_up / g_22_up;

      // Lower cell edge
      const BoutReal c_down = 0.5 * (Kin[i] + K_down[iym]); // K at the lower boundary
      const BoutReal J_down =
          0.5 * (coord->J()[i] + coord->J().ydown()[iym]); // Jacobian at boundary
      const BoutReal g_22_down = 0.5 * (coord->g_22()[i] + coord->g_22().ydown()[iym]);
      const BoutReal gradient_down =
          2. * (fin[i] - f_down[iym]) / (coord->dy()[i] + coord->dy().ydown()[iym]);

      const BoutReal flux_down = c_down * J_down * gradient_down / g_22_down;

      result[i] = (flux_up - flux_down) / (coord->dy()[i] * coord->J()[i]);
    }

    return result;
  }

  // Calculate in field-aligned coordinates
  const auto& K = toFieldAligned(Kin, "RGN_NOX");
  const auto& f = toFieldAligned(fin, "RGN_NOX");

  Field3D result{zeroFrom(f)};
  flow_ylow = zeroFrom(f);

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    // Calculate flux at upper surface

    const auto iyp = i.yp();
    const auto iym = i.ym();

    if (bndry_flux || mesh->periodicY(i.x()) || !mesh->lastY(i.x())
        || (i.y() != mesh->yend)) {
      BoutReal flux = 0.0;

      if constexpr (conduction_method == bout::ConductionMethod::Original) {
        const BoutReal c = 0.5 * (K[i] + K[iyp]); // K at the upper boundary
        const BoutReal J =
            0.5 * (coord->J()[i] + coord->J()[iyp]); // Jacobian at boundary
        const BoutReal g_22 = 0.5 * (coord->g_22()[i] + coord->g_22()[iyp]);

        const BoutReal gradient =
            2. * (f[iyp] - f[i]) / (coord->dy()[i] + coord->dy()[iyp]);

        flux = c * J * gradient / g_22;
      } else if constexpr (conduction_method == bout::ConductionMethod::ProductJK) {
        // Intended to reduce sensitivity of result to K in small cells
        const BoutReal cJ =
            0.5 * (K[i] * coord->J()[i] + K[iyp] * coord->J()[iyp]); // K * J at boundary
        const BoutReal g_22 = 0.5 * (coord->g_22()[i] + coord->g_22()[iyp]);

        const BoutReal gradient =
            2. * (f[iyp] - f[i]) / (coord->dy()[i] + coord->dy()[iyp]);

        flux = cJ * gradient / g_22;
      } else if constexpr (conduction_method == bout::ConductionMethod::Harmonic) {
        // Harmonic average (serial resistance)
        const BoutReal cond_i =
            K[i] * coord->J()[i] / (coord->g_22()[i] * coord->dy()[i]);
        const BoutReal cond_iyp =
            K[iyp] * coord->J()[iyp] / (coord->g_22()[iyp] * coord->dy()[iyp]);
        const BoutReal denom = cond_i + cond_iyp;

        // Harmonic mean: series resistance of two half-cells
        const BoutReal C_edge =
            (std::abs(denom) > std::numeric_limits<BoutReal>::epsilon())
                ? 2.0 * cond_i * cond_iyp / denom
                : 0.0;

        flux = C_edge * (f[iyp] - f[i]);
      }
      result[i] += flux / (coord->dy()[i] * coord->J()[i]);
    }

    // Calculate flux at lower surface
    if (bndry_flux || mesh->periodicY(i.x()) || !mesh->firstY(i.x())
        || (i.y() != mesh->ystart)) {
      BoutReal flux = 0.0;

      if constexpr (conduction_method == bout::ConductionMethod::Original) {
        const BoutReal c = 0.5 * (K[i] + K[iym]); // K at the lower boundary
        const BoutReal J =
            0.5 * (coord->J()[i] + coord->J()[iym]); // Jacobian at boundary
        const BoutReal g_22 = 0.5 * (coord->g_22()[i] + coord->g_22()[iym]);

        const BoutReal gradient =
            2. * (f[i] - f[iym]) / (coord->dy()[i] + coord->dy()[iym]);

        flux = c * J * gradient / g_22;
      } else if constexpr (conduction_method == bout::ConductionMethod::ProductJK) {
        const BoutReal cJ =
            0.5 * (K[i] * coord->J()[i] + K[iym] * coord->J()[iym]); // K * J at boundary
        const BoutReal g_22 = 0.5 * (coord->g_22()[i] + coord->g_22()[iym]);

        const BoutReal gradient =
            2. * (f[i] - f[iym]) / (coord->dy()[i] + coord->dy()[iym]);

        flux = cJ * gradient / g_22;
      } else if constexpr (conduction_method == bout::ConductionMethod::Harmonic) {
        const BoutReal cond_i =
            K[i] * coord->J()[i] / (coord->g_22()[i] * coord->dy()[i]);
        const BoutReal cond_iym =
            K[iym] * coord->J()[iym] / (coord->g_22()[iym] * coord->dy()[iym]);
        const BoutReal denom = cond_i + cond_iym;

        const BoutReal C_edge =
            (std::abs(denom) > std::numeric_limits<BoutReal>::epsilon())
                ? 2.0 * cond_i * cond_iym / denom
                : 0.0;

        flux = C_edge * (f[i] - f[iym]);
      }

      result[i] -= flux / (coord->dy()[i] * coord->J()[i]);
      flow_ylow[i] = -flux * coord->dx()[i] * coord->dz()[i];
    }
  }

  // Shifted to field aligned coordinates, so need to shift back
  result = fromFieldAligned(result, "RGN_NOBNDRY");
  flow_ylow = fromFieldAligned(flow_ylow);

  return result;
}
} // namespace

Field3D Div_par_K_Grad_par_mod(const Field3DParallel& Kin, const Field3DParallel& fin,
                               Field3D& flow_ylow, bool bndry_flux,
                               bout::ConductionMethod method) {
  using enum bout::ConductionMethod;

  switch (method) {
  case Original:
    return Div_par_K_Grad_par_mod_impl<Original>(Kin, fin, flow_ylow, bndry_flux);
  case ProductJK:
    return Div_par_K_Grad_par_mod_impl<ProductJK>(Kin, fin, flow_ylow, bndry_flux);
  case Harmonic:
    return Div_par_K_Grad_par_mod_impl<Harmonic>(Kin, fin, flow_ylow, bndry_flux);
  }
  throw BoutException(
      "Unknown method `{}` - choose from `Original`, `ProductJK` or `Harmonic`.",
      toString(method));
}
/*******************************************************************************
* Delp2
* perpendicular Laplacian operator
*******************************************************************************/

bout::FieldMetric Delp2(const Field2D& f, CELL_LOC outloc, [[maybe_unused]] bool useFFT) {
  const auto& coords = *f.getCoordinates(outloc);
  return coords.G1() * DDX(f, outloc) + coords.g11() * D2DX2(f, outloc);
}

Field3D Delp2(const Field3D& f, CELL_LOC outloc, bool useFFT) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  ASSERT1(f.getLocation() == outloc);
  const auto* mesh = f.getMesh();

  if (mesh->GlobalNx == 1 && mesh->GlobalNz == 1) {
    // copy mesh, location, etc
    return f * 0;
  }
  ASSERT2(mesh->xstart > 0); // Need at least one guard cell;

  Field3D result{emptyFrom(f).setLocation(outloc)};

  if (useFFT and not bout::build::use_metric_3d and mesh->getNZPE() == 1) {
    const int ncz = mesh->LocalNz;

    // Allocate memory
    auto ft = Matrix<dcomplex>(mesh->LocalNx, (ncz / 2) + 1);
    auto delft = Matrix<dcomplex>(mesh->LocalNx, (ncz / 2) + 1);

    // Loop over y indices
    // Note: should not include y-guard or y-boundary points here as that would
    // use values from corner cells in dx, which may not be initialised.
    for (int jy = mesh->ystart; jy <= mesh->yend; jy++) {

      // Take forward FFT

      for (int jx = 0; jx < mesh->LocalNx; jx++) {
        rfft(&f(jx, jy, 0), ncz, &ft(jx, 0));
      }

      // Loop over kz
      for (int jz = 0; jz <= ncz / 2; jz++) {

        // No smoothing in the x direction
        for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
          // Perform x derivative

          dcomplex a;
          dcomplex b;
          dcomplex c;
          laplace_tridag_coefs(jx, jy, jz, a, b, c, nullptr, nullptr, outloc);

          delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
        }
      }

      // Reverse FFT
      for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {

        irfft(&delft(jx, 0), ncz, &result(jx, jy, 0));
      }
    }
  } else {
    const auto& coords = *f.getCoordinates(outloc);
    result = coords.G1() * DDX(f, outloc) + coords.G3() * DDZ(f, outloc)
             + coords.g11() * D2DX2(f, outloc) + coords.g33() * D2DZ2(f, outloc)
             + 2 * coords.g13() * D2DXDZ(f, outloc);
  }

  ASSERT2(result.getLocation() == outloc);

  return result;
}

FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc, bool useFFT) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  ASSERT1(f.getLocation() == outloc);
  const auto* mesh = f.getMesh();

  if (mesh->GlobalNx == 1 && mesh->GlobalNz == 1) {
    // copy mesh, location, etc
    return f * 0;
  }
  ASSERT2(mesh->xstart > 0); // Need at least one guard cell

  FieldPerp result{emptyFrom(f).setLocation(outloc)};

  const int jy = f.getIndex();
  result.setIndex(jy);

  if (useFFT and mesh->getNZPE() == 1) {
    const int ncz = mesh->LocalNz;

    // Allocate memory
    auto ft = Matrix<dcomplex>(mesh->LocalNx, (ncz / 2) + 1);
    auto delft = Matrix<dcomplex>(mesh->LocalNx, (ncz / 2) + 1);

    // Take forward FFT
    for (int jx = 0; jx < mesh->LocalNx; jx++) {
      rfft(&f(jx, 0), ncz, &ft(jx, 0));
    }

    // Loop over kz
    for (int jz = 0; jz <= ncz / 2; jz++) {

      // No smoothing in the x direction
      for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
        // Perform x derivative

        dcomplex a;
        dcomplex b;
        dcomplex c;
        laplace_tridag_coefs(jx, jy, jz, a, b, c);

        delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
      }
    }

    // Reverse FFT
    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
      irfft(&delft(jx, 0), ncz, &result(jx, 0));
    }

  } else {
    throw BoutException("Non-fourier Delp2 not currently implented for FieldPerp.");
    // Would be the following but don't have standard derivative operators for FieldPerps
    // yet
    // result = G1 * ::DDX(f, outloc) + G3 * ::DDZ(f, outloc) + g11 * ::D2DX2(f, outloc)
    //          + g33 * ::D2DZ2(f, outloc) + 2 * g13 * ::D2DXDZ(f, outloc);
  };

  return result;
}

/*******************************************************************************
* LaplacePerp
* Full perpendicular Laplacian operator on scalar field
*
* Laplace_perp = Laplace - Laplace_par
*******************************************************************************/

bout::FieldMetric Laplace_perp(const Field2D& f, CELL_LOC outloc,
                               const std::string& dfdy_boundary_condition,
                               const std::string& dfdy_region) {
  return Laplace(f, outloc, dfdy_boundary_condition, dfdy_region)
         - Laplace_par(f, outloc);
}

Field3D Laplace_perp(const Field3D& f, CELL_LOC outloc,
                     const std::string& dfdy_boundary_condition,
                     const std::string& dfdy_region) {
  return Laplace(f, outloc, dfdy_boundary_condition, dfdy_region)
         - Laplace_par(f, outloc);
}

/*******************************************************************************
 * LaplacePar
 * Full parallel Laplacian operator on scalar field
 *
 * LaplacePar(f) = Div( b (b dot Grad(f)) )
 *
 *******************************************************************************/

bout::FieldMetric Laplace_par(const Field2D& f, CELL_LOC outloc) {
  const auto& coords = *f.getCoordinates(outloc);
  return D2DY2(f, outloc) / coords.g_22()
         + DDY(bout::FieldMetric{coords.J() / coords.g_22()}, outloc) * DDY(f, outloc)
               / coords.J();
}

Field3D Laplace_par(const Field3D& f, CELL_LOC outloc) {
  const auto& coords = *f.getCoordinates(outloc);
  return D2DY2(f, outloc) / coords.g_22()
         + DDY(coords.J().asField3DParallel() / coords.g_22(), outloc) * DDY(f, outloc)
               / coords.J();
}

/*******************************************************************************
* Laplacian
* Full Laplacian operator on scalar field
*******************************************************************************/

bout::FieldMetric Laplace(const Field2D& f, CELL_LOC outloc,
                          const std::string& dfdy_boundary_condition,
                          const std::string& dfdy_region) {
  const auto& coords = *f.getCoordinates(outloc);

  return coords.G1() * DDX(f, outloc) + coords.G2() * DDY(f, outloc)
         + coords.g11() * D2DX2(f, outloc) + coords.g22() * D2DY2(f, outloc)
         + 2.0 * coords.g12()
               * D2DXDY(f, outloc, "DEFAULT", "RGN_NOBNDRY", dfdy_boundary_condition,
                        dfdy_region);
}

Field3D Laplace(const Field3D& f, CELL_LOC outloc,
                const std::string& dfdy_boundary_condition,
                const std::string& dfdy_region) {
  const auto& coords = *f.getCoordinates(outloc);

  return coords.G1() * DDX(f, outloc) + coords.G2() * DDY(f, outloc)
         + coords.G3() * DDZ(f, outloc) + coords.g11() * D2DX2(f, outloc)
         + coords.g22() * D2DY2(f, outloc) + coords.g33() * D2DZ2(f, outloc)
         + 2.0
               * (coords.g12()
                      * D2DXDY(f, outloc, "DEFAULT", "RGN_NOBNDRY",
                               dfdy_boundary_condition, dfdy_region)
                  + coords.g13() * D2DXDZ(f, outloc) + coords.g23() * D2DYDZ(f, outloc));
}

/*******************************************************************************
 * Laplace_perpXY
 * Inverse of Laplacian operator in LaplaceXY solver
 *******************************************************************************/

Field2D Laplace_perpXY(const Field2D& A, const Field2D& f) {
#if BOUT_USE_METRIC_3D
  throw BoutException("Coordinates::Laplace_perpXY for 3D metric not implemented");
#else
  const auto& coords = *f.getCoordinates();

  Field2D result;
  result.allocate();
  for (auto i : result.getRegion(RGN_NOBNDRY)) {
    result[i] = 0.;

    // outer x boundary
    const auto outer_x_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.xp()]); };
    const BoutReal outer_x_A = outer_x_avg(A);
    const BoutReal outer_x_J = outer_x_avg(coords.J());
    const BoutReal outer_x_g11 = outer_x_avg(coords.g11());
    const BoutReal outer_x_dx = outer_x_avg(coords.dx());
    const BoutReal outer_x_value = outer_x_A * outer_x_J * outer_x_g11
                                   / (coords.J()[i] * outer_x_dx * coords.dx()[i]);
    result[i] += outer_x_value * (f[i.xp()] - f[i]);

    // inner x boundary
    const auto inner_x_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.xm()]); };
    const BoutReal inner_x_A = inner_x_avg(A);
    const BoutReal inner_x_J = inner_x_avg(coords.J());
    const BoutReal inner_x_g11 = inner_x_avg(coords.g11());
    const BoutReal inner_x_dx = inner_x_avg(coords.dx());
    const BoutReal inner_x_value = inner_x_A * inner_x_J * inner_x_g11
                                   / (coords.J()[i] * inner_x_dx * coords.dx()[i]);
    result[i] += inner_x_value * (f[i.xm()] - f[i]);

    // upper y boundary
    const auto upper_y_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.yp()]); };
    const BoutReal upper_y_A = upper_y_avg(A);
    const BoutReal upper_y_J = upper_y_avg(coords.J());
    const BoutReal upper_y_g_22 = upper_y_avg(coords.g_22());
    const BoutReal upper_y_g23 = upper_y_avg(coords.g23());
    const BoutReal upper_y_g_23 = upper_y_avg(coords.g_23());
    const BoutReal upper_y_dy = upper_y_avg(coords.dy());
    const BoutReal upper_y_value =
        -upper_y_A * upper_y_J * upper_y_g23 * upper_y_g_23
        / (upper_y_g_22 * coords.J()[i] * upper_y_dy * coords.dy()[i]);
    result[i] += upper_y_value * (f[i.yp()] - f[i]);

    // lower y boundary
    const auto lower_y_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.ym()]); };
    const BoutReal lower_y_A = lower_y_avg(A);
    const BoutReal lower_y_J = lower_y_avg(coords.J());
    const BoutReal lower_y_g_22 = lower_y_avg(coords.g_22());
    const BoutReal lower_y_g23 = lower_y_avg(coords.g23());
    const BoutReal lower_y_g_23 = lower_y_avg(coords.g_23());
    const BoutReal lower_y_dy = lower_y_avg(coords.dy());
    const BoutReal lower_y_value =
        -lower_y_A * lower_y_J * lower_y_g23 * lower_y_g_23
        / (lower_y_g_22 * coords.J()[i] * lower_y_dy * coords.dy()[i]);
    result[i] += lower_y_value * (f[i.ym()] - f[i]);
  }

  return result;
#endif
}

/*******************************************************************************
* b0xGrad_dot_Grad
* Terms of form b0 x Grad(phi) dot Grad(A)
* Used for ExB terms and perturbed B field using A_||
*******************************************************************************/

bout::FieldMetric b0xGrad_dot_Grad(const Field2D& phi, const Field2D& A,
                                   CELL_LOC outloc) {

  if (outloc == CELL_DEFAULT) {
    outloc = A.getLocation();
  }

  ASSERT1(phi.getMesh() == A.getMesh());

  const Coordinates* metric = phi.getCoordinates(outloc);

  // Calculate phi derivatives
  const bout::FieldMetric dpdx = DDX(phi, outloc);
  const bout::FieldMetric dpdy = DDY(phi, outloc);

  // Calculate advection velocity
  const bout::FieldMetric vx = -metric->g_23() * dpdy;
  const bout::FieldMetric vy = metric->g_23() * dpdx;

  // Upwind A using these velocities
  bout::FieldMetric result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc);
  result /= metric->J() * sqrt(metric->g_22());

  ASSERT1(result.getLocation() == outloc);

#if BOUT_USE_TRACK
  result.name = "b0xGrad_dot_Grad(" + phi.name + "," + A.name + ")";
#endif
  return result;
}

Field3D b0xGrad_dot_Grad(const Field2D& phi, const Field3D& A, CELL_LOC outloc) {

  if (outloc == CELL_DEFAULT) {
    outloc = A.getLocation();
  }

  ASSERT1(phi.getMesh() == A.getMesh());

  const Mesh* mesh = phi.getMesh();

  const Coordinates* metric = phi.getCoordinates(outloc);

  // Calculate phi derivatives
  const bout::FieldMetric dpdx = DDX(phi, outloc);
  const bout::FieldMetric dpdy = DDY(phi, outloc);

  // Calculate advection velocity
  const bout::FieldMetric vx = -metric->g_23() * dpdy;
  const bout::FieldMetric vy = metric->g_23() * dpdx;
  bout::FieldMetric vz = metric->g_12() * dpdy - metric->g_22() * dpdx;

  if (mesh->IncIntShear) {
    // BOUT-06 style differencing
    vz += metric->IntShiftTorsion() * vx;
  }

  // Upwind A using these velocities

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc) + VDDZ(vz, A, outloc);

  result /= (metric->J() * sqrt(metric->g_22()));

#if BOUT_USE_TRACK
  result.name = "b0xGrad_dot_Grad(" + phi.name + "," + A.name + ")";
#endif

  ASSERT2(result.getLocation() == outloc);

  return result;
}

Field3D b0xGrad_dot_Grad(const Field3D& p, const Field2D& A, CELL_LOC outloc) {

  if (outloc == CELL_DEFAULT) {
    outloc = A.getLocation();
  }

  ASSERT1(p.getMesh() == A.getMesh());

  const Coordinates* metric = p.getCoordinates(outloc);

  // Calculate phi derivatives
  const Field3D dpdx = DDX(p, outloc);
  const Field3D dpdy = DDY(p, outloc);
  const Field3D dpdz = DDZ(p, outloc);

  // Calculate advection velocity
  const Field3D vx = metric->g_22() * dpdz - metric->g_23() * dpdy;
  const Field3D vy = metric->g_23() * dpdx - metric->g_12() * dpdz;

  // Upwind A using these velocities

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc);

  result /= (metric->J() * sqrt(metric->g_22()));

#if BOUT_USE_TRACK
  result.name = "b0xGrad_dot_Grad(" + p.name + "," + A.name + ")";
#endif

  ASSERT2(result.getLocation() == outloc);

  return result;
}

Field3D b0xGrad_dot_Grad(const Field3D& phi, const Field3D& A, CELL_LOC outloc) {

  if (outloc == CELL_DEFAULT) {
    outloc = A.getLocation();
  }

  ASSERT1(phi.getMesh() == A.getMesh());

  const Mesh* mesh = phi.getMesh();

  const Coordinates* metric = phi.getCoordinates(outloc);

  // Calculate phi derivatives
  const Field3D dpdx = DDX(phi, outloc);
  const Field3D dpdy = DDY(phi, outloc);
  const Field3D dpdz = DDZ(phi, outloc);

  // Calculate advection velocity
  const Field3D vx = metric->g_22() * dpdz - metric->g_23() * dpdy;
  const Field3D vy = metric->g_23() * dpdx - metric->g_12() * dpdz;
  Field3D vz = metric->g_12() * dpdy - metric->g_22() * dpdx;

  if (mesh->IncIntShear) {
    // BOUT-06 style differencing
    vz += metric->IntShiftTorsion() * vx;
  }

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc) + VDDZ(vz, A, outloc);

  result /= (metric->J() * sqrt(metric->g_22()));

#if BOUT_USE_TRACK
  result.name = "b0xGrad_dot_Grad(" + phi.name + "," + A.name + ")";
#endif

  ASSERT2(result.getLocation() == outloc);

  return result;
}

/*******************************************************************************
 * Poisson bracket
 * Terms of form b0 x Grad(f) dot Grad(g) / B = [f, g]
 *******************************************************************************/

bout::FieldMetric bracket(const Field2D& f, const Field2D& g, BRACKET_METHOD method,
                          CELL_LOC outloc, Solver* UNUSED(solver)) {

  ASSERT1_FIELDS_COMPATIBLE(f, g);
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation());

  bout::FieldMetric result{emptyFrom(f)};

  if ((method == BRACKET_SIMPLE) || (method == BRACKET_ARAKAWA)) {
    // Use a subset of terms for comparison to BOUT-06
    result = 0.0;
    result.setLocation(outloc);
  } else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g, outloc) / f.getCoordinates(outloc)->Bxy();
  }
  return result;
}

Field3D bracket(const Field3D& f, const Field2D& g, BRACKET_METHOD method,
                CELL_LOC outloc, Solver* solver) {

  ASSERT1_FIELDS_COMPATIBLE(f, g);
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation());

  [[maybe_unused]] const Mesh* mesh = f.getMesh();

  Field3D result{emptyFrom(f).setLocation(outloc)};

  const Coordinates* metric = f.getCoordinates(outloc);

  switch (method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)

    if (solver == nullptr) {
      throw BoutException("CTU method requires access to the solver");
    }

#if not(BOUT_USE_METRIC_3D)
    const int ncz = mesh->LocalNz;

    for (int x = mesh->xstart; x <= mesh->xend; x++) {
      for (int y = mesh->ystart; y <= mesh->yend; y++) {
        for (int z = 0; z < ncz; z++) {
          const int zm = (z - 1 + ncz) % ncz;
          const int zp = (z + 1) % ncz;

          BoutReal gp;
          BoutReal gm;

          // Vx = DDZ(f)
          const BoutReal vx = (f(x, y, zp) - f(x, y, zm)) / (2. * metric->dz(x, y, z));

          // Set stability condition
          solver->setMaxTimestep(metric->dx(x, y, z) / (fabs(vx) + 1e-16));

          // X differencing
          if (vx > 0.0) {
            gp = g(x, y);

            gm = g(x - 1, y);

          } else {
            gp = g(x + 1, y);

            gm = g(x, y);
          }

          result(x, y, z) = vx * (gp - gm) / metric->dx(x, y, z);
        }
      }
    }
#else
    throw BoutException("BRACKET_CTU not valid with 3D metrics yet.");
#endif
    break;
  }
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow. Here as a test

#if not(BOUT_USE_METRIC_3D)
    const int ncz = mesh->LocalNz;

    BOUT_FOR(j2D, result.getRegion2D("RGN_NOBNDRY")) {
      // Get constants for this iteration
      const BoutReal spacingFactor = 1.0 / (12 * metric->dz()[j2D] * metric->dx()[j2D]);
      const int jy = j2D.y();
      const int jx = j2D.x();
      const int xm = jx - 1;
      const int xp = jx + 1;

      // Extract relevant Field2D values
      const BoutReal gxm = g(xm, jy);
      const BoutReal gc = g(jx, jy);
      const BoutReal gxp = g(xp, jy);

      // Index Field3D as 2D to get start of z data block
      const auto fxm = f(xm, jy);
      const auto fc = f(jx, jy);
      const auto fxp = f(xp, jy);

      // Here we split the loop over z into three parts; the first value, the middle block
      // and the last value
      // this is to allow the loop used in the middle block to vectorise.

      // The first value
      {
        const int jzp = 1;
        const int jzm = ncz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = 2 * (fc[jzp] - fc[jzm]) * (gxp - gxm);

        // J+x
        const BoutReal Jpx = (gxp * (fxp[jzp] - fxp[jzm])) - (gxm * (fxm[jzp] - fxm[jzm]))
                             + (gc * (fxp[jzm] - fxp[jzp] - fxm[jzm] + fxm[jzp]));

        result(jx, jy, 0) = (Jpp + Jpx) * spacingFactor;
      }

      // The middle block
      for (int jz = 1; jz < mesh->LocalNz - 1; jz++) {
        const int jzp = jz + 1;
        const int jzm = jz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = 2 * (fc[jzp] - fc[jzm]) * (gxp - gxm);

        // J+x
        const BoutReal Jpx = (gxp * (fxp[jzp] - fxp[jzm])) - (gxm * (fxm[jzp] - fxm[jzm]))
                             + (gc * (fxp[jzm] - fxp[jzp] - fxm[jzm] + fxm[jzp]));

        result(jx, jy, jz) = (Jpp + Jpx) * spacingFactor;
      }

      // The last value
      {
        const int jzp = 0;
        const int jzm = ncz - 2;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = 2 * (fc[jzp] - fc[jzm]) * (gxp - gxm);

        // J+x
        const BoutReal Jpx = (gxp * (fxp[jzp] - fxp[jzm])) - (gxm * (fxm[jzp] - fxm[jzm]))
                             + (gc * (fxp[jzm] - fxp[jzp] - fxm[jzm] + fxm[jzp]));

        result(jx, jy, ncz - 1) = (Jpp + Jpx) * spacingFactor;
      }
    }
#else
    throw BoutException("BRACKET_ARAKAWA not valid with 3D metrics yet.");
#endif

    break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f, outloc), g, outloc);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g, outloc) / metric->Bxy();
  }
  }
  return result;
}

Field3D bracket(const Field2D& f, const Field3D& g, BRACKET_METHOD method,
                CELL_LOC outloc, Solver* solver) {

  ASSERT1_FIELDS_COMPATIBLE(f, g);
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation())

  Mesh* mesh = f.getMesh();

  Field3D result(mesh);

  switch (method) {
  case BRACKET_CTU:
    throw BoutException("Bracket method CTU is not yet implemented for [2d,3d] fields.");
    break;
  case BRACKET_ARAKAWA:
    // It is symmetric, therefore we can return -[3d,2d]
    return -bracket(g, f, method, outloc, solver);
    break;
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDZ(Field3D{-DDX(f, outloc)}, g, outloc);
    break;
  }
  default: {
    // Use full expression with all terms
    const Coordinates* metric = f.getCoordinates(outloc);
    result = b0xGrad_dot_Grad(f, g, outloc) / metric->Bxy();
  }
  }

  return result;
}

Field3D bracket(const Field3D& f, const Field3D& g, BRACKET_METHOD method,
                CELL_LOC outloc, [[maybe_unused]] Solver* solver) {
  ASSERT1_FIELDS_COMPATIBLE(f, g);
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation());

  Mesh* mesh = f.getMesh();

  Field3D result{emptyFrom(f).setLocation(outloc)};

  const Coordinates* metric = f.getCoordinates(outloc);

  if (mesh->GlobalNx == 1 || mesh->GlobalNz == 1) {
    result = 0;
    result.setLocation(outloc);
    return result;
  }

  switch (method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)
#if not(BOUT_USE_METRIC_3D)
    if (solver == nullptr) {
      throw BoutException("CTU method requires access to the solver");
    }

    // Get current timestep
    const BoutReal dt = solver->getCurrentTimestep();

    FieldPerp vx(mesh);
    FieldPerp vz(mesh);
    vx.allocate();
    vx.setLocation(outloc);
    vz.allocate();
    vz.setLocation(outloc);

    const int ncz = mesh->LocalNz;
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      for (int x = 1; x <= mesh->LocalNx - 2; x++) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          const int zm = (z - 1 + ncz) % ncz;
          const int zp = (z + 1) % ncz;

          // Vx = DDZ(f)
          vx(x, z) = (f(x, y, zp) - f(x, y, zm)) / (2. * metric->dz(x, y, z));
          // Vz = -DDX(f)
          vz(x, z) = (f(x - 1, y, z) - f(x + 1, y, z))
                     / ((0.5 * metric->dx(x - 1, y)) + metric->dx(x, y)
                        + (0.5 * metric->dx(x + 1, y)));

          // Set stability condition
          solver->setMaxTimestep(fabs(metric->dx(x, y)) / (fabs(vx(x, z)) + 1e-16));
          solver->setMaxTimestep(metric->dz(x, y) / (fabs(vz(x, z)) + 1e-16));
        }
      }

      // Simplest form: use cell-centered velocities (no divergence included so not flux conservative)

      for (int x = mesh->xstart; x <= mesh->xend; x++) {
        for (int z = 0; z < ncz; z++) {
          const int zm = (z - 1 + ncz) % ncz;
          const int zp = (z + 1) % ncz;

          BoutReal gp;
          BoutReal gm;

          // X differencing
          if (vx(x, z) > 0.0) {
            gp = g(x, y, z)
                 + ((0.5 * dt / metric->dz(x, y))
                    * ((vz(x, z) > 0) ? vz(x, z) * (g(x, y, zm) - g(x, y, z))
                                      : vz(x, z) * (g(x, y, z) - g(x, y, zp))));

            gm = g(x - 1, y, z)
                 + ((0.5 * dt / metric->dz(x, y))
                    * ((vz(x, z) > 0) ? vz(x, z) * (g(x - 1, y, zm) - g(x - 1, y, z))
                                      : vz(x, z) * (g(x - 1, y, z) - g(x - 1, y, zp))));

          } else {
            gp = g(x + 1, y, z)
                 + ((0.5 * dt / metric->dz(x, y))
                    * ((vz(x, z) > 0) ? vz(x, z) * (g(x + 1, y, zm) - g(x + 1, y, z))
                                      : vz[x][z] * (g(x + 1, y, z) - g(x + 1, y, zp))));

            gm = g(x, y, z)
                 + ((0.5 * dt / metric->dz(x, y))
                    * ((vz(x, z) > 0) ? vz(x, z) * (g(x, y, zm) - g(x, y, z))
                                      : vz(x, z) * (g(x, y, z) - g(x, y, zp))));
          }

          result(x, y, z) = vx(x, z) * (gp - gm) / metric->dx(x, y);

          // Z differencing
          if (vz(x, z) > 0.0) {
            gp = g(x, y, z)
                 + ((0.5 * dt / metric->dx(x, y))
                    * ((vx[x][z] > 0) ? vx[x][z] * (g(x - 1, y, z) - g(x, y, z))
                                      : vx[x][z] * (g(x, y, z) - g(x + 1, y, z))));

            gm = g(x, y, zm)
                 + ((0.5 * dt / metric->dx(x, y))
                    * ((vx(x, z) > 0) ? vx(x, z) * (g(x - 1, y, zm) - g(x, y, zm))
                                      : vx(x, z) * (g(x, y, zm) - g(x + 1, y, zm))));
          } else {
            gp = g(x, y, zp)
                 + ((0.5 * dt / metric->dx(x, y))
                    * ((vx(x, z) > 0) ? vx(x, z) * (g(x - 1, y, zp) - g(x, y, zp))
                                      : vx(x, z) * (g(x, y, zp) - g(x + 1, y, zp))));

            gm = g(x, y, z)
                 + ((0.5 * dt / metric->dx(x, y))
                    * ((vx(x, z) > 0) ? vx(x, z) * (g(x - 1, y, z) - g(x, y, z))
                                      : vx(x, z) * (g(x, y, z) - g(x + 1, y, z))));
          }

          result(x, y, z) += vz(x, z) * (gp - gm) / metric->dz(x, y);
        }
      }
    }
#else
    throw BoutException("BRACKET_CTU not valid with 3D metrics yet.");
#endif
    break;
  }
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow
    const int ncz = mesh->LocalNz;

    // We need to discard const qualifier in order to manipulate
    // storage array directly
    Field3D f_temp = f;
    Field3D g_temp = g;

    BOUT_FOR(j2D, result.getRegion2D("RGN_NOBNDRY")) {
#if not(BOUT_USE_METRIC_3D)
      const BoutReal spacingFactor = 1.0 / (12 * metric->dz()[j2D] * metric->dx()[j2D]);
#endif
      const int jy = j2D.y();
      const int jx = j2D.x();
      const int xm = jx - 1;
      const int xp = jx + 1;

      const auto Fxm = f_temp(xm, jy);
      const auto Fx = f_temp(jx, jy);
      const auto Fxp = f_temp(xp, jy);
      const auto Gxm = g_temp(xm, jy);
      const auto Gx = g_temp(jx, jy);
      const auto Gxp = g_temp(xp, jy);

      // Here we split the loop over z into three parts; the first value, the middle block
      // and the last value
      // this is to allow the loop used in the middle block to vectorise.

      {
        const int jz = 0;
        const int jzp = 1;
        const int jzm = ncz - 1;
#if BOUT_USE_METRIC_3D
        const BoutReal spacingFactor =
            1.0 / (12 * metric->dz(jx, jy, jz) * metric->dx(jx, jy, jz));
#endif

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = (((Fx[jzp] - Fx[jzm]) * (Gxp[jz] - Gxm[jz]))
                              - ((Fxp[jz] - Fxm[jz]) * (Gx[jzp] - Gx[jzm])));

        // J+x
        const BoutReal Jpx =
            ((Gxp[jz] * (Fxp[jzp] - Fxp[jzm])) - (Gxm[jz] * (Fxm[jzp] - Fxm[jzm]))
             - (Gx[jzp] * (Fxp[jzp] - Fxm[jzp])) + (Gx[jzm] * (Fxp[jzm] - Fxm[jzm])));

        // Jx+
        const BoutReal Jxp =
            ((Gxp[jzp] * (Fx[jzp] - Fxp[jz])) - (Gxm[jzm] * (Fxm[jz] - Fx[jzm]))
             - (Gxm[jzp] * (Fx[jzp] - Fxm[jz])) + (Gxp[jzm] * (Fxp[jz] - Fx[jzm])));

        result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
      }

      for (int jz = 1; jz < mesh->LocalNz - 1; jz++) {
#if BOUT_USE_METRIC_3D
        const BoutReal spacingFactor =
            1.0 / (12 * metric->dz(jx, jy, jz) * metric->dx(jx, jy, jz));
#endif
        const int jzp = jz + 1;
        const int jzm = jz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = (((Fx[jzp] - Fx[jzm]) * (Gxp[jz] - Gxm[jz]))
                              - ((Fxp[jz] - Fxm[jz]) * (Gx[jzp] - Gx[jzm])));

        // J+x
        const BoutReal Jpx =
            ((Gxp[jz] * (Fxp[jzp] - Fxp[jzm])) - (Gxm[jz] * (Fxm[jzp] - Fxm[jzm]))
             - (Gx[jzp] * (Fxp[jzp] - Fxm[jzp])) + (Gx[jzm] * (Fxp[jzm] - Fxm[jzm])));

        // Jx+
        const BoutReal Jxp =
            ((Gxp[jzp] * (Fx[jzp] - Fxp[jz])) - (Gxm[jzm] * (Fxm[jz] - Fx[jzm]))
             - (Gxm[jzp] * (Fx[jzp] - Fxm[jz])) + (Gxp[jzm] * (Fxp[jz] - Fx[jzm])));

        result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
      }

      {
        const int jz = ncz - 1;
        const int jzp = 0;
        const int jzm = ncz - 2;
#if BOUT_USE_METRIC_3D
        const BoutReal spacingFactor =
            1.0 / (12 * metric->dz(jx, jy, jz) * metric->dx(jx, jy, jz));
#endif

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = (((Fx[jzp] - Fx[jzm]) * (Gxp[jz] - Gxm[jz]))
                              - ((Fxp[jz] - Fxm[jz]) * (Gx[jzp] - Gx[jzm])));

        // J+x
        const BoutReal Jpx =
            ((Gxp[jz] * (Fxp[jzp] - Fxp[jzm])) - (Gxm[jz] * (Fxm[jzp] - Fxm[jzm]))
             - (Gx[jzp] * (Fxp[jzp] - Fxm[jzp])) + (Gx[jzm] * (Fxp[jzm] - Fxm[jzm])));

        // Jx+
        const BoutReal Jxp =
            ((Gxp[jzp] * (Fx[jzp] - Fxp[jz])) - (Gxm[jzm] * (Fxm[jz] - Fx[jzm]))
             - (Gxm[jzp] * (Fx[jzp] - Fxm[jz])) + (Gxp[jzm] * (Fxp[jz] - Fx[jzm])));

        result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
      }
    }
    break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f, outloc), g, outloc) + VDDZ(-DDX(f, outloc), g, outloc);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g, outloc) / metric->Bxy();
  }
  }

  return result;
}
