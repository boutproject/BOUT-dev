/**************************************************************************
* Various differential operators defined on BOUT grid
*
**************************************************************************
* Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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

#include <globals.hxx>
#include <bout/solver.hxx>
#include <difops.hxx>
#include <vecops.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <fft.hxx>
#include <msg_stack.hxx>
#include <bout/assert.hxx>

#include <invert_laplace.hxx> // Delp2 uses same coefficients as inversion code

#include <interpolation.hxx>
#include <unused.hxx>

#include <cmath>

/*******************************************************************************
* Grad_par
* The parallel derivative along unperturbed B-field
*******************************************************************************/

const Field2D Grad_par(const Field2D &var, CELL_LOC outloc, const std::string &method) {
  return var.getCoordinates(outloc)->Grad_par(var, outloc, method);
}

const Field2D Grad_par(const Field2D &var, const std::string &method, CELL_LOC outloc) {
  return var.getCoordinates(outloc)->Grad_par(var, outloc, method);
}

const Field3D Grad_par(const Field3D &var, CELL_LOC outloc, const std::string &method) {
  return var.getCoordinates(outloc)->Grad_par(var, outloc, method);
}

const Field3D Grad_par(const Field3D &var, const std::string &method, CELL_LOC outloc) {
  return var.getCoordinates(outloc)->Grad_par(var, outloc, method);
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

const Field3D Grad_parP(const Field3D &apar, const Field3D &f) {
  ASSERT1(areFieldsCompatible(apar, f));
  ASSERT1(f.hasParallelSlices());

  Mesh *mesh = apar.getMesh();

  Field3D result{emptyFrom(f)};

  Coordinates *metric = apar.getCoordinates();

  Field3D gys{emptyFrom(f)};

  // Need Y derivative everywhere
  BOUT_FOR(i, f.getRegion("RGN_NOY")) {
    int x = i.x(), y = i.y(), z = i.z();
    gys[i] =
        (f.yup()(x, y + 1, z) - f.ydown()(x, y - 1, z))
        / (0.5 * metric->dy(x, y + 1) + metric->dy(x, y) + 0.5 * metric->dy(x, y - 1));
  }

  BOUT_FOR(i, f.getRegion("RGN_NOBNDRY")) {
    BoutReal by = 1. / sqrt(metric->g_22[i]);

    // Z indices zm and zp
    auto ixm = i.xm(), ixp = i.xp();
    auto iym = i.ym(), iyp = i.yp();
    auto izm = i.zm(), izp = i.zp();

    auto ixmym = ixm.ym(), ixmyp = ixm.yp();
    auto ixpym = ixp.ym(), ixpyp = ixp.yp();

    auto ixpzm = ixp.zm(), ixpzp = ixp.zp();
    auto ixmzm = ixm.zm(), ixmzp = ixm.zp();

    auto iypzm = iyp.zm(), iypzp = iyp.zp();
    auto iymzm = iym.zm(), iymzp = iym.zp();

    // bx = -DDZ(apar)
    BoutReal bx = (apar[izm] - apar[izp]) / (2. * metric->dz);
    // bz = DDX(f)
    BoutReal bz = (apar[ixp] - apar[ixm])
                  / (0.5 * metric->dx[ixm] + metric->dx[i] + 0.5 * metric->dx[ixp]);

    // Now calculate (bx*d/dx + by*d/dy + bz*d/dz) f

    // Length dl for predictor
    BoutReal dl = fabs(metric->dx[i]) / (fabs(bx) + 1e-16);
    dl = BOUTMIN(dl, fabs(metric->dy[i]) / (fabs(by) + 1e-16));
    dl = BOUTMIN(dl, metric->dz / (fabs(bz) + 1e-16));

    BoutReal fp, fm;

    // X differencing
    fp = f[ixp] + (0.25 * dl / metric->dz) * bz * (f[ixpzm] - f[ixpzp])
         - 0.5 * dl * by * gys[ixp];

    fm = f[ixm] + (0.25 * dl / metric->dz) * bz * (f[ixmzm] - f[ixmzp])
         - 0.5 * dl * by * gys[ixm];

    result[i] =
        bx * (fp - fm) / (0.5 * metric->dx[ixm] + metric->dx[i] + 0.5 * metric->dx[ixp]);

    // Z differencing

    fp = f[izp] + (0.25 * dl / metric->dx[i]) * bx * (f[ixmzp] - f[ixpzp])
         - 0.5 * dl * by * gys[izp];

    fm = f[izm] + (0.25 * dl / metric->dx[i]) * bx * (f[ixmzm] - f[ixpzm])
         - 0.5 * dl * by * gys[izm];

    result[i] += bz * (fp - fm) / (2. * metric->dz);

    // Y differencing

    fp = f.yup()[iyp]
         - 0.5 * dl * bx * (f.yup()[ixpyp] - f.yup()[ixmyp])
               / (0.5 * metric->dx[ixm] + metric->dx[i] + 0.5 * metric->dx[ixp])

         + (0.25 * dl / metric->dz) * bz * (f.yup()[iypzm] - f.yup()[iypzp]);

    fm = f.ydown()[iym]
         - 0.5 * dl * bx * (f.ydown()[ixpym] - f.ydown()[ixmym])
               / (0.5 * metric->dx[ixm] + metric->dx[i] + 0.5 * metric->dx[ixp])
         + (0.25 * dl / metric->dz) * bz * (f.ydown()[iymzm] - f.ydown()[iymzp]);

    result[i] +=
        by * (fp - fm) / (0.5 * metric->dy[iym] + metric->dy[i] + 0.5 * metric->dy[iyp]);
  }
  
  ASSERT2(result.getLocation() == f.getLocation());

  return result;
}

/*******************************************************************************
* Vpar_Grad_par
* vparallel times the parallel derivative along unperturbed B-field
*******************************************************************************/

const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f, CELL_LOC outloc, const std::string &method) {
  return f.getCoordinates(outloc)->Vpar_Grad_par(v, f, outloc, method);
}

const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f, const std::string &method, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Vpar_Grad_par(v, f, outloc, method);
}

const Field3D Vpar_Grad_par(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method) {
  return f.getCoordinates(outloc)->Vpar_Grad_par(v, f, outloc, method);
}

const Field3D Vpar_Grad_par(const Field3D &v, const Field3D &f, const std::string &method, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Vpar_Grad_par(v, f, outloc, method);
}

/*******************************************************************************
* Div_par
* parallel divergence operator B \partial_{||} (F/B)
*******************************************************************************/
const Field2D Div_par(const Field2D &f, CELL_LOC outloc, const std::string &method) {
  return f.getCoordinates(outloc)->Div_par(f, outloc, method);
}

const Field2D Div_par(const Field2D &f, const std::string &method, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Div_par(f, outloc, method);
}

const Field3D Div_par(const Field3D &f, CELL_LOC outloc, const std::string &method) {
  return f.getCoordinates(outloc)->Div_par(f, outloc, method);
}

const Field3D Div_par(const Field3D &f, const std::string &method, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Div_par(f, outloc, method);
}

const Field3D Div_par(const Field3D& f, const Field3D& v) {
  ASSERT1(areFieldsCompatible(f, v));
  ASSERT1(f.hasParallelSlices());
  ASSERT1(v.hasParallelSlices());

  // Parallel divergence, using velocities at cell boundaries
  // Note: Not guaranteed to be flux conservative
  Mesh* mesh = f.getMesh();

  Field3D result{emptyFrom(f)};

  Coordinates *coord = f.getCoordinates();

  BOUT_FOR(i, f.getRegion("RGN_NOBNDRY")) {
    const auto iym = i.ym(), iyp = i.yp();

    // Value of f and v at left cell face
    BoutReal fL = 0.5 * (f[i] + f.ydown()[iym]);
    BoutReal vL = 0.5 * (v[i] + v.ydown()[iym]);

    BoutReal fR = 0.5 * (f[i] + f.yup()[iyp]);
    BoutReal vR = 0.5 * (v[i] + v.yup()[iyp]);

    // Calculate flux at right boundary (y+1/2)
    BoutReal fluxRight = fR * vR * (coord->J[i] + coord->J[iyp])
                         / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[iyp]));

    // Calculate at left boundary (y-1/2)
    BoutReal fluxLeft = fL * vL * (coord->J[i] + coord->J[iym])
                        / (sqrt(coord->g_22[i]) + sqrt(coord->g_22[iym]));

    result[i] = (fluxRight - fluxLeft) / (coord->dy[i] * coord->J[i]);
  }

  return result;
}

//////// Flux methods

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method) {
  Coordinates *metric = f.getCoordinates(outloc);

  Field2D Bxy_floc = f.getCoordinates()->Bxy;

  if (!f.hasParallelSlices()) {
    return metric->Bxy*FDDY(v, f/Bxy_floc, outloc, method)/sqrt(metric->g_22);
  }

  // Need to modify yup and ydown fields
  Field3D f_B = f / Bxy_floc;
  // Distinct fields
  f_B.splitParallelSlices();
  f_B.yup() = f.yup() / Bxy_floc;
  f_B.ydown() = f.ydown() / Bxy_floc;
  return metric->Bxy*FDDY(v, f_B, outloc, method)/sqrt(metric->g_22);
}

const Field3D Div_par_flux(const Field3D &v, const Field3D &f, const std::string &method, CELL_LOC outloc) {
  return Div_par_flux(v,f, outloc, method);
}

/*******************************************************************************
* Parallel derivatives converting between left and cell centred
* NOTE: These are a quick hack to test if this works. The whole staggered grid
*       thing needs to be thought through.
*******************************************************************************/

const Field3D Grad_par_CtoL(const Field3D &var) {
  ASSERT1(var.getLocation() == CELL_CENTRE);

  Mesh *mesh = var.getMesh();
  Field3D result{emptyFrom(var).setLocation(CELL_YLOW)};

  Coordinates *metric = var.getCoordinates(CELL_YLOW);

  if (var.hasParallelSlices()) {
    // NOTE: Need to calculate one more point than centred vars
    BOUT_FOR(i, var.getRegion("RGN_NOY")) {
      const auto iym = i.ym();
      result[i] = 2. * (var[i] - var.ydown()[iym])
                  / (metric->dy[i] * sqrt(metric->g_22[i])
                     + metric->dy[iym] * sqrt(metric->g_22[iym]));
    }
  } else {
    // No yup/ydown fields, so transform to cell centred
    Field3D var_fa = toFieldAligned(var, RGN_NOX);

    BOUT_FOR(i, var.getRegion("RGN_NOY")) {
      const auto iym = i.ym();
      result[i] =
          2. * (var_fa[i] - var_fa[iym]) / (metric->dy[i] * sqrt(metric->g_22[i])
                                            + metric->dy[iym] * sqrt(metric->g_22[iym]));
    }

    result = fromFieldAligned(result, RGN_NOBNDRY);
  }

  return result;
}

const Field2D Grad_par_CtoL(const Field2D &var) {
  Field2D result{emptyFrom(var).setLocation(CELL_YLOW)};

  Coordinates *metric = var.getCoordinates(CELL_YLOW);

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    result[i] = (var[i] - var[i.ym()]) / (metric->dy[i] * sqrt(metric->g_22[i]));
  }

  return result;
}


const Field3D Vpar_Grad_par_LCtoC(const Field3D &v, const Field3D &f, REGION region) {
  ASSERT1(v.getMesh() == f.getMesh());
  ASSERT1(v.getLocation() == CELL_YLOW);
  ASSERT1(f.getLocation() == CELL_CENTRE);

  Field3D result{emptyFrom(f).setLocation(CELL_CENTRE)};

  bool vUseParallelSlices = v.hasParallelSlices();
  bool fUseParallelSlices = f.hasParallelSlices();

  if (vUseParallelSlices && fUseParallelSlices) {
    // Both v and f have up/down fields
    BOUT_OMP(parallel) {
      stencil fval, vval;
      BOUT_FOR_INNER(i, result.getRegion(region)) {
        vval.m = v.ydown()[i.ym()];
        vval.c = v[i];
        vval.p = v.yup()[i.yp()];

        fval.m = f.ydown()[i.ym()];
        fval.c = f[i];
        fval.p = f.yup()[i.yp()];

        // Left side
        result[i] = (vval.c >= 0.0) ? vval.c * fval.m : vval.c * fval.c;
        // Right side
        result[i] -= (vval.p >= 0.0) ? vval.p * fval.c : vval.p * fval.p;
      }
    }
  }
  else {
    // Both must shift to field aligned
    // (even if one of v and f has yup/ydown fields, it doesn't make sense to
    // multiply them with one in field-aligned and one in non-field-aligned
    // coordinates)
    Field3D v_fa = toFieldAligned(v, RGN_NOX);
    Field3D f_fa = toFieldAligned(f, RGN_NOX);

    BOUT_OMP(parallel) {
      stencil fval, vval;
      BOUT_FOR_INNER(i, result.getRegion(region)) {
        fval.m = f_fa[i.ym()];
        fval.c = f_fa[i];
        fval.p = f_fa[i.yp()];

        vval.m = v_fa[i.ym()];
        vval.c = v_fa[i];
        vval.p = v_fa[i.yp()];

        // Left side
        result[i] = (vval.c >= 0.0) ? vval.c * fval.m : vval.c * fval.c;
        // Right side
        result[i] -= (vval.p >= 0.0) ? vval.p * fval.c : vval.p * fval.p;
      }

      result = fromFieldAligned(result, region);
    }
  }

  return result;
}

const Field3D Grad_par_LtoC(const Field3D &var) {
  const auto varMesh = var.getMesh();

  if (varMesh->StaggerGrids) {
    ASSERT1(var.getLocation() == CELL_YLOW);
  }

  Field3D result{emptyFrom(var).setLocation(CELL_CENTRE)};

  Coordinates *metric = var.getCoordinates(CELL_CENTRE);

  if (var.hasParallelSlices()) {
    BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
      result[i] = (var.yup()[i.yp()] - var[i]) / (metric->dy[i]*sqrt(metric->g_22[i]));
    }
  } else {
    // No yup/ydown field, so transform to field aligned

    Field3D var_fa = toFieldAligned(var, RGN_NOX);

    BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
      result[i] = (var_fa[i.yp()] - var_fa[i]) / (metric->dy[i]*sqrt(metric->g_22[i]));
    }
    result = fromFieldAligned(result, RGN_NOBNDRY);
  }

  return result;
}

const Field2D Grad_par_LtoC(const Field2D &var) {
  Field2D result{emptyFrom(var).setLocation(CELL_CENTRE)};

  Coordinates *metric = var.getCoordinates(CELL_CENTRE);

  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    result[i] = (var[i.yp()] - var[i]) / (metric->dy[i] * sqrt(metric->g_22[i]));
  }

  return result;
}

const Field2D Div_par_LtoC(const Field2D &var) {
  return var.getCoordinates(CELL_CENTRE)->Bxy *
         Grad_par_LtoC(var / var.getCoordinates(CELL_YLOW)->Bxy);
}

const Field3D Div_par_LtoC(const Field3D &var) {
  Mesh* mesh = var.getMesh();

  Field3D result{emptyFrom(var).setLocation(CELL_CENTRE)};

  Coordinates *metric = var.getCoordinates(CELL_CENTRE);

  // NOTE: Need to calculate one more point than centred vars
  for (int jx = 0; jx < mesh->LocalNx; jx++) {
    for (int jy = 0; jy < mesh->LocalNy - 1; jy++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        result(jx, jy, jz) = metric->Bxy(jx, jy) * 2. *
                             (var.yup()(jx, jy + 1, jz) / metric->Bxy(jx, jy + 1) -
                              var(jx, jy, jz) / metric->Bxy(jx, jy)) /
                             (metric->dy(jx, jy) * sqrt(metric->g_22(jx, jy)) +
                              metric->dy(jx, jy - 1) * sqrt(metric->g_22(jx, jy - 1)));
      }
    }
  }

  return result;
}

const Field2D Div_par_CtoL(const Field2D &var) {
  return var.getCoordinates(CELL_CENTRE)->Bxy *
         Grad_par_CtoL(var / var.getCoordinates(CELL_YLOW)->Bxy);
}

const Field3D Div_par_CtoL(const Field3D &var) {
  Mesh* mesh = var.getMesh();

  Field3D result{emptyFrom(var).setLocation(CELL_CENTRE)};

  Coordinates *metric = var.getCoordinates(CELL_CENTRE);

  // NOTE: Need to calculate one more point than centred vars
  for (int jx = 0; jx < mesh->LocalNx; jx++) {
    for (int jy = 1; jy < mesh->LocalNy; jy++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        result(jx, jy, jz) = metric->Bxy(jx, jy) * 2. *
                             (var(jx, jy, jz) / metric->Bxy(jx, jy) -
                              var.ydown()(jx, jy - 1, jz) / metric->Bxy(jx, jy - 1)) /
                             (metric->dy(jx, jy) * sqrt(metric->g_22(jx, jy)) +
                              metric->dy(jx, jy - 1) * sqrt(metric->g_22(jx, jy - 1)));
      }
    }
  }

  return result;
}

/*******************************************************************************
* Grad2_par2
* second parallel derivative
*
* (b dot Grad)(b dot Grad)
*
* Note: For parallel Laplacian use LaplacePar
*******************************************************************************/

const Field2D Grad2_par2(const Field2D &f, CELL_LOC outloc, const std::string &method) {
  return f.getCoordinates(outloc)->Grad2_par2(f, outloc, method);
}

const Field3D Grad2_par2(const Field3D &f, CELL_LOC outloc, const std::string &method) {
  return f.getCoordinates(outloc)->Grad2_par2(f, outloc, method);
}

/*******************************************************************************
* Div_par_K_Grad_par
* Parallel divergence of diffusive flux, K*Grad_par
*******************************************************************************/

const Field2D Div_par_K_Grad_par(BoutReal kY, const Field2D &f, CELL_LOC outloc) {
  return kY*Grad2_par2(f, outloc);
}

const Field3D Div_par_K_Grad_par(BoutReal kY, const Field3D &f, CELL_LOC outloc) {
  return kY*Grad2_par2(f, outloc);
}

const Field2D Div_par_K_Grad_par(const Field2D &kY, const Field2D &f, CELL_LOC outloc) {
  return interp_to(kY, outloc)*Grad2_par2(f, outloc) + Div_par(kY, outloc)*Grad_par(f, outloc);
}

const Field3D Div_par_K_Grad_par(const Field2D &kY, const Field3D &f, CELL_LOC outloc) {
  return interp_to(kY, outloc)*Grad2_par2(f, outloc) + Div_par(kY, outloc)*Grad_par(f, outloc);
}

const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field2D &f, CELL_LOC outloc) {
  return interp_to(kY, outloc)*Grad2_par2(f, outloc) + Div_par(kY, outloc)*Grad_par(f, outloc);
}

const Field3D Div_par_K_Grad_par(const Field3D &kY, const Field3D &f, CELL_LOC outloc) {
  return interp_to(kY, outloc)*Grad2_par2(f, outloc) + Div_par(kY, outloc)*Grad_par(f, outloc);
}

/*******************************************************************************
* Delp2
* perpendicular Laplacian operator
*******************************************************************************/

const Field2D Delp2(const Field2D& f, CELL_LOC outloc, bool useFFT) {
  return f.getCoordinates(outloc)->Delp2(f, outloc, useFFT);
}

const Field3D Delp2(const Field3D& f, CELL_LOC outloc, bool useFFT) {
  return f.getCoordinates(outloc)->Delp2(f, outloc, useFFT);
}

const FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc, bool useFFT) {
  return f.getCoordinates(outloc)->Delp2(f, outloc, useFFT);
}

/*******************************************************************************
* LaplacePerp
* Full perpendicular Laplacian operator on scalar field
*
* Laplace_perp = Laplace - Laplace_par
*******************************************************************************/

const Field2D Laplace_perp(const Field2D &f, CELL_LOC outloc) {
  return Laplace(f, outloc) - Laplace_par(f, outloc);
}

const Field3D Laplace_perp(const Field3D &f, CELL_LOC outloc) {
  return Laplace(f, outloc) - Laplace_par(f, outloc);
}

/*******************************************************************************
* LaplacePar
* Full parallel Laplacian operator on scalar field
*
* LaplacePar(f) = Div( b (b dot Grad(f)) ) 
*
*******************************************************************************/

const Field2D Laplace_par(const Field2D &f, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Laplace_par(f, outloc);
}

const Field3D Laplace_par(const Field3D &f, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Laplace_par(f, outloc);
}

/*******************************************************************************
* Laplacian
* Full Laplacian operator on scalar field
*******************************************************************************/

const Field2D Laplace(const Field2D &f, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Laplace(f, outloc);
}

const Field3D Laplace(const Field3D &f, CELL_LOC outloc) {
  return f.getCoordinates(outloc)->Laplace(f, outloc);
}

/*******************************************************************************
* b0xGrad_dot_Grad
* Terms of form b0 x Grad(phi) dot Grad(A)
* Used for ExB terms and perturbed B field using A_||
*******************************************************************************/

const Field2D b0xGrad_dot_Grad(const Field2D &phi, const Field2D &A, CELL_LOC outloc) {
  
  TRACE("b0xGrad_dot_Grad( Field2D , Field2D )");
  
  if (outloc == CELL_DEFAULT) outloc = A.getLocation();

  ASSERT1(phi.getMesh() == A.getMesh());

  Coordinates *metric = phi.getCoordinates(outloc);

  // Calculate phi derivatives
  Field2D dpdx = DDX(phi, outloc);
  Field2D dpdy = DDY(phi, outloc);
  
  // Calculate advection velocity
  Field2D vx = -metric->g_23*dpdy;
  Field2D vy = metric->g_23*dpdx;

  // Upwind A using these velocities
  Field2D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc);
  result /= metric->J*sqrt(metric->g_22);

  ASSERT1(result.getLocation() == outloc);

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field2D &phi, const Field3D &A, CELL_LOC outloc) {
  TRACE("b0xGrad_dot_Grad( Field2D , Field3D )");

  if (outloc == CELL_DEFAULT) outloc = A.getLocation();

  ASSERT1(phi.getMesh() == A.getMesh());

  Mesh *mesh = phi.getMesh();

  Coordinates *metric = phi.getCoordinates(outloc);

  // Calculate phi derivatives
  Field2D dpdx = DDX(phi, outloc);
  Field2D dpdy = DDY(phi, outloc);

  // Calculate advection velocity
  Field2D vx = -metric->g_23 * dpdy;
  Field2D vy = metric->g_23 * dpdx;
  Field2D vz = metric->g_12 * dpdy - metric->g_22 * dpdx;

  if(mesh->IncIntShear) {
    // BOUT-06 style differencing
    vz += metric->IntShiftTorsion * vx;
  }

  // Upwind A using these velocities

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc) + VDDZ(vz, A, outloc);

  result /= (metric->J*sqrt(metric->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif

  ASSERT2(result.getLocation() == outloc);
  
  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &p, const Field2D &A, CELL_LOC outloc) {
  TRACE("b0xGrad_dot_Grad( Field3D , Field2D )");

  if (outloc == CELL_DEFAULT) outloc = A.getLocation();

  ASSERT1(p.getMesh() == A.getMesh());

  Coordinates *metric = p.getCoordinates(outloc);

  // Calculate phi derivatives
  Field3D dpdx = DDX(p, outloc);
  Field3D dpdy = DDY(p, outloc);
  Field3D dpdz = DDZ(p, outloc);

  // Calculate advection velocity
  Field3D vx = metric->g_22 * dpdz - metric->g_23 * dpdy;
  Field3D vy = metric->g_23 * dpdx - metric->g_12 * dpdz;

  // Upwind A using these velocities

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc);

  result /= (metric->J*sqrt(metric->g_22));
  
#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+p.name+","+A.name+")";
#endif
  
  ASSERT2(result.getLocation() == outloc);

  return result;
}

const Field3D b0xGrad_dot_Grad(const Field3D &phi, const Field3D &A, CELL_LOC outloc) {
  TRACE("b0xGrad_dot_Grad( Field3D , Field3D )");

  if (outloc == CELL_DEFAULT) outloc = A.getLocation();

  ASSERT1(phi.getMesh() == A.getMesh());

  Mesh * mesh = phi.getMesh();

  Coordinates *metric = phi.getCoordinates(outloc);

  // Calculate phi derivatives
  Field3D dpdx = DDX(phi, outloc);
  Field3D dpdy = DDY(phi, outloc);
  Field3D dpdz = DDZ(phi, outloc);

  // Calculate advection velocity
  Field3D vx = metric->g_22 * dpdz - metric->g_23 * dpdy;
  Field3D vy = metric->g_23 * dpdx - metric->g_12 * dpdz;
  Field3D vz = metric->g_12 * dpdy - metric->g_22 * dpdx;

  if(mesh->IncIntShear) {
    // BOUT-06 style differencing
    vz += metric->IntShiftTorsion * vx;
  }

  Field3D result = VDDX(vx, A, outloc) + VDDY(vy, A, outloc) + VDDZ(vz, A, outloc);

  result /=  (metric->J*sqrt(metric->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif

  ASSERT2(result.getLocation() == outloc);

  return result;
}

/*******************************************************************************
 * Poisson bracket
 * Terms of form b0 x Grad(f) dot Grad(g) / B = [f, g]
 *******************************************************************************/

const Field2D bracket(const Field2D &f, const Field2D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *UNUSED(solver)) {
  TRACE("bracket(Field2D, Field2D)");

  ASSERT1(areFieldsCompatible(f, g));
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation());

  Field2D result{emptyFrom(f)};

  if( (method == BRACKET_SIMPLE) || (method == BRACKET_ARAKAWA)) {
    // Use a subset of terms for comparison to BOUT-06
    result = 0.0;
    result.setLocation(outloc);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g, outloc) / f.getCoordinates(outloc)->Bxy;
  }
  return result;
}

const Field3D bracket(const Field3D &f, const Field2D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *solver) {
  TRACE("bracket(Field3D, Field2D)");

  ASSERT1(areFieldsCompatible(f, g));
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation());

  Mesh *mesh = f.getMesh();

  Field3D result{emptyFrom(f).setLocation(outloc)};

  Coordinates *metric = f.getCoordinates(outloc);

  switch(method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)
    
    if(!solver)
      throw BoutException("CTU method requires access to the solver");
    
    BOUT_FOR(i, f.getRegion("RGN_NOBNDRY")) {
      const auto izm = i.zm(), izp = i.zp();

      BoutReal gp, gm;

      // Vx = DDZ(f)
      BoutReal vx = (f[izp] - f[izm]) / (2. * metric->dz);

      // Set stability condition
      solver->setMaxTimestep(metric->dx[i] / (fabs(vx) + 1e-16));

      // X differencing
      if (vx > 0.0) {
        gp = g[i];

        gm = g[i.xm()];

      } else {
        gp = g[i.xp()];

        gm = g[i];
      }

      result[i] = vx * (gp - gm) / metric->dx[i];
    }

    break;
  }
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow. Here as a test

    const BoutReal fac = 1.0 / (12 * metric->dz);

#ifdef BOUT_HAS_Z_GUARD_CELLS_IMPLEMENTED
    // Not making change yet as without z-guard cells this new version
    // won't vectorise.
    BOUT_FOR(i, result.getRegion2D("RGN_NOBNDRY")) {
      // Get constants for this iteration
      const BoutReal spacingFactor = fac / metric->dx[i];

      const auto xm = i.xm(), xp = i.xp();

      // Extract relevant Field2D values
      const BoutReal gxm = gp[xm], gc = g[i], gxp = g[xp];

      const auto jStart = mesh->ind2Dto3D(i, mesh->zstart);
      const auto jEnd = mesh->ind2Dto3D(i, mesh->zend);
      for (auto j = jStart; j <= jEnd; j++) {
        const auto zp = j.zp(), xmzp = zp.xm(), xpzp = zp.xp();
        const auto zm = j.zm(), xmzm = zm.xm(), xpzm = zm.xp();

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = 2 * (f[zp] - f[zm]) * (gxp - gxm);

        // J+x
        const BoutReal Jpx = gxp * (f[xpzp] - f[xpzm]) - gxm * (f[xmzp] - f[xmzm])
                             + gc * (f[xpzm] - f[xpzp] - f[xmzm] + f[xmzp]);

        result[j] = (Jpp + Jpx) * spacingFactor;
      }
    }
#else
    // Following loop tagged for replacement
    const int ncz = mesh->LocalNz;

    BOUT_FOR(j2D, result.getRegion2D("RGN_NOBNDRY")) {
      // Get constants for this iteration
      const BoutReal spacingFactor = fac / metric->dx[j2D];
      const int jy = j2D.y(), jx = j2D.x();
      const int xm = jx - 1, xp = jx + 1;

      // Extract relevant Field2D values
      const BoutReal gxm = g(xm, jy), gc = g(jx, jy), gxp = g(xp, jy);

      // Index Field3D as 2D to get start of z data block
      const auto fxm = f(xm, jy), fc = f(jx, jy), fxp = f(xp, jy);

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
        const BoutReal Jpx = gxp * (fxp[jzp] - fxp[jzm]) - gxm * (fxm[jzp] - fxm[jzm]) +
                             gc * (fxp[jzm] - fxp[jzp] - fxm[jzm] + fxm[jzp]);

        result(jx, jy, 0) = (Jpp + Jpx) * spacingFactor;
      }

      // The middle block
      for (int jz = 1; jz < mesh->LocalNz - 1; jz++) {
        const int jzp = jz + 1;
        const int jzm = jz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = 2 * (fc[jzp] - fc[jzm]) * (gxp - gxm);

        // J+x
        const BoutReal Jpx = gxp * (fxp[jzp] - fxp[jzm]) - gxm * (fxm[jzp] - fxm[jzm]) +
                             gc * (fxp[jzm] - fxp[jzp] - fxm[jzm] + fxm[jzp]);

        result(jx, jy, jz) = (Jpp + Jpx) * spacingFactor;
      }

      // The last value
      {
        const int jzp = 0;
        const int jzm = ncz - 2;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = 2 * (fc[jzp] - fc[jzm]) * (gxp - gxm);

        // J+x
        const BoutReal Jpx = gxp * (fxp[jzp] - fxp[jzm]) - gxm * (fxm[jzp] - fxm[jzm]) +
                             gc * (fxp[jzm] - fxp[jzp] - fxm[jzm] + fxm[jzp]);

        result(jx, jy, ncz - 1) = (Jpp + Jpx) * spacingFactor;
      }
    }
#endif
    break;
  }
  case BRACKET_ARAKAWA_OLD: {
   const BoutReal partialFactor = 1.0/(12 * metric->dz);
    BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
      // Being lazy in migrating to a BOUT_FOR loop from a nested loop
      // just calculate individual indices for each direction -- don't have
      // to change the main loop body but this is probably inefficient.
      const int jx = i.x(), jy = i.y(), jz = i.z();
      const int jzp = i.zp().ind;
      const int jzm = i.zm().ind;

      const BoutReal spacingFactor = partialFactor / metric->dx(jx, jy);

      // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
      BoutReal Jpp =
          ((f(jx, jy, jzp) - f(jx, jy, jzm)) * (g(jx + 1, jy) - g(jx - 1, jy))
           - (f(jx + 1, jy, jz) - f(jx - 1, jy, jz)) * (g(jx, jy) - g(jx, jy)));

      // J+x
      BoutReal Jpx = (g(jx + 1, jy) * (f(jx + 1, jy, jzp) - f(jx + 1, jy, jzm))
                      - g(jx - 1, jy) * (f(jx - 1, jy, jzp) - f(jx - 1, jy, jzm))
                      - g(jx, jy) * (f(jx + 1, jy, jzp) - f(jx - 1, jy, jzp))
                      + g(jx, jy) * (f(jx + 1, jy, jzm) - f(jx - 1, jy, jzm)));

      // Jx+
      BoutReal Jxp = (g(jx + 1, jy) * (f(jx, jy, jzp) - f(jx + 1, jy, jz))
                      - g(jx - 1, jy) * (f(jx - 1, jy, jz) - f(jx, jy, jzm))
                      - g(jx - 1, jy) * (f(jx, jy, jzp) - f(jx - 1, jy, jz))
                      + g(jx + 1, jy) * (f(jx + 1, jy, jz) - f(jx, jy, jzm)));

      result[i] = (Jpp + Jpx + Jxp) * spacingFactor;
    }

    break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f, outloc), g, outloc);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g, outloc) / metric->Bxy;
  }
  }
  return result;
}

const Field3D bracket(const Field2D &f, const Field3D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *solver) {
  TRACE("bracket(Field2D, Field3D)");

  ASSERT1(areFieldsCompatible(f, g));
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation())

  Mesh *mesh = f.getMesh();

  Field3D result(mesh);

  switch(method) {
  case BRACKET_CTU:
    throw BoutException("Bracket method CTU is not yet implemented for [2d,3d] fields.");
    break;
  case BRACKET_ARAKAWA: 
    // It is symmetric, therefore we can return -[3d,2d]
    return -bracket(g,f,method,outloc,solver);
    break;
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDZ(-DDX(f, outloc), g, outloc);
    break;
  }
  default: {
    // Use full expression with all terms
    Coordinates *metric = f.getCoordinates(outloc);
    result = b0xGrad_dot_Grad(f, g, outloc) / metric->Bxy;
  }
  }
  
  return result;
}

const Field3D bracket(const Field3D &f, const Field3D &g, BRACKET_METHOD method,
                      CELL_LOC outloc, Solver *solver) {
  TRACE("Field3D, Field3D");

  ASSERT1(areFieldsCompatible(f, g));
  if (outloc == CELL_DEFAULT) {
    outloc = g.getLocation();
  }
  ASSERT1(outloc == g.getLocation());

  Mesh *mesh = f.getMesh();

  Field3D result{emptyFrom(f).setLocation(outloc)};

  Coordinates *metric = f.getCoordinates(outloc);

  if (mesh->GlobalNx == 1 || mesh->GlobalNz == 1) {
    result=0;
    result.setLocation(outloc);
    return result;
  }
  
  switch(method) {
  case BRACKET_CTU: {
    // First order Corner Transport Upwind method
    // P.Collela JCP 87, 171-200 (1990)
    
    if(!solver)
      throw BoutException("CTU method requires access to the solver");

    // Get current timestep
    BoutReal dt = solver->getCurrentTimestep();
    
    FieldPerp vx(mesh), vz(mesh);
    vx.allocate();
    vx.setLocation(outloc);
    vz.allocate();
    vz.setLocation(outloc);

#ifdef BOUT_HAS_Z_GUARD_CELLS_IMPLEMENTED
    throw BoutException("Bracket_CTU needs updating to account for z-guard cells.");
#endif

    int ncz = mesh->LocalNz;
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      for(int x=1;x<=mesh->LocalNx-2;x++) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          // Following two lines need updating for z with guard-cells
          int zm = (z - 1 + ncz) % ncz;
          int zp = (z + 1) % ncz;
          
          // Vx = DDZ(f)
          vx(x,z) = (f(x,y,zp) - f(x,y,zm))/(2.*metric->dz);
          // Vz = -DDX(f)
          vz(x,z) = (f(x-1,y,z) - f(x+1,y,z))/(0.5*metric->dx(x-1,y) + metric->dx(x,y) + 0.5*metric->dx(x+1,y));
          
          // Set stability condition
          solver->setMaxTimestep(fabs(metric->dx(x,y)) / (fabs(vx(x,z)) + 1e-16));
          solver->setMaxTimestep(metric->dz / (fabs(vz(x,z)) + 1e-16));
        }
      }
      
      // Simplest form: use cell-centered velocities (no divergence included so not flux conservative)
      
      for(int x=mesh->xstart;x<=mesh->xend;x++)
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          // Following two lines need updating for z with guard-cells
          int zm = (z - 1 + ncz) % ncz;
          int zp = (z + 1) % ncz;

          BoutReal gp, gm;

          // X differencing
          if (vx(x, z) > 0.0) {
            gp = g(x, y, z) +
                 (0.5 * dt / metric->dz) * ((vz(x, z) > 0)
                                                ? vz(x, z) * (g(x, y, zm) - g(x, y, z))
                                                : vz(x, z) * (g(x, y, z) - g(x, y, zp)));

            gm = g(x - 1, y, z) +
                 (0.5 * dt / metric->dz) *
                     ((vz(x, z) > 0) ? vz(x, z) * (g(x - 1, y, zm) - g(x - 1, y, z))
                                     : vz(x, z) * (g(x - 1, y, z) - g(x - 1, y, zp)));

          } else {
            gp = g(x + 1, y, z) +
                 (0.5 * dt / metric->dz) *
                     ((vz(x, z) > 0) ? vz(x, z) * (g(x + 1, y, zm) - g(x + 1, y, z))
                                     : vz[x][z] * (g(x + 1, y, z) - g(x + 1, y, zp)));

            gm = g(x, y, z) +
                 (0.5 * dt / metric->dz) * ((vz(x, z) > 0)
                                                ? vz(x, z) * (g(x, y, zm) - g(x, y, z))
                                                : vz(x, z) * (g(x, y, z) - g(x, y, zp)));
          }

          result(x, y, z) = vx(x, z) * (gp - gm) / metric->dx(x, y);

          // Z differencing
          if (vz(x, z) > 0.0) {
            gp = g(x, y, z) +
                 (0.5 * dt / metric->dx(x, y)) *
                     ((vx[x][z] > 0) ? vx[x][z] * (g(x - 1, y, z) - g(x, y, z))
                                     : vx[x][z] * (g(x, y, z) - g(x + 1, y, z)));

            gm = g(x, y, zm) +
                 (0.5 * dt / metric->dx(x, y)) *
                     ((vx(x, z) > 0) ? vx(x, z) * (g(x - 1, y, zm) - g(x, y, zm))
                                     : vx(x, z) * (g(x, y, zm) - g(x + 1, y, zm)));
          } else {
            gp = g(x, y, zp) +
                 (0.5 * dt / metric->dx(x, y)) *
                     ((vx(x, z) > 0) ? vx(x, z) * (g(x - 1, y, zp) - g(x, y, zp))
                                     : vx(x, z) * (g(x, y, zp) - g(x + 1, y, zp)));

            gm = g(x, y, z) +
                 (0.5 * dt / metric->dx(x, y)) *
                     ((vx(x, z) > 0) ? vx(x, z) * (g(x - 1, y, z) - g(x, y, z))
                                     : vx(x, z) * (g(x, y, z) - g(x + 1, y, z)));
          }

          result(x, y, z) += vz(x, z) * (gp - gm) / metric->dz;
        }
    }
    break;
  }
  case BRACKET_ARAKAWA: {
    // Arakawa scheme for perpendicular flow
    const BoutReal partialFactor = 1.0/(12 * metric->dz);

#ifdef BOUT_HAS_Z_GUARD_CELLS_IMPLEMENTED
    // Not making change yet as without z-guard cells this new version
    // won't vectorise.
    BOUT_FOR(i, result.getRegion2D("RGN_NOBNDRY")) {
      // Get constants for this iteration
      const BoutReal spacingFactor = fac / metric->dx[i];

      const auto xm = i.xm(), xp = i.xp();

      const auto jStart = mesh->ind2Dto3D(i, mesh->zstart);
      const auto jEnd = mesh->ind2Dto3D(i, mesh->zend);

      for (auto j = jStart; j <= jEnd; j++) {
        const auto zp = j.zp(), xmzp = zp.xm(), xpzp = zp.xp();
        const auto zm = j.zm(), xmzm = zm.xm(), xpzm = zm.xp();

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp =
            ((f[zp] - f[zm]) * (g[xp] - g[xm]) - (f[xp] - f[xm]) * (g[zp] - g[zm]));

        // J+x
        const BoutReal Jpx =
            (g[xp] * (f[xpzp] - f[xpzm]) - g[xm] * (f[xmzp] - f[xmzm])
             - g[zp] * (f[xpzp] - f[xmzp]) + g[zm] * (f[xpzm] - f[xmzm]));

        // Jx+
        const BoutReal Jxp = (g[xpzp] * (f[zp] - f[xp]) - g[xmzm] * (f[xm] - f[zm])
                              - g[xmzp] * (f[zp] - f[xm]) + g[xpzm] * (f[xp] - f[zm]));

        result[j] = (Jpp + Jpx + Jxp) * spacingFactor;
      }
#else
    // We need to discard const qualifier in order to manipulate
    // storage array directly
    Field3D f_temp = f;
    Field3D g_temp = g;

    const int ncz = mesh->LocalNz;

    BOUT_FOR(j2D, result.getRegion2D("RGN_NOBNDRY")) {
      const BoutReal spacingFactor = partialFactor / metric->dx[j2D];
      const int jy = j2D.y(), jx = j2D.x();
      const int xm = jx - 1, xp = jx + 1;

      const auto Fxm = f_temp(xm, jy), Fx = f_temp(jx, jy), Fxp = f_temp(xp, jy);
      const auto Gxm = g_temp(xm, jy), Gx = g_temp(jx, jy), Gxp = g_temp(xp, jy);

      // Here we split the loop over z into three parts; the first value, the middle block
      // and the last value
      // this is to allow the loop used in the middle block to vectorise.

      {
        const int jz = 0;
        const int jzp = 1;
        const int jzm = ncz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = ((Fx[jzp] - Fx[jzm]) * (Gxp[jz] - Gxm[jz]) -
                              (Fxp[jz] - Fxm[jz]) * (Gx[jzp] - Gx[jzm]));

        // J+x
        const BoutReal Jpx =
            (Gxp[jz] * (Fxp[jzp] - Fxp[jzm]) - Gxm[jz] * (Fxm[jzp] - Fxm[jzm]) -
             Gx[jzp] * (Fxp[jzp] - Fxm[jzp]) + Gx[jzm] * (Fxp[jzm] - Fxm[jzm]));

        // Jx+
        const BoutReal Jxp =
            (Gxp[jzp] * (Fx[jzp] - Fxp[jz]) - Gxm[jzm] * (Fxm[jz] - Fx[jzm]) -
             Gxm[jzp] * (Fx[jzp] - Fxm[jz]) + Gxp[jzm] * (Fxp[jz] - Fx[jzm]));

        result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
      }

      for (int jz = 1; jz < mesh->LocalNz - 1; jz++) {
        const int jzp = jz + 1;
        const int jzm = jz - 1;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = ((Fx[jzp] - Fx[jzm]) * (Gxp[jz] - Gxm[jz]) -
                              (Fxp[jz] - Fxm[jz]) * (Gx[jzp] - Gx[jzm]));

        // J+x
        const BoutReal Jpx =
            (Gxp[jz] * (Fxp[jzp] - Fxp[jzm]) - Gxm[jz] * (Fxm[jzp] - Fxm[jzm]) -
             Gx[jzp] * (Fxp[jzp] - Fxm[jzp]) + Gx[jzm] * (Fxp[jzm] - Fxm[jzm]));

        // Jx+
        const BoutReal Jxp =
            (Gxp[jzp] * (Fx[jzp] - Fxp[jz]) - Gxm[jzm] * (Fxm[jz] - Fx[jzm]) -
             Gxm[jzp] * (Fx[jzp] - Fxm[jz]) + Gxp[jzm] * (Fxp[jz] - Fx[jzm]));

        result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
      }

      {
        const int jz = ncz - 1;
        const int jzp = 0;
        const int jzm = ncz - 2;

        // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
        const BoutReal Jpp = ((Fx[jzp] - Fx[jzm]) * (Gxp[jz] - Gxm[jz]) -
                              (Fxp[jz] - Fxm[jz]) * (Gx[jzp] - Gx[jzm]));

        // J+x
        const BoutReal Jpx =
            (Gxp[jz] * (Fxp[jzp] - Fxp[jzm]) - Gxm[jz] * (Fxm[jzp] - Fxm[jzm]) -
             Gx[jzp] * (Fxp[jzp] - Fxm[jzp]) + Gx[jzm] * (Fxp[jzm] - Fxm[jzm]));

        // Jx+
        const BoutReal Jxp =
            (Gxp[jzp] * (Fx[jzp] - Fxp[jz]) - Gxm[jzm] * (Fxm[jz] - Fx[jzm]) -
             Gxm[jzp] * (Fx[jzp] - Fxm[jz]) + Gxp[jzm] * (Fxp[jz] - Fx[jzm]));

        result(jx, jy, jz) = (Jpp + Jpx + Jxp) * spacingFactor;
      }
    }
#endif
      break;
    }
  case BRACKET_ARAKAWA_OLD: {
    // Arakawa scheme for perpendicular flow
    // Just use main Arakawa method rather than maintaining two versions
    return bracket(f, g, BRACKET_ARAKAWA, outloc, solver);
    break;
  }
  case BRACKET_SIMPLE: {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(f, outloc), g, outloc) + VDDZ(-DDX(f, outloc), g, outloc);
    break;
  }
  default: {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(f, g, outloc) / metric->Bxy;
  }
  }
  
  return result;
}
