/**************************************************************************
 * Functions to interpolate between cell locations (e.g. lower Y and centred)
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
#include <interpolation.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <stencils.hxx>
#include <unused.hxx>
#include <bout/scorepwrapper.hxx>
#include <bout/indexoffset.hxx>

/// Perform interpolation between centre -> shifted or vice-versa
/*!
  Interpolate using 4th-order staggered formula

  @param[in] s  Input stencil. mm -> -3/2, m -> -1/2, p -> +1/2, pp -> +3/2
*/
BoutReal interp(const stencil &s) { return (9. * (s.m + s.p) - s.mm - s.pp) / 16.; }

/*!
  Interpolate between different cell locations

  NOTE: This requires communication if the result is required in guard cells
  NOTE: Since corner guard cells cannot be communicated, it never makes sense
  to calculate interpolation in guard cells. If guard cell values are required,
  we must communicate (unless interpolating in z). Since mesh->communicate()
  communicates both x- and y-guard cells by default, there is no difference
  between RGN_ALL, RGN_NOX and RGN_NOY.

  @param[in]   var  Input variable
  @param[in]   loc  Location of output values
  @param[in]   region  Region where output will be calculated
*/
const Field3D interp_to(const Field3D &var, CELL_LOC loc, REGION region) {

  SCOREP0();
  Mesh *fieldmesh = var.getMesh();
  Field3D result(fieldmesh);

  if ((loc != CELL_CENTRE && loc != CELL_DEFAULT) && (mesh->StaggerGrids == false)) {
    throw BoutException("Asked to interpolate, but StaggerGrids is disabled!");
  }
  if (fieldmesh->StaggerGrids && (var.getLocation() != loc)) {

    // Staggered grids enabled, and need to perform interpolation
    TRACE("Interpolating %s -> %s", strLocation(var.getLocation()), strLocation(loc));

    if (region != RGN_NOBNDRY) {
      // result is requested in some boundary region(s)
      result = var; // NOTE: This is just for boundaries. FIX!
    }
    result.allocate();

    // Cell location of the input field
    CELL_LOC location = var.getLocation();

    if ((location == CELL_CENTRE) || (loc == CELL_CENTRE)) {
      // Going between centred and shifted

      stencil s;
      CELL_LOC dir;

      // Get the non-centre location for interpolation direction
      dir = (loc == CELL_CENTRE) ? location : loc;

      switch (dir) {
      case CELL_XLOW: {
        ASSERT0(mesh->xstart >= 2); // At least 2 boundary cells needed for interpolation in x-direction
        BOUT_OMP(parallel)
        {
          stencil s;
          IndexOffset<Ind3D> offset(*mesh);

          //for (const auto &i : result.region(RGN_NOBNDRY)) {
          BLOCK_REGION_LOOP_PARALLEL_SECTION(result.getMesh()->getRegion3D(region), i,

            // Set stencils
            s.mm = var[offset.xmm(i)];
            s.m = var[offset.xm(i)];
            s.c = var[i];
            s.p = var[offset.xp(i)];
            s.pp = var[offset.xpp(i)];

            if ((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
              // Producing a stencil centred around a lower X value
              s.pp = s.p;
              s.p = s.c;
            } else if (location == CELL_XLOW) {
              // Stencil centred around a cell centre
              s.mm = s.m;
              s.m = s.c;
            }

            result[i] = interp(s);
          );
        }
        break;
      }
      case CELL_YLOW: {
        ASSERT0(mesh->ystart >= 2); // At least 2 boundary cells needed for interpolation in y-direction

        if (var.hasYupYdown() && ((&var.yup() != &var) || (&var.ydown() != &var))) {
          // Field "var" has distinct yup and ydown fields which
          // will be used to calculate a derivative along
          // the magnetic field
          throw BoutException("At the moment, fields with yup/ydown cannot use interp_to.\n"
                              "If we implement a 3-point stencil for interpolate or double-up\n"
                              "/double-down fields, then we can use this case.");
          BOUT_OMP(parallel)
          {
            stencil s;
            IndexOffset<Ind3D> offset(*mesh);

            //for (const auto &i : result.region(RGN_NOBNDRY)) {
            BLOCK_REGION_LOOP_PARALLEL_SECTION(result.getMesh()->getRegion3D(region), i,
              // Set stencils
              s.mm = nan("");
              s.m = var.ydown()[offset.ym(i)];
              s.c = var[i];
              s.p = var.yup()[offset.yp(i)];
              s.pp = nan("");

              if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
                // Producing a stencil centred around a lower Y value
                s.pp = s.p;
                s.p = s.c;
              } else if (location == CELL_YLOW) {
                // Stencil centred around a cell centre
                s.mm = s.m;
                s.m = s.c;
              }

              result[i] = interp(s);
            );
          }
        } else {
          // var has no yup/ydown fields, so we need to shift into field-aligned
          // coordinates

          Field3D var_fa = fieldmesh->toFieldAligned(var);
          Field3D result_fa(fieldmesh);
          result_fa.allocate();
          if (fieldmesh->ystart > 1) {

            // More than one guard cell, so set pp and mm values
            // This allows higher-order methods to be used
            BOUT_OMP(parallel)
            {
              stencil s;
              IndexOffset<Ind3D> offset(*mesh);

              //for (const auto &i : result.region(RGN_NOBNDRY)) {
              BLOCK_REGION_LOOP_PARALLEL_SECTION(result.getMesh()->getRegion3D("RGN_NOY"), i,
                // Set stencils
                s.mm = var_fa[offset.ymm(i)];
                s.m = var_fa[offset.ym(i)];
                s.c = var_fa[i];
                s.p = var_fa[offset.yp(i)];
                s.pp = var_fa[offset.ypp(i)];

                if (location == CELL_CENTRE) {
                  // Producing a stencil centred around a lower Y value
                  s.pp = s.p;
                  s.p  = s.c;
                } else {
                  // Stencil centred around a cell centre
                  s.mm = s.m;
                  s.m = s.c;
                }

                result_fa[i] = interp(s);
              );
            }
          } else {
            // Only one guard cell, so no pp or mm values
            // Note: at the moment we cannot reach this case because of the
            // 'ASSERT0(mesh->ystart >=2)' above, but if we implement a 3-point
            // stencil for interp, then this will be useful
            BOUT_OMP(parallel)
            {
              stencil s;
              IndexOffset<Ind3D> offset(*mesh);

              //for (const auto &i : result.region(RGN_NOBNDRY)) {
              BLOCK_REGION_LOOP_PARALLEL_SECTION(result.getMesh()->getRegion3D(region), i,
                // Set stencils
                s.mm = nan("");
                s.m = var_fa[offset.ym(i)];
                s.c = var_fa[i];
                s.p = var_fa[offset.yp(i)];
                s.pp = nan("");

                if (location == CELL_CENTRE) {
                  // Producing a stencil centred around a lower Y value
                  s.pp = s.p;
                  s.p = s.c;
                } else {
                  // Stencil centred around a cell centre
                  s.mm = s.m;
                  s.m = s.c;
                }

                result_fa[i] = interp(s);
              );
            }
          }
          
          result = fieldmesh->fromFieldAligned(result_fa);
        }
        break;
      }
      case CELL_ZLOW: {
        for (const auto &i : result.region(region)) {
          s.c = var[i];
          s.p = var[i.zp()];
          s.m = var[i.zm()];
          s.pp = var[i.offset(0, 0, 2)];
          s.mm = var[i.offset(0, 0, -2)];

          if (location == CELL_CENTRE) {
            // Producing a stencil centred around a lower Z value
            s.pp = s.p;
            s.p = s.c;
          } else {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m = s.c;
          }

          result[i] = interp(s);
        }
        break;
      }
      default: {
        // This should never happen
        throw BoutException("Unsupported direction of interpolation\n"
                            " - don't know how to interpolate to %s",strLocation(loc));
      }
      };

      if ((dir != CELL_ZLOW) && (region != RGN_NOBNDRY)) {
        fieldmesh->communicate(result);
      }

    } else {
      // Shifted -> shifted
      // For now, shift to centre then to final location loc
      // We probably should not rely on this, but it might work if one of the
      // shifts is in the z-direction where guard cells aren't needed.
      result = interp_to(interp_to(var, CELL_CENTRE), loc, region);
    }
    result.setLocation(loc);

    return result;
  }

  // Nothing to do - just return unchanged
  // Copying into result to return as returning var may increase the number of
  // references to the var data whilst returning result doesn't
  result = var;
  return result;
}

const Field2D interp_to(const Field2D &var, CELL_LOC UNUSED(loc), REGION UNUSED(region)) {
  SCOREP0();
  // Currently do nothing
  return var;
}

void printLocation(const Field3D &var) { output.write(strLocation(var.getLocation())); }

const char *strLocation(CELL_LOC loc) {
  switch (loc) {
  case CELL_CENTRE: {
    return " Cell centred";
  }
  case CELL_XLOW: {
    return " Lower X";
  }
  case CELL_YLOW: {
    return " Lower Y";
  }
  case CELL_ZLOW: {
    return " Lower Z";
  }
  default: { return " Default (Unknown)"; }
  };
}

// 4-point Lagrangian interpolation
// offset must be between 0 and 1
BoutReal lagrange_4pt(BoutReal v2m, BoutReal vm, BoutReal vp, BoutReal v2p,
                      BoutReal offset) {
  return -offset * (offset - 1.0) * (offset - 2.0) * v2m / 6.0 +
         0.5 * (offset * offset - 1.0) * (offset - 2.0) * vm -
         0.5 * offset * (offset + 1.0) * (offset - 2.0) * vp +
         offset * (offset * offset - 1.0) * v2p / 6.0;
}

BoutReal lagrange_4pt(BoutReal v[], BoutReal offset) {
  return lagrange_4pt(v[0], v[1], v[2], v[3], offset);
}

const Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                          const Field3D &delta_z) {
  TRACE("Interpolating 3D field");

  Mesh *mesh = f.getMesh();
  ASSERT1(mesh == delta_x.getMesh());
  ASSERT1(mesh == delta_z.getMesh());
  Field3D result(mesh);
  result.allocate();

  // Loop over output grid points
  for (int jx = 0; jx < mesh->LocalNx; jx++) {
    for (int jy = 0; jy < mesh->LocalNy; jy++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Need to get value of f at
        // [jx + delta_x[jx][jy][jz]][jy][jz + delta_z[jx][jy][jz]]

        // get lower (rounded down) index
        int jxmnew = static_cast<int>(delta_x(jx, jy, jz));
        int jzmnew = static_cast<int>(delta_z(jx, jy, jz));
        // and the distance from this point
        BoutReal xs = delta_x(jx, jy, jz) - static_cast<BoutReal>(jxmnew);
        BoutReal zs = delta_z(jx, jy, jz) - static_cast<BoutReal>(jzmnew);
        // Get new lower index
        jxmnew += jx;
        jzmnew += jz;

        // Check bounds. If beyond bounds just constant
        if (jxmnew < 0) {
          jxmnew = 0;
          xs = 0.0;
        } else if (jxmnew >= (mesh->LocalNx - 1)) {
          // Want to always be able to use [jxnew] and [jxnew+1]
          jxmnew = mesh->LocalNx - 2;
          xs = 1.0;
        }

        int jx2mnew = (jxmnew == 0) ? 0 : (jxmnew - 1);
        int jxpnew = jxmnew + 1;
        int jx2pnew = (jxmnew == (mesh->LocalNx - 2)) ? jxpnew : (jxpnew + 1);

        int ncz = mesh->LocalNz;

        // Get the 4 Z points
        jzmnew = ((jzmnew % ncz) + ncz) % ncz;
        int jzpnew = (jzmnew + 1) % ncz;
        int jz2pnew = (jzmnew + 2) % ncz;
        int jz2mnew = (jzmnew - 1 + ncz) % ncz;

        // Now have 4 indices for X and Z to interpolate

        // Interpolate in Z first
        BoutReal xvals[4];

        xvals[0] = lagrange_4pt(f(jx2mnew, jy, jz2mnew), f(jx2mnew, jy, jzmnew),
                                f(jx2mnew, jy, jzpnew), f(jx2mnew, jy, jz2pnew), zs);
        xvals[1] = lagrange_4pt(f(jxmnew, jy, jz2mnew), f(jxmnew, jy, jzmnew),
                                f(jxmnew, jy, jzpnew), f(jxmnew, jy, jz2pnew), zs);
        xvals[2] = lagrange_4pt(f(jxpnew, jy, jz2mnew), f(jxpnew, jy, jzmnew),
                                f(jxpnew, jy, jzpnew), f(jxpnew, jy, jz2pnew), zs);
        xvals[3] = lagrange_4pt(f(jx2pnew, jy, jz2mnew), f(jx2pnew, jy, jzmnew),
                                f(jx2pnew, jy, jzpnew), f(jx2pnew, jy, jz2pnew), zs);
        // Then in X
        result(jx, jy, jz) = lagrange_4pt(xvals, xs);
      }
    }
  }
  return result;
}

const Field3D interpolate(const Field2D &f, const Field3D &delta_x,
                          const Field3D &UNUSED(delta_z)) {
  SCOREP0();
  return interpolate(f, delta_x);
}

const Field3D interpolate(const Field2D &f, const Field3D &delta_x) {
  SCOREP0();
  TRACE("interpolate(Field2D, Field3D)");

  Mesh *mesh = f.getMesh();
  ASSERT1(mesh == delta_x.getMesh());
  Field3D result(mesh);
  result.allocate();

  // Loop over output grid points
  for (int jx = 0; jx < mesh->LocalNx; jx++) {
    for (int jy = 0; jy < mesh->LocalNy; jy++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Need to get value of f at
        // [jx + delta_x[jx][jy][jz]][jy][jz + delta_z[jx][jy][jz]]

        // get lower (rounded down) index
        int jxnew = static_cast<int>(delta_x(jx, jy, jz));
        // and the distance from this point
        BoutReal xs = delta_x(jx, jy, jz) - static_cast<BoutReal>(jxnew);
        // Get new lower index
        jxnew += jx;

        // Check bounds. If beyond bounds just constant
        if (jxnew < 0) {
          jxnew = 0;
          xs = 0.0;
        } else if (jxnew >= (mesh->LocalNx - 1)) {
          // Want to always be able to use [jxnew] and [jxnew+1]
          jxnew = mesh->LocalNx - 2;
          xs = 1.0;
        }
        // Interpolate in X
        result(jx, jy, jz) = f(jxnew, jy) * (1.0 - xs) + f(jxnew + 1, jy) * xs;
      }
    }
  }
  return result;
}
