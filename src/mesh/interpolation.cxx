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

  Mesh *fieldmesh = var.getMesh();
  Field3D result(fieldmesh);

  if ((loc != CELL_CENTRE && loc != CELL_DEFAULT) && (fieldmesh->StaggerGrids == false)) {
    throw BoutException("Asked to interpolate, but StaggerGrids is disabled!");
  }
  if (fieldmesh->StaggerGrids && (var.getLocation() != loc)) {

    // Staggered grids enabled, and need to perform interpolation
    TRACE("Interpolating %s -> %s", strLocation(var.getLocation()), strLocation(loc));

    if (region != RGN_NOBNDRY) {
      // result is requested in some boundary region(s)
      result = var; // NOTE: This is just for boundaries. FIX!
    }
    // NOTE: invalidateGuards() is called in Field3D::alloctate() if the data
    // block is not already allocated, so will be called here if
    // region==RGN_NOBNDRY
    result.allocate();

    // Cell location of the input field
    CELL_LOC location = var.getLocation();

    if ((location == CELL_CENTRE) || (loc == CELL_CENTRE)) {
      // Going between centred and shifted
      CELL_LOC dir;

      // Get the non-centre location for interpolation direction
      dir = (loc == CELL_CENTRE) ? location : loc;

      switch (dir) {
      case CELL_XLOW: {
        // At least 2 boundary cells needed for interpolation in x-direction
        ASSERT0(fieldmesh->xstart >= 2);

        BOUT_OMP(parallel) {
          stencil s;
          BOUT_FOR_INNER(i, result.getRegion("RGN_NOBNDRY")) {
            // Set stencils
            s.mm = var[i.xmm()];
            s.m = var[i.xm()];
            s.c = var[i];
            s.p = var[i.xp()];
            s.pp = var[i.xpp()];

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
          }
        }
        break;
      }
      case CELL_YLOW: {
        // At least 2 boundary cells needed for interpolation in y-direction
        ASSERT0(fieldmesh->ystart >= 2);

        if (var.hasYupYdown() && ((&var.yup() != &var) || (&var.ydown() != &var))) {
          // Field "var" has distinct yup and ydown fields which
          // will be used to calculate a derivative along
          // the magnetic field
          throw BoutException("At the moment, fields with yup/ydown cannot use interp_to.\n"
                              "If we implement a 3-point stencil for interpolate or double-up\n"
                              "/double-down fields, then we can use this case.");
          BOUT_OMP(parallel) {
            stencil s;
            BOUT_FOR_INNER(i, result.getRegion("RGN_NOBNDRY")) {
              // Set stencils
              s.m = var.ydown()[i.ym()];
              s.c = var[i];
              s.p = var.yup()[i.yp()];

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
            }
          }
        } else {
          // var has no yup/ydown fields, so we need to shift into field-aligned
          // coordinates

          Field3D var_fa = fieldmesh->toFieldAligned(var);
          if (region != RGN_NOBNDRY) {
            // repeat the hack above for boundary points
            // this avoids a duplicate toFieldAligned call if we had called
            // result = toFieldAligned(result)
            // to get the boundary cells
            //
            // result is requested in some boundary region(s)
            result = var_fa; // NOTE: This is just for boundaries. FIX!
            result.allocate();
          }
          result.setCoordinateSystem(CoordinateSystem::FieldAligned);
          if (fieldmesh->ystart > 1) {

            // More than one guard cell, so set pp and mm values
            // This allows higher-order methods to be used
            BOUT_OMP(parallel) {
              stencil s;
              BOUT_FOR_INNER(i, result.getRegion("RGN_NOBNDRY")) {
                // Set stencils
                s.mm = var_fa[i.ymm()];
                s.m = var_fa[i.ym()];
                s.c = var_fa[i];
                s.p = var_fa[i.yp()];
                s.pp = var_fa[i.ypp()];

                if (location == CELL_CENTRE) {
                  // Producing a stencil centred around a lower Y value
                  s.pp = s.p;
                  s.p = s.c;
                } else {
                  // Stencil centred around a cell centre
                  s.mm = s.m;
                  s.m = s.c;
                }

                result[i] = interp(s);
              }
            }
          } else {
            // Only one guard cell, so no pp or mm values
            // Note: at the moment we cannot reach this case because of the
            // 'ASSERT0(fieldmesh->ystart >=2)' above, but if we implement a 3-point
            // stencil for interp, then this will be useful
            BOUT_OMP(parallel) {
              stencil s;
              BOUT_FOR_INNER(i, result.getRegion("RGN_NOBNDRY")) {
                // Set stencils
                s.m = var_fa[i.ym()];
                s.c = var_fa[i];
                s.p = var_fa[i.yp()];

                if (location == CELL_CENTRE) {
                  // Producing a stencil centred around a lower Y value
                  s.pp = s.p;
                  s.p = s.c;
                } else {
                  // Stencil centred around a cell centre
                  s.mm = s.m;
                  s.m = s.c;
                }

                result[i] = interp(s);
              }
            }
          }
          
          result = fieldmesh->fromFieldAligned(result, RGN_NOBNDRY);
        }
        break;
      }
      case CELL_ZLOW: {
        /// Convert REGION enum to a Region string identifier
        const auto region_str = REGION_STRING(region);

        BOUT_OMP(parallel) {
          stencil s;
          BOUT_FOR_INNER(i, result.getRegion(region_str)) {
            s.mm = var[i.zmm()];
            s.m = var[i.zm()];
            s.c = var[i];
            s.p = var[i.zp()];
            s.pp = var[i.zpp()];

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

const Field2D interp_to(const Field2D &var, CELL_LOC loc, REGION region) {

  Mesh *fieldmesh = var.getMesh();
  Field2D result(fieldmesh);

  if ((loc != CELL_CENTRE && loc != CELL_DEFAULT) && (fieldmesh->StaggerGrids == false)) {
    throw BoutException("Asked to interpolate, but StaggerGrids is disabled!");
  }
  if (fieldmesh->StaggerGrids && (var.getLocation() != loc)) {

    // Staggered grids enabled, and need to perform interpolation
    TRACE("Interpolating %s -> %s", strLocation(var.getLocation()), strLocation(loc));

    if (region != RGN_NOBNDRY) {
      // result is requested in some boundary region(s)
      result = var; // NOTE: This is just for boundaries. FIX!
    }
    // NOTE: invalidateGuards() is called in Field3D::alloctate() if the data
    // block is not already allocated, so will be called here if
    // region==RGN_NOBNDRY
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
        ASSERT0(fieldmesh->xstart >= 2); // At least 2 boundary cells needed for interpolation in x-direction

        BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {

          // Set stencils
          s.c = var[i];
          s.p = var[i.xp()];
          s.m = var[i.xm()];
          s.pp = var[i.offset(2, 0, 0)];
          s.mm = var[i.offset(-2, 0, 0)];

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
        }
        break;
      }
      case CELL_YLOW: {
        ASSERT0(fieldmesh->ystart >= 2); // At least 2 boundary cells needed for interpolation in y-direction

        if (fieldmesh->ystart > 1) {

          // More than one guard cell, so set pp and mm values
          // This allows higher-order methods to be used
          BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
            // Set stencils
            s.c = var[i];
            s.p = var[i.yp()];
            s.m = var[i.ym()];
            s.pp = var[i.offset(0, 2, 0)];
            s.mm = var[i.offset(0, -2, 0)];

            if (location == CELL_CENTRE) {
              // Producing a stencil centred around a lower Y value
              s.pp = s.p;
              s.p  = s.c;
              } else {
                // Stencil centred around a cell centre
                s.mm = s.m;
                s.m = s.c;
              }

            result[i] = interp(s);
          }
        } else {
          // Only one guard cell, so no pp or mm values
          // Note: at the moment we cannot reach this case because of the
          // 'ASSERT0(fieldmesh->ystart >=2)' above, but if we implement a 3-point
          // stencil for interp, then this will be useful
          BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
            // Set stencils
            s.c = var[i];
            s.p = var[i.yp()];
            s.m = var[i.ym()];

            if (location == CELL_CENTRE) {
              // Producing a stencil centred around a lower Y value
              s.pp = s.p;
              s.p = s.c;
            } else {
              // Stencil centred around a cell centre
              s.mm = s.m;
              s.m = s.c;
            }

            result[i] = interp(s);
          }
        }
        break;
      }
      case CELL_ZLOW: {
        // Nothing to do for Field2D as Field2D is constant in z-direction
        result = var;
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

void printLocation(const Field3D &var) { output.write(strLocation(var.getLocation())); }
void printLocation(const Field2D &var) { output.write(strLocation(var.getLocation())); }

const char *strLocation(CELL_LOC loc) { return CELL_LOC_STRING(loc).c_str(); }

const Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                          const Field3D &delta_z) {
  TRACE("Interpolating 3D field");
  Lagrange4pt interpolateMethod{f.getMesh()};
  return interpolateMethod.interpolate(f, delta_x, delta_z);
}

const Field3D interpolate(const Field2D &f, const Field3D &delta_x,
                          const Field3D &UNUSED(delta_z)) {
  return interpolate(f, delta_x);
}

const Field3D interpolate(const Field2D &f, const Field3D &delta_x) {
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
