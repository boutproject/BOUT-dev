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
#include <stencils.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <unused.hxx>

/// Perform interpolation between centre -> shifted or vice-versa
/*!
  Interpolate using 4th-order staggered formula
  
  @param[in] s  Input stencil. mm -> -3/2, m -> -1/2, p -> +1/2, pp -> +3/2
*/
BoutReal interp(const stencil &s)
{
  return ( 9.*(s.m + s.p) - s.mm - s.pp ) / 16.;
}

/*!
  Interpolate between different cell locations
  
  NOTE: This requires communication

  @param[in]   var  Input variable
  @param[in]   loc  Location of output values
*/
const Field3D interp_to(const Field3D &var, CELL_LOC loc)
{
  REGION region = RGN_INTERP;

  if(mesh->StaggerGrids && (var.getLocation() != loc)) {

    // Staggered grids enabled, and need to perform interpolation
    TRACE("Interpolating %s -> %s", strLocation(var.getLocation()), strLocation(loc));

    Field3D result;

    result = var; // NOTE: This is just for boundaries. FIX!
    result.allocate();
    result.setLocation(loc); // Set the result location
    
    if((var.getLocation() == CELL_CENTRE) || (loc == CELL_CENTRE)) {
      // Going between centred and shifted
      
      bindex bx;
      CELL_LOC dir; 
      
      // Get the non-centre location for interpolation direction
      dir = (loc == CELL_CENTRE) ? var.getLocation() : loc;

      switch(dir) {
      case CELL_XLOW: {
        if (loc == CELL_XLOW) {
          for(const auto &i : result.region(region)) {
            // Producing a stencil centred around a lower X value
            stencil s;
            s.p = var[i];
            s.m = var[i.xm()];
            s.pp = var[i.xp()];
            s.mm = var[i.offset(-2,0,0)];
            result[i] = interp(s);
          }
        } else {
          for(const auto &i : result.region(region)) {
            // Producing a stencil centred around a cell centre
            stencil s;
            s.p = var[i.xp()];
            s.m = var[i];
            s.pp = var[i.offset(2,0,0)];
            s.mm = var[i.xm()];
            result[i] = interp(s);
          }
        }
	break;
	// Need to communicate in X
      }
      case CELL_YLOW: {
        if(var.hasYupYdown() &&
            ( (&var.yup() != &var) || (&var.ydown() != &var) )) {
          // Field "var" has distinct yup and ydown fields which
          // will should be used to calculate interpolation along
          // the magnetic field
          throw BoutException("interp_to not implemented for fields with yup/ydown");
        } else {
          Field3D var_fa = mesh->toFieldAligned(var);
          if (loc == CELL_YLOW) {
            for(const auto &i : result.region(region)) {
              // Producing a stencil centred around a lower Y value
              stencil s;
              s.p = var[i];
              s.m = var[i.ym()];
              s.pp = var[i.yp()];
              s.mm = var[i.offset(0,-2,0)];
              result[i] = interp(s);
            }
          } else {
            for(const auto &i : result.region(region)) {
              // Producing a stencil centred around a cell centre
              stencil s;
              s.p = var[i.yp()];
              s.m = var[i];
              s.pp = var[i.offset(0,2,0)];
              s.mm = var[i.ym()];
              result[i] = interp(s);
            }
          }
        }
	break;
	// Need to communicate in Y
      }
      case CELL_ZLOW: {
	if (loc == CELL_ZLOW) {
          for(const auto &i : result.region(region)) {
            // Producing a stencil centred around a lower Z value
            stencil s;
            s.p = var[i];
            s.m = var[i.zm()];
            s.pp = var[i.zp()];
            s.mm = var[i.offset(0,0,-2)];
            result[i] = interp(s);
          }
        } else {
          for(const auto &i : result.region(region)) {
            // Producing a stencil centred around a cell centre
            stencil s;
            s.p = var[i.zp()];
            s.m = var[i];
            s.pp = var[i.offset(0,0,2)];
            s.mm = var[i.zm()];
            result[i] = interp(s);
          }
        }
	break;
      }
      default: {
	// This should never happen
	throw BoutException("Don't know what to do");
      }
      };
      
      if(dir != CELL_ZLOW) {
	// COMMUNICATION
	
	mesh->communicate(result);

	// BOUNDARIES

      }

    }else {
      // Shifted -> shifted
      // For now, shift to centre then to shifted
      
      result = interp_to( interp_to(var, CELL_CENTRE) , loc);
    }
    result.setLocation(loc);

    return result;
  }
  
  // Nothing to do - just return unchanged
  return var;
}

const Field2D interp_to(const Field2D &var, CELL_LOC loc) {
  // Throw exception if something needs to be done
  if (loc == CELL_DEFAULT || var.getLocation() == loc || (loc == CELL_ZLOW && var.getLocation() == CELL_CENTRE) || (loc == CELL_CENTRE && var.getLocation() == CELL_ZLOW)) {
    // Nothing needs to be done for Field2D if var is already at loc, or if the interpolation would be in the z-direction (since a Field2D is axi-symmetric)
    return var;
  } else {
    throw BoutException("interp_to is not currently implemented for Field2D unless nothing needs to be done");
  }
}

void printLocation(const Field3D &var) {
  output.write(strLocation(var.getLocation()));
}

const char* strLocation(CELL_LOC loc) {
  switch(loc) {
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
  default: {
    return " Default (Unknown)";
  }
  };
}

// 4-point Lagrangian interpolation
// offset must be between 0 and 1
BoutReal lagrange_4pt(BoutReal v2m, BoutReal vm, BoutReal vp, BoutReal v2p, BoutReal offset)
{
  return -offset*(offset-1.0)*(offset-2.0)*v2m/6.0
    + 0.5*(offset*offset - 1.0)*(offset-2.0)*vm
    - 0.5*offset*(offset+1.0)*(offset-2.0)*vp
    + offset*(offset*offset - 1.0)*v2p/6.0;
}

BoutReal lagrange_4pt(BoutReal v[], BoutReal offset)
{
  return lagrange_4pt(v[0], v[1], v[2], v[3], offset);
}

const Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                          const Field3D &delta_z) {
  TRACE("Interpolating 3D field");

  Field3D result;
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

const Field3D interpolate(const Field2D &f, const Field3D &delta_x, const Field3D &UNUSED(delta_z)) {
  return interpolate(f, delta_x);
}

const Field3D interpolate(const Field2D &f, const Field3D &delta_x) {
  TRACE("interpolate(Field2D, Field3D)");

  Field3D result;
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
