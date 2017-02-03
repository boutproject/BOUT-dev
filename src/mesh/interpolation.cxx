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
  if(mesh->StaggerGrids && (var.getLocation() != loc)) {
    
    //output.write("\nINTERPOLATING %s -> %s\n", strLocation(var.getLocation()), strLocation(loc));

    // Staggered grids enabled, and need to perform interpolation

#ifdef CHECK
    msg_stack.push("Interpolating %s -> %s", strLocation(var.getLocation()), strLocation(loc));
#endif

    Field3D result;

    result = var; // NOTE: This is just for boundaries. FIX!
    result.allocate();
    
    if((var.getLocation() == CELL_CENTRE) || (loc == CELL_CENTRE)) {
      // Going between centred and shifted
      
      bindex bx;
      stencil s;
      CELL_LOC dir; 
      
      // Get the non-centre location for interpolation direction
      dir = (loc == CELL_CENTRE) ? var.getLocation() : loc;

      switch(dir) {
      case CELL_XLOW: {
	start_index(&bx, RGN_NOX);
	do {
	  var.setXStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = interp(s);
	}while(next_index3(&bx));
	break;
	// Need to communicate in X
      }
      case CELL_YLOW: {
	start_index(&bx, RGN_NOY);
	do {
	  var.setYStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = interp(s);
	}while(next_index3(&bx));
	break;
	// Need to communicate in Y
      }
      case CELL_ZLOW: {
	start_index(&bx, RGN_NOZ);
	do {
	  var.setZStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = interp(s);
	}while(next_index3(&bx));
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

#ifdef CHECK
    msg_stack.pop();
#endif

    return result;
  }
  
  // Nothing to do - just return unchanged
  return var;
}

const Field2D interp_to(const Field2D &var, CELL_LOC UNUSED(loc)) {
  // Currently do nothing
  return var;
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

const Field3D interpolate(const Field3D &f, const Field3D &delta_x, const Field3D &delta_z) {
  TRACE("Interpolating 3D field");
  
  Field3D result;
  result.allocate();

  // Loop over output grid points
  for(int jx=0;jx<mesh->LocalNx;jx++)
    for(int jy=0;jy<mesh->LocalNy;jy++)
      for(int jz=0;jz<mesh->LocalNz;jz++) {
	// Need to get value of f at 
	// [jx + delta_x[jx][jy][jz]][jy][jz + delta_z[jx][jy][jz]]

	// get lower (rounded down) index
	int jxmnew = (int) delta_x(jx,jy,jz);
	int jzmnew = (int) delta_z(jx,jy,jz);
	// and the distance from this point
	BoutReal xs = delta_x(jx,jy,jz) - ((BoutReal) jxmnew);
	BoutReal zs = delta_z(jx,jy,jz) - ((BoutReal) jzmnew);
	// Get new lower index
	jxmnew += jx; jzmnew += jz;
	
	// Check bounds. If beyond bounds just constant
	if(jxmnew < 0) {
	  jxmnew = 0;
	  xs = 0.0;
	}else if(jxmnew >= (mesh->LocalNx-1)) {
	  // Want to always be able to use [jxnew] and [jxnew+1]
	  jxmnew = mesh->LocalNx-2; 
	  xs = 1.0;
	}

	int jx2mnew = (jxmnew == 0) ? 0 : (jxmnew - 1);
	int jxpnew = jxmnew + 1;
	int jx2pnew = (jxmnew == (mesh->LocalNx-2)) ? jxpnew : (jxpnew + 1);

	int ncz = mesh->LocalNz;

	// Get the 4 Z points
	jzmnew = ((jzmnew % ncz) + ncz) % ncz;
	int jzpnew = (jzmnew + 1) % ncz;
	int jz2pnew = (jzmnew + 2) % ncz;
	int jz2mnew = (jzmnew - 1 + ncz) % ncz;

	// Now have 4 indices for X and Z to interpolate
	
	// Interpolate in Z first
	BoutReal xvals[4];
	
	xvals[0] = lagrange_4pt(f(jx2mnew,jy,jz2mnew),
				f(jx2mnew,jy,jzmnew),
				f(jx2mnew,jy,jzpnew),
				f(jx2mnew,jy,jz2pnew),
				zs);
	xvals[1] = lagrange_4pt(f(jxmnew,jy,jz2mnew),
				f(jxmnew,jy,jzmnew),
				f(jxmnew,jy,jzpnew),
				f(jxmnew,jy,jz2pnew),
				zs);
	xvals[2] = lagrange_4pt(f(jxpnew,jy,jz2mnew),
				f(jxpnew,jy,jzmnew),
				f(jxpnew,jy,jzpnew),
				f(jxpnew,jy,jz2pnew),
				zs);
	xvals[3] = lagrange_4pt(f(jx2pnew,jy,jz2mnew),
				f(jx2pnew,jy,jzmnew),
				f(jx2pnew,jy,jzpnew),
				f(jx2pnew,jy,jz2pnew),
				zs);
	// Then in X
	result(jx,jy,jz) = lagrange_4pt(xvals, xs);
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
  for(int jx=0;jx<mesh->LocalNx;jx++)
    for(int jy=0;jy<mesh->LocalNy;jy++)
      for(int jz=0;jz<mesh->LocalNz;jz++) {
	// Need to get value of f at 
	// [jx + delta_x[jx][jy][jz]][jy][jz + delta_z[jx][jy][jz]]
	
	// get lower (rounded down) index
	int jxnew = (int) delta_x(jx,jy,jz);
	// and the distance from this point
	BoutReal xs = delta_x(jx,jy,jz) - ((BoutReal) jxnew);
	// Get new lower index
	jxnew += jx;
	
	// Check bounds. If beyond bounds just constant
	if(jxnew < 0) {
	  jxnew = 0;
	  xs = 0.0;
	}else if(jxnew >= (mesh->LocalNx-1)) {
	  // Want to always be able to use [jxnew] and [jxnew+1]
	  jxnew = mesh->LocalNx-2; 
	  xs = 1.0;
	}
	// Interpolate in X
	result(jx,jy,jz) = f(jxnew,jy)*(1.0 - xs) + f(jxnew+1,jy)*xs;
      }

  return result;
}
