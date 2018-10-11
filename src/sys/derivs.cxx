/**************************************************************************
 * Basic derivative methods
 *
 *
 * Four kinds of differencing methods:
 *
 * 1. First derivative DD*
 *    Central differencing e.g. Div(f)
 *
 * 2. Second derivatives D2D*2
 *    Central differencing e.g. Delp2(f)
 *
 * 3. Upwinding VDD*
 *    Terms like v*Grad(f)
 *
 * 4. Flux methods FDD* (e.g. flux conserving, limiting)
 *    Div(v*f)
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
#include <derivs.hxx>
#include <stencils.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <interpolation.hxx>
#include <bout/constants.hxx>
#include <msg_stack.hxx>

#include <cmath>
#include <string.h>
#include <stdlib.h>

#include <output.hxx>

/*******************************************************************************
 * First central derivatives
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

const Field3D DDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  Field3D result =  f.getMesh()->indexDDX(f,outloc, method, region);
  Coordinates *coords = f.getCoordinates(outloc);
  result /= coords->dx;

  if(f.getMesh()->IncIntShear) {
    // Using BOUT-06 style shifting
    result += coords->IntShiftTorsion * DDZ(f, outloc, method, region);
  }

  ASSERT2(((outloc == CELL_DEFAULT) && (result.getLocation() == f.getLocation())) ||
          (result.getLocation() == outloc));

  return result;
}

const Field2D DDX(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getCoordinates(outloc)->DDX(f, outloc, method, region);
}

////////////// Y DERIVATIVE /////////////////

const Field3D DDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexDDY(f, outloc, method, region) / f.getCoordinates(outloc)->dy;
}

const Field2D DDY(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getCoordinates(outloc)->DDY(f, outloc, method, region);
}

////////////// Z DERIVATIVE /////////////////

const Field3D DDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexDDZ(f, outloc, method, region) / f.getCoordinates(outloc)->dz;
}

const Field2D DDZ(const Field2D &f, CELL_LOC UNUSED(outloc), DIFF_METHOD UNUSED(method),
                  REGION UNUSED(region)) {
  auto tmp = Field2D(0., f.getMesh());
  tmp.setLocation(f.getLocation());
  return tmp;
}

const Vector3D DDZ(const Vector3D &v, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  Vector3D result(v.x.getMesh());

  ASSERT2(v.x.getMesh()==v.y.getMesh());
  ASSERT2(v.x.getMesh()==v.z.getMesh());
  Coordinates *metric = v.x.getCoordinates(outloc);

  if(v.covariant){
    // From equation (2.6.32) in D'Haeseleer
    result.x = DDZ(v.x, outloc, method, region) - v.x*metric->G1_13 - v.y*metric->G2_13 - v.z*metric->G3_13;
    result.y = DDZ(v.y, outloc, method, region) - v.x*metric->G1_23 - v.y*metric->G2_23 - v.z*metric->G3_23;
    result.z = DDZ(v.z, outloc, method, region) - v.x*metric->G1_33 - v.y*metric->G2_33 - v.z*metric->G3_33;
    result.covariant = true;
  }
  else{
    // From equation (2.6.31) in D'Haeseleer
    result.x = DDZ(v.x, outloc, method, region) + v.x*metric->G1_13 + v.y*metric->G1_23 + v.z*metric->G1_33;
    result.y = DDZ(v.y, outloc, method, region) + v.x*metric->G2_13 + v.y*metric->G2_23 + v.z*metric->G2_33;
    result.z = DDZ(v.z, outloc, method, region) + v.x*metric->G3_13 + v.y*metric->G3_23 + v.z*metric->G3_33;
    result.covariant = false;
  }

  ASSERT2(((outloc == CELL_DEFAULT) && (result.getLocation() == v.getLocation())) ||
          (result.getLocation() == outloc));

  return result;
}

const Vector2D DDZ(const Vector2D &v, CELL_LOC UNUSED(outloc), DIFF_METHOD UNUSED(method),
                   REGION UNUSED(region)) {
  Vector2D result(v.x.getMesh());

  result.covariant = v.covariant;

  // Vector 2D is constant in the z direction
  // Gx_y3 contains z-derivatives (where G is the Christoffel symbol of the
  // second kind, and x and y in {1, 2, 3})
  result.x = 0.;
  result.y = 0.;
  result.z = 0.;

  return result;
}

/*******************************************************************************
 * 2nd derivative
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

const Field3D D2DX2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  Coordinates *coords = f.getCoordinates(outloc);

  Field3D result = f.getMesh()->indexD2DX2(f, outloc, method, region) / SQ(coords->dx);
  
  if(coords->non_uniform) {
    // Correction for non-uniform f.getMesh()
    result += coords->d1_dx * f.getMesh()->indexDDX(f, outloc, DIFF_DEFAULT, region)/coords->dx;
  }

  ASSERT2(((outloc == CELL_DEFAULT) && (result.getLocation() == f.getLocation())) ||
          (result.getLocation() == outloc));

  return result;
}

const Field2D D2DX2(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  Coordinates *coords = f.getCoordinates(outloc);

  Field2D result = f.getMesh()->indexD2DX2(f, outloc, method, region) / SQ(coords->dx);

  if(coords->non_uniform) {
    // Correction for non-uniform f.getMesh()
    result += coords->d1_dx * f.getMesh()->indexDDX(f, outloc, DIFF_DEFAULT, region) / coords->dx;
  }

  return result;
}

////////////// Y DERIVATIVE /////////////////

const Field3D D2DY2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  Coordinates *coords = f.getCoordinates(outloc);

  Field3D result = f.getMesh()->indexD2DY2(f, outloc, method, region) / SQ(coords->dy);

  if(coords->non_uniform) {
    // Correction for non-uniform f.getMesh()
    result += coords->d1_dy * f.getMesh()->indexDDY(f, outloc, DIFF_DEFAULT, region) / coords->dy;
  }

  ASSERT2(((outloc == CELL_DEFAULT) && (result.getLocation() == f.getLocation())) ||
          (result.getLocation() == outloc));

  return result;
}

const Field2D D2DY2(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  Coordinates *coords = f.getCoordinates(outloc);

  Field2D result = f.getMesh()->indexD2DY2(f, outloc, method, region) / SQ(coords->dy);
  if(coords->non_uniform) {
    // Correction for non-uniform f.getMesh()
    result += coords->d1_dy * f.getMesh()->indexDDY(f, outloc, DIFF_DEFAULT, region) / coords->dy;
  }
  
  return result;
}

////////////// Z DERIVATIVE /////////////////

const Field3D D2DZ2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexD2DZ2(f, outloc, method, region) /
         SQ(f.getCoordinates(outloc)->dz);
}

const Field2D D2DZ2(const Field2D &f, CELL_LOC UNUSED(outloc), DIFF_METHOD UNUSED(method),
                    REGION UNUSED(region)) {
  auto tmp = Field2D(0., f.getMesh());
  tmp.setLocation(f.getLocation());
  return tmp;
}

/*******************************************************************************
 * Fourth derivatives
 *******************************************************************************/

const Field3D D4DX4(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexD4DX4(f, outloc, method, region) /
         SQ(SQ(f.getCoordinates(outloc)->dx));
}

const Field2D D4DX4(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexD4DX4(f, outloc, method, region) /
         SQ(SQ(f.getCoordinates(outloc)->dx));
}

const Field3D D4DY4(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexD4DY4(f, outloc, method, region) /
         SQ(SQ(f.getCoordinates(outloc)->dy));
}

const Field2D D4DY4(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexD4DY4(f, outloc, method, region) /
         SQ(SQ(f.getCoordinates(outloc)->dy));
}

const Field3D D4DZ4(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexD4DZ4(f, outloc, method, region) /
         SQ(SQ(f.getCoordinates(outloc)->dz));
}

const Field2D D4DZ4(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexD4DZ4(f, outloc, method, region) /
         SQ(SQ(f.getCoordinates(outloc)->dz));
}

/*******************************************************************************
 * Mixed derivatives
 *******************************************************************************/

/*!
 * Mixed derivative in X and Y
 *
 * This first takes derivatives in X, then in Y.
 *
 * ** Applies boundary_condition in Y, communicates
 */
const Field2D D2DXDY(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                     REGION region, string boundary_condition) {
  Field2D dfdy = DDY(f, outloc, method, RGN_NOBNDRY);
  dfdy.applyBoundary(boundary_condition);
  f.getMesh()->communicate(dfdy);
  return DDX(dfdy, outloc, method, region);
}

/*!
 * Mixed derivative in X and Y
 *
 * This first takes derivatives in X, then in Y.
 *
 * ** Applies boundary_condition in Y, communicates
 */
const Field3D D2DXDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                     REGION region, string boundary_condition) {
  Field3D dfdy = DDY(f, outloc, method, RGN_NOBNDRY);
  dfdy.applyBoundary(boundary_condition);
  f.getMesh()->communicate(dfdy);
  return DDX(dfdy, outloc, method, region);
}

/*!
 * Mixed derivative in Y and X
 *
 * This first takes derivatives in Y, then in X.
 *
 * ** Applies boundary_condition in X, communicates
 */
const Field2D D2DYDX(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                     REGION region, string boundary_condition) {
  Field2D dfdx = DDX(f, outloc, method, RGN_NOBNDRY);
  dfdx.applyBoundary(boundary_condition);
  f.getMesh()->communicate(dfdx);
  return DDY(dfdx, outloc, method, region);
}

/*!
 * Mixed derivative in Y and X
 *
 * This first takes derivatives in Y, then in X.
 *
 * ** Applies boundary_condition in X, communicates
 */
const Field3D D2DYDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                     REGION region, string boundary_condition) {
  // if free_o3 boundary condition is used, it might be better to have the
  // boundary condition after the communication because there might not be
  // enough grid points for the 3-point stencil. Then we would need to call
  // calcYupYdown again somehow.
  Field3D dfdx = DDX(f, outloc, method, RGN_NOBNDRY);
  dfdx.applyBoundary(boundary_condition);
  f.getMesh()->communicate(dfdx);
  return DDY(dfdx, outloc, method, region);
}

const Field2D D2DXDZ(const Field2D &f, CELL_LOC UNUSED(outloc),
                     DIFF_METHOD UNUSED(method), REGION UNUSED(region),
                     string UNUSED(boundary_condition)) {
  auto tmp = Field2D(0., f.getMesh());
  tmp.setLocation(f.getLocation());
  return tmp;
}

/// X-Z mixed derivative
const Field3D D2DXDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                     REGION region, string UNUSED(boundary_condition)) {
  // Take derivative in Z, including in X boundaries. Then take derivative in X
  // Maybe should average results of DDX(DDZ) and DDZ(DDX)?
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  // region specifies what the combined derivative should return
  // Therefore we need to add the X boundary to the inner derivative
  // RGN_NOY and RGN_NOZ include the X boundary, therefore we need to
  // throw - or add communication code.
  REGION region_inner;
  switch (region){
  case RGN_NOBNDRY:
    region_inner = RGN_NOY;
    break;
  case RGN_NOX:
    region_inner = RGN_ALL;
    break;
  default:
    throw BoutException("Unhandled region case in D2DXDZ");
  }

  return DDX(DDZ(f, outloc,method, region_inner),outloc,method,region);
}

const Field2D D2DZDX(const Field2D &f, CELL_LOC UNUSED(outloc),
                     DIFF_METHOD UNUSED(method), REGION UNUSED(region),
                     string UNUSED(boundary_condition)) {
  auto tmp = Field2D(0., f.getMesh());
  tmp.setLocation(f.getLocation());
  return tmp;
}

/// Z-X mixed derivative
const Field3D D2DZDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                     REGION region, string UNUSED(boundary_condition)) {
  // Take derivative in X, then take derivative in Z
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());

  return DDZ(DDX(f, outloc, method, region), outloc, method, region);
}

const Field2D D2DYDZ(const Field2D &f, CELL_LOC UNUSED(outloc),
                     DIFF_METHOD UNUSED(method), REGION UNUSED(region),
                     string UNUSED(boundary_condition)) {
  auto tmp = Field2D(0., f.getMesh());
  tmp.setLocation(f.getLocation());
  return tmp;
}

const Field3D D2DYDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                     REGION region, string UNUSED(boundary_condition)) {
  // Take derivative in Z, including in Y boundaries. Then take derivative in Y
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  // region specifies what the combined derivative should return
  // Therefore we need to add the Y boundary to the inner derivative
  // RGN_NOX and RGN_NOZ include the Y boundary, therefore we need to
  // throw - or add communication code.
  REGION region_inner;
  switch (region){
  case RGN_NOBNDRY:
    region_inner = RGN_NOX;
    break;
  case RGN_NOY:
    region_inner = RGN_ALL;
    break;
  default:
    throw BoutException("Unhandled region case in D2DYDZ");
  }

  Field3D dfdz = DDZ(f, outloc, method, region_inner);

  if (f.hasYupYdown()) {
    // Need to set yup/ydown fields for dfdz
    dfdz.splitYupYdown();
    dfdz.yup() = DDZ(f.yup(), outloc, method, region_inner);
    dfdz.ydown() = DDZ(f.ydown(), outloc, method, region_inner);
  }

  return DDY(dfdz, outloc, method, region);
}

const Field2D D2DZDY(const Field2D &f, CELL_LOC UNUSED(outloc),
                     DIFF_METHOD UNUSED(method), REGION UNUSED(region),
                     string UNUSED(boundary_condition)) {
  auto tmp = Field2D(0., f.getMesh());
  tmp.setLocation(f.getLocation());
  return tmp;
}

const Field3D D2DZDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                     REGION region, string UNUSED(boundary_condition)) {
  // Take derivative in Y, then take derivative in Z
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());

  return DDZ(DDY(f, outloc,method, region), outloc, method, region);
}

/*******************************************************************************
 * Advection schemes
 *
 * Jan 2009  - Re-written to use Set*Stencil routines
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

/// Special case where both arguments are 2D. Output location ignored for now
const Field2D VDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexVDDX(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dx;
}

/// General version for 2 or 3-D objects
const Field3D VDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexVDDX(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dx;
}

////////////// Y DERIVATIVE /////////////////

// special case where both are 2D
const Field2D VDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexVDDY(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dy;
}

// general case
const Field3D VDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexVDDY(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dy;
}

////////////// Z DERIVATIVE /////////////////

// special case where both are 2D
const Field2D VDDZ(const Field2D &v, const Field2D &UNUSED(f), CELL_LOC UNUSED(outloc),
                   DIFF_METHOD UNUSED(method), REGION UNUSED(region)) {
  // Should we take location from v or f?
  auto tmp = Field2D(0., v.getMesh());
  tmp.setLocation(v.getLocation());
  return tmp;
}

// Note that this is zero because no compression is included
const Field2D VDDZ(const Field3D &v, const Field2D &UNUSED(f), CELL_LOC UNUSED(outloc),
                   DIFF_METHOD UNUSED(method), REGION UNUSED(region)) {
  // Should we take location from v or f?
  auto tmp = Field2D(0., v.getMesh());
  tmp.setLocation(v.getLocation());
  return tmp;
}

// general case
const Field3D VDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexVDDZ(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dz;
}

/*******************************************************************************
 * Flux conserving schemes
 *******************************************************************************/
const Field2D FDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexFDDX(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dx;
}

const Field3D FDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexFDDX(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dx;
}

/////////////////////////////////////////////////////////////////////////

const Field2D FDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexFDDY(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dy;
}

const Field3D FDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexFDDY(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dy;
}

/////////////////////////////////////////////////////////////////////////

const Field2D FDDZ(const Field2D &v, const Field2D &UNUSED(f), CELL_LOC UNUSED(outloc),
                   DIFF_METHOD UNUSED(method), REGION UNUSED(region)) {
  // Should we take location from v or f?
  auto tmp = Field2D(0., v.getMesh());
  tmp.setLocation(v.getLocation());
  return tmp;
}

const Field3D FDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return f.getMesh()->indexFDDZ(v, f, outloc, method, region) /
         f.getCoordinates(outloc)->dz;
}
