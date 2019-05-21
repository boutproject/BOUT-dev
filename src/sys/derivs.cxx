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

#include <output.hxx>
#include <unused.hxx>

/*******************************************************************************
 * First central derivatives
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

const Field3D DDX(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  Field3D result = bout::derivatives::index::DDX(f, outloc, method, region);
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

const Field2D DDX(const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return f.getCoordinates(outloc)->DDX(f, outloc, method, region);
}

////////////// Y DERIVATIVE /////////////////

const Field3D DDY(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::DDY(f, outloc, method, region)
         / f.getCoordinates(outloc)->dy;
}

const Field2D DDY(const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return f.getCoordinates(outloc)->DDY(f, outloc, method, region);
}

////////////// Z DERIVATIVE /////////////////

const Field3D DDZ(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::DDZ(f, outloc, method, region)
         / f.getCoordinates(outloc)->dz;
}

const Field2D DDZ(const Field2D &f, CELL_LOC UNUSED(outloc), const std::string &UNUSED(method),
                  REGION UNUSED(region)) {
  auto tmp = Field2D(0., f.getMesh());
  tmp.setLocation(f.getLocation());
  return tmp;
}

const Vector3D DDZ(const Vector3D &v, CELL_LOC outloc, const std::string &method, REGION region) {
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

const Vector2D DDZ(const Vector2D &v, CELL_LOC UNUSED(outloc), const std::string &UNUSED(method),
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

const Field3D D2DX2(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  Coordinates *coords = f.getCoordinates(outloc);

  Field3D result =
      bout::derivatives::index::D2DX2(f, outloc, method, region) / SQ(coords->dx);

  if(coords->non_uniform) {
    // Correction for non-uniform f.getMesh()
    result += coords->d1_dx * bout::derivatives::index::DDX(f, outloc, "DEFAULT", region)
              / coords->dx;
  }

  ASSERT2(((outloc == CELL_DEFAULT) && (result.getLocation() == f.getLocation())) ||
          (result.getLocation() == outloc));

  return result;
}

const Field2D D2DX2(const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  Coordinates *coords = f.getCoordinates(outloc);

  Field2D result =
      bout::derivatives::index::D2DX2(f, outloc, method, region) / SQ(coords->dx);

  if(coords->non_uniform) {
    // Correction for non-uniform f.getMesh()
    result += coords->d1_dx * bout::derivatives::index::DDX(f, outloc, "DEFAULT", region)
              / coords->dx;
  }

  return result;
}

////////////// Y DERIVATIVE /////////////////

const Field3D D2DY2(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  Coordinates *coords = f.getCoordinates(outloc);

  Field3D result =
      bout::derivatives::index::D2DY2(f, outloc, method, region) / SQ(coords->dy);

  if(coords->non_uniform) {
    // Correction for non-uniform f.getMesh()
    result += coords->d1_dy * bout::derivatives::index::DDY(f, outloc, "DEFAULT", region)
              / coords->dy;
  }

  ASSERT2(((outloc == CELL_DEFAULT) && (result.getLocation() == f.getLocation())) ||
          (result.getLocation() == outloc));

  return result;
}

const Field2D D2DY2(const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  Coordinates *coords = f.getCoordinates(outloc);

  Field2D result =
      bout::derivatives::index::D2DY2(f, outloc, method, region) / SQ(coords->dy);
  if(coords->non_uniform) {
    // Correction for non-uniform f.getMesh()
    result += coords->d1_dy * bout::derivatives::index::DDY(f, outloc, "DEFAULT", region)
              / coords->dy;
  }
  
  return result;
}

////////////// Z DERIVATIVE /////////////////

const Field3D D2DZ2(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::D2DZ2(f, outloc, method, region)
         / SQ(f.getCoordinates(outloc)->dz);
}

const Field2D D2DZ2(const Field2D &f, CELL_LOC outloc, const std::string &UNUSED(method),
                    REGION UNUSED(region)) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return zeroFrom(f).setLocation(outloc);
}

/*******************************************************************************
 * Fourth derivatives
 *******************************************************************************/

const Field3D D4DX4(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::D4DX4(f, outloc, method, region)
         / SQ(SQ(f.getCoordinates(outloc)->dx));
}

const Field2D D4DX4(const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::D4DX4(f, outloc, method, region)
         / SQ(SQ(f.getCoordinates(outloc)->dx));
}

const Field3D D4DY4(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::D4DY4(f, outloc, method, region)
         / SQ(SQ(f.getCoordinates(outloc)->dy));
}

const Field2D D4DY4(const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::D4DY4(f, outloc, method, region)
         / SQ(SQ(f.getCoordinates(outloc)->dy));
}

const Field3D D4DZ4(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::D4DZ4(f, outloc, method, region)
         / SQ(SQ(f.getCoordinates(outloc)->dz));
}

const Field2D D4DZ4(const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::D4DZ4(f, outloc, method, region)
         / SQ(SQ(f.getCoordinates(outloc)->dz));
}

/*******************************************************************************
 * Mixed derivatives
 *******************************************************************************/

/*!
 * Mixed derivative in X and Y
 *
 * This first takes derivatives in Y, then in X.
 *
 * ** Communicates and applies boundary in X.
 */
const Field2D D2DXDY(const Field2D& f, CELL_LOC outloc, const std::string& method,
    REGION region, const std::string& dfdy_boundary_condition) {

  // If staggering in x, take y-derivative at f's location.
  const auto y_location =
    (outloc == CELL_XLOW or f.getLocation() == CELL_XLOW) ? CELL_DEFAULT : outloc;

  Field2D dfdy = DDY(f, y_location, method, region);

  // Set x-guard cells and x-boundary cells before calculating DDX
  f.getMesh()->communicate(dfdy);
  dfdy.applyBoundary(dfdy_boundary_condition);

  return DDX(dfdy, outloc, method, region);
}

/*!
 * Mixed derivative in X and Y
 *
 * This first takes derivatives in Y, then in X.
 *
 * ** Communicates and applies boundary in X.
 */
const Field3D D2DXDY(const Field3D& f, CELL_LOC outloc, const std::string& method,
    REGION region, const std::string& dfdy_boundary_condition) {

  // If staggering in x, take y-derivative at f's location.
  const auto y_location =
    (outloc == CELL_XLOW or f.getLocation() == CELL_XLOW) ? CELL_DEFAULT : outloc;

  Field3D dfdy = DDY(f, y_location, method, region);

  // Set x-guard cells and x-boundary cells before calculating DDX
  f.getMesh()->communicate(dfdy);
  dfdy.applyBoundary(dfdy_boundary_condition);

  return DDX(dfdy, outloc, method, region);
}

const Field2D D2DXDZ(const Field2D &f, CELL_LOC outloc,
                     const std::string &UNUSED(method), REGION UNUSED(region)) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return zeroFrom(f).setLocation(outloc);
}

/// X-Z mixed derivative
const Field3D D2DXDZ(const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  // If staggering in z, take x-derivative at f's location.
  const auto x_location =
    (outloc == CELL_ZLOW or f.getLocation() == CELL_ZLOW) ? CELL_DEFAULT : outloc;

  return DDZ(DDX(f, x_location, method, region), outloc, method, region);
}

const Field2D D2DYDZ(const Field2D &f, CELL_LOC outloc,
                     const std::string &UNUSED(method), REGION UNUSED(region)) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return zeroFrom(f).setLocation(outloc);
}

const Field3D D2DYDZ(const Field3D& f, CELL_LOC outloc, const std::string& method,
                     REGION region) {
  // If staggering in z, take y-derivative at f's location.
  const auto y_location =
    (outloc == CELL_ZLOW or f.getLocation() == CELL_ZLOW) ? CELL_DEFAULT : outloc;

  return DDZ(DDY(f, y_location, method, region), outloc, method, region);
}

/*******************************************************************************
 * Advection schemes
 *
 * Jan 2009  - Re-written to use Set*Stencil routines
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

/// Special case where both arguments are 2D. Output location ignored for now
const Field2D VDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::VDDX(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dx;
}

/// General version for 2 or 3-D objects
const Field3D VDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::VDDX(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dx;
}

////////////// Y DERIVATIVE /////////////////

// special case where both are 2D
const Field2D VDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::VDDY(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dy;
}

// general case
const Field3D VDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::VDDY(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dy;
}

////////////// Z DERIVATIVE /////////////////

// special case where both are 2D
const Field2D VDDZ(const Field2D &UNUSED(v), const Field2D &f, CELL_LOC outloc,
                   const std::string &UNUSED(method), REGION UNUSED(region)) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return zeroFrom(f).setLocation(outloc);
}

// Note that this is zero because no compression is included
const Field2D VDDZ(const Field3D &UNUSED(v), const Field2D &f, CELL_LOC outloc,
                   const std::string &UNUSED(method), REGION UNUSED(region)) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return zeroFrom(f).setLocation(outloc);
}

// general case
const Field3D VDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::VDDZ(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dz;
}

/*******************************************************************************
 * Flux conserving schemes
 *******************************************************************************/
const Field2D FDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::FDDX(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dx;
}

const Field3D FDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::FDDX(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dx;
}

/////////////////////////////////////////////////////////////////////////

const Field2D FDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::FDDY(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dy;
}

const Field3D FDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::FDDY(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dy;
}

/////////////////////////////////////////////////////////////////////////

const Field2D FDDZ(const Field2D &UNUSED(v), const Field2D &f, CELL_LOC outloc,
                   const std::string &UNUSED(method), REGION UNUSED(region)) {
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  return zeroFrom(f).setLocation(outloc);
}

const Field3D FDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc, const std::string &method, REGION region) {
  return bout::derivatives::index::FDDZ(v, f, outloc, method, region)
         / f.getCoordinates(outloc)->dz;
}
