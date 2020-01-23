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
#include <interpolation_xz.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <unused.hxx>

void printLocation(const Field3D& var) {
  output << toString(var.getLocation());
}
void printLocation(const Field2D& var) {
  output << toString(var.getLocation());
}

const char* strLocation(CELL_LOC loc) { return toString(loc).c_str(); }

const Field3D interpolate(const Field3D &f, const Field3D &delta_x,
                          const Field3D &delta_z) {
  TRACE("Interpolating 3D field");
  XZLagrange4pt interpolateMethod{f.getMesh()};
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
  Field3D result{emptyFrom(delta_x)};

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

void XZInterpolationFactory::ensureRegistered() {}

namespace {
RegisterXZInterpolation<XZHermiteSpline> registerinterphermitespline{"hermitespline"};
RegisterXZInterpolation<XZMonotonicHermiteSpline> registerinterpmonotonichermitespline{
    "monotonichermitespline"};
RegisterInterpolation<HermiteSplineOnlyZ> registerinterphermitesplineonlyz{"hermitesplineonlyz"};
RegisterXZInterpolation<XZLagrange4pt> registerinterplagrange4pt{"lagrange4pt"};
RegisterXZInterpolation<XZBilinear> registerinterpbilinear{"bilinear"};
} // namespace

constexpr decltype(XZInterpolationFactory::type_name) XZInterpolationFactory::type_name;
constexpr decltype(XZInterpolationFactory::section_name) XZInterpolationFactory::section_name;
constexpr decltype(XZInterpolationFactory::option_name) XZInterpolationFactory::option_name;
constexpr decltype(XZInterpolationFactory::default_type) XZInterpolationFactory::default_type;
