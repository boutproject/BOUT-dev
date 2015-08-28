/**************************************************************
 * Gyro-averaging operators
 *
 *
 * 2010-09-03 Ben Dudson <bd512@york.ac.uk>
 *    * Initial version, simple averaging operator
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
 **************************************************************/

#include <globals.hxx>
#include <difops.hxx>
#include <gyro_average.hxx>
#include <invert_laplace.hxx>

/// Approximate G(f) = f + rho^2*Delp2(f) using Taylor expansion
const Field3D gyroTaylor0(const Field3D &f, const Field3D &rho) {
  return f + SQ(rho) * Delp2(f);
}

/// Pade approximation G_0 = (1 - rho^2*Delp2)g = f
const Field3D gyroPade0(const Field3D &f, BoutReal rho, int flags) {
  /// Have to use Z average of rho for efficient inversion

  Field2D a = 1.0;
  Field2D d = -rho*rho;

  /// Invert, leaving boundaries unchanged
  return invert_laplace(f, flags, &a, NULL, &d);
}

/// Pade approximation G_0 = (1 - rho^2*Delp2)g = f
const Field3D gyroPade0(const Field3D &f, const Field2D &rho, int flags) {
  /// Have to use Z average of rho for efficient inversion

  Field2D a = 1.0;
  Field2D d = -rho*rho;

  /// Invert, leaving boundaries unchanged
  return invert_laplace(f, flags, &a, NULL, &d);
}

/// Pade approximation G_0 = (1 - rho^2*Delp2)g = f
const Field3D gyroPade0(const Field3D &f, const Field3D &rho, int flags) {
  /// Have to use Z average of rho for efficient inversion
  return gyroPade0(f, rho.DC(), flags);
}

/// Pade approximation G_1 = (1 - 0.5*rho^2*Delp2)g = f
const Field3D gyroPade1(const Field3D &f, BoutReal rho, int flags) {
  Field2D a = 1.0;
  Field2D d = -0.5*rho*rho;

  /// Invert, leaving boundaries unchanged
  return invert_laplace(f, flags, &a, NULL, &d);
}

/// Pade approximation G_1 = (1 - 0.5*rho^2*Delp2)g = f
const Field3D gyroPade1(const Field3D &f, const Field2D &rho, int flags) {
  Field2D a = 1.0;
  Field2D d = -0.5*rho*rho;

  /// Invert, leaving boundaries unchanged
  return invert_laplace(f, flags, &a, NULL, &d);
}

/// Pade approximation G_1 = (1 - 0.5*rho^2*Delp2)g = f
const Field3D gyroPade1(const Field3D &f, const Field3D &rho, int flags) {
  /// Have to use Z average of rho for efficient inversion
  return gyroPade1(f, rho.DC(), flags);
}

const Field2D gyroPade1(const Field2D &f, const Field2D &rho, int flags) {
  // Very inefficient implementation
  Field3D tmp = f;
  tmp = gyroPade1(tmp, rho, flags);
  return tmp.DC();
}

/// Pade approximation G_2 = (1 - 0.5*rho^2*Delp2)g = f
const Field3D gyroPade2(const Field3D &f, BoutReal rho, int flags) {
  Field3D result = gyroPade1(gyroPade1(f, rho, flags), rho, flags);
  mesh->communicate(result);
  result = 0.5*rho*rho*Delp2( result );
  result.applyBoundary("dirichlet");
  return result;
}

/// Pade approximation G_2 = (1 - 0.5*rho^2*Delp2)g = f
const Field3D gyroPade2(const Field3D &f, const Field2D &rho, int flags) {
  Field3D result = gyroPade1(gyroPade1(f, rho, flags), rho, flags);
  mesh->communicate(result);
  result = 0.5*rho*rho*Delp2( result );
  result.applyBoundary("dirichlet");
  return result;
}

/// Pade approximation G_2 = (1 - 0.5*rho^2*Delp2)g = f
const Field3D gyroPade2(const Field3D &f, const Field3D &rho, int flags) {
  /// Have to use Z average of rho for efficient inversion
  return gyroPade2(f, rho.DC(), flags);
}

