/*!************************************************************
 * \file gyro_average.hxx
 * 
 * Gyro-averaging operators
 *
 *
 * 2010-09-03 Ben Dudson <bd512@york.ac.uk>
 *    * Initial version, simple averaging operator
 * 
 **************************************************************
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

#ifndef __GYRO_AVERAGE_H__
#define __GYRO_AVERAGE_H__

#include "field3d.hxx"
#include "invert_laplace.hxx"

/// INVERT_BNDRY_ONE | INVERT_IN_RHS | INVERT_OUT_RHS; uses old-style
/// Laplacian inversion flags
constexpr int GYRO_FLAGS = INVERT_BNDRY_ONE + INVERT_RHS;

/// Gyro-average using Taylor series approximation
///
///     \f$ \Gamma(f) = f + \rho^2 \nabla_\perp^2(f)\f$
///
/// Note: Faster, but less robust than Pade approximations
///
/// @param[in]  f  The field to gyro-average
/// @param[in] rho Gyro-radius
Field3D gyroTaylor0(const Field3D& f, const Field3D& rho);

/// Gyro-average using Pade approximation
///
///     \f$ \Gamma_0 = (1 - \rho^2 \nabla_\perp^2)g = f\f$
///
/// NOTE: Uses Z average of rho for efficient inversion
///
/// @param[in] f   The field to gyro-average
/// @param[in] rho  Gyro-radius
/// @param[in] flags  Flags to be passed to the Laplacian inversion operator
Field3D gyroPade0(const Field3D& f, const Field3D& rho, int inner_boundary_flags, int outer_boundary_flags);
Field3D gyroPade0(const Field3D& f, const Field2D& rho, int inner_boundary_flags, int outer_boundary_flags);
Field3D gyroPade0(const Field3D& f, BoutReal rho, int inner_boundary_flags, int outer_boundary_flags);

// Can replace these with default arguments to versions above once the deprecated versions
// below are removed
inline Field3D gyroPade0(const Field3D& f, const Field3D& rho) {
  return gyroPade0(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}
inline Field3D gyroPade0(const Field3D& f, const Field2D& rho) {
  return gyroPade0(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}
inline Field3D gyroPade0(const Field3D& f, BoutReal rho) {
  return gyroPade0(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}

[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field3D gyroPade0(const Field3D& f, const Field3D& rho, int flags);
[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field3D gyroPade0(const Field3D& f, const Field2D& rho, int flags);
[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field3D gyroPade0(const Field3D& f, BoutReal rho, int flags);

/// Pade approximation \f$Gamma_1 = (1 - \frac{1}{2} \rho^2 \nabla_\perp^2)g = f\f$
///
/// Note: Have to use Z average of rho for efficient inversion
/// 
/// @param[in] f   The field to gyro-average
/// @param[in] rho  Gyro-radius
/// @param[in] flags  Flags to be passed to the Laplacian inversion operator
Field3D gyroPade1(const Field3D& f, const Field3D& rho, int inner_boundary_flags, int outer_boundary_flags);
Field3D gyroPade1(const Field3D& f, const Field2D& rho, int inner_boundary_flags, int outer_boundary_flags);
Field3D gyroPade1(const Field3D& f, BoutReal rho, int inner_boundary_flags, int outer_boundary_flags);
Field2D gyroPade1(const Field2D& f, const Field2D& rho, int inner_boundary_flags, int outer_boundary_flags);

// Can replace these with default arguments to versions above once the deprecated versions
// below are removed
inline Field3D gyroPade1(const Field3D& f, const Field3D& rho) {
  return gyroPade1(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}
inline Field3D gyroPade1(const Field3D& f, const Field2D& rho) {
  return gyroPade1(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}
inline Field3D gyroPade1(const Field3D& f, BoutReal rho) {
  return gyroPade1(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}
inline Field2D gyroPade1(const Field2D& f, const Field2D& rho) {
  return gyroPade1(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}

[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field3D gyroPade1(const Field3D& f, const Field3D& rho, int flags);
[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field3D gyroPade1(const Field3D& f, const Field2D& rho, int flags);
[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field3D gyroPade1(const Field3D& f, BoutReal rho, int flags);
[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field2D gyroPade1(const Field2D& f, const Field2D& rho, int flags);

/// Pade approximation 
///
/// \f[
///    \Gamma_2(f) = \frac{1}{2}\rho^2 \nabla_\perp^2 ( 1 - \frac{1}{2} \rho^2 \nabla_\perp^2)^{-1}\Gamma_1(f)
/// \f]
///
/// Note: Have to use Z average of rho for efficient inversion
///
/// @param[in] f   The field to gyro-average
/// @param[in] rho  Gyro-radius
/// @param[in] flags  Flags to be passed to the Laplacian inversion operator
Field3D gyroPade2(const Field3D& f, const Field3D& rho, int inner_boundary_flags, int outer_boundary_flags);
Field3D gyroPade2(const Field3D& f, const Field2D& rho, int inner_boundary_flags, int outer_boundary_flags);
Field3D gyroPade2(const Field3D& f, BoutReal rho, int inner_boundary_flags, int outer_boundary_flags);

// Can replace these with default arguments to versions above once the deprecated versions
// below are removed
inline Field3D gyroPade2(const Field3D& f, const Field3D& rho) {
  return gyroPade2(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}
inline Field3D gyroPade2(const Field3D& f, const Field2D& rho) {
  return gyroPade2(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}
inline Field3D gyroPade2(const Field3D& f, BoutReal rho) {
  return gyroPade2(f, rho, GYRO_FLAGS, GYRO_FLAGS);
}

[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field3D gyroPade2(const Field3D& f, const Field3D& rho, int flags);
[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field3D gyroPade2(const Field3D& f, const Field2D& rho, int flags);
[[gnu::deprecated("Please use version with separate inner_boundary_flags and outer_boundary_flags")]]
Field3D gyroPade2(const Field3D& f, BoutReal rho, int flags);

#endif // __GYRO_AVERAGE_H__
