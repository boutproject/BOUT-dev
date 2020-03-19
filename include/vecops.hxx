/*!************************************************************************
 * \file vecops.hxx
 *
 * Operators on vector objects
 * B.Dudson, October 2007
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

#ifndef __VECOPS_H__
#define __VECOPS_H__

class Field2D;
class Field3D;
class Vector2D;
class Vector3D;

#include "bout/deprecated.hxx"
#include "bout_types.hxx"
#include "bout/coordinates.hxx"
// Those are needed because we implement functions here.
// They can be dropped if we remove the deprecated wrappers.
#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"


/// Gradient of scalar field \p f, returning a covariant vector
///
/// All locations supported
///
/// @param[in] f  The field to differentiate
/// @param[in] outloc  The location where the result is desired
///                    By default this is the same location as the input \p f
/// @param[in] method  The method to use. The default is set in the options.
const Vector2D Grad(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                    const std::string& method = "DEFAULT");
const Vector3D Grad(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                    const std::string& method = "DEFAULT");

/// Perpendicular gradient of scalar field \p f
///
/// outloc must be either CELL_DEFAULT or f.getLocation() --> argument can be removed
///
/// result.x = df/dx - g_12/(JB)^2 df/dy
/// result.y = 0
/// result.z = df/dz - g_23/(JB)^2 df/dy
///
/// @param[in] f  The field to differentiate
/// @param[in] outloc  The cell location where the result is desired
/// @param[in] method  The method to use. The default is set in the options.
///
const Vector3D Grad_perp(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT");
inline const Vector3D Grad_perp(const Field3D& f, CELL_LOC outloc, DIFF_METHOD method) {
  return Grad_perp(f, outloc, toString(method));
}
const Vector2D Grad_perp(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT");

/// Divergence of a vector \p v, returning a scalar
///
/// All locations except `CELL_VSHIFT` supported. Note that if \p v is
/// at `CELL_VSHIFT`, then \p outloc must be `CELL_CENTRE`
///
/// @param[in] v  The vector to differentiate
/// @param[in] outloc  The cell location where the result is desired
/// @param[in] method  The method to use. The default is set in the options.
///
const Coordinates::metric_field_type Div(const Vector2D &v, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT");                                         
const Field3D Div(const Vector3D &v, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT");                  

const Coordinates::metric_field_type Div(const Vector2D &v, const Field2D &f,
		  CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT");

const Field3D Div(const Vector3D& v, const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT");
DEPRECATED(inline const Field3D Div(const Vector3D& v, const Field3D& f,
                                    const std::string& method,
                                    CELL_LOC outloc = CELL_DEFAULT)) {
  return Div(v, f, outloc, method);
}
inline const Field3D Div(const Vector3D& v, const Field3D& f, CELL_LOC outloc,
                         DIFF_METHOD method = DIFF_DEFAULT) {
  return Div(v, f, outloc, toString(method));
}
DEPRECATED(inline const Field3D Div(const Vector3D& v, const Field3D& f,
                                    DIFF_METHOD method, CELL_LOC outloc = CELL_DEFAULT)) {
  return Div(v, f, outloc, toString(method));
}

/// Curl of a vector
///
/// Does not currently support any output locations. \p v must not be
/// at `CELL_VSHIFT`
///
/// We can't support VSHIFT here as, e.g. DDY can't produce an output
/// at CELL_XLOW unless the input field is at CELL_XLOW, but then that
/// field will also be needed at CELL_YLOW, for example for another
/// component.
///
/// @param[in] v  The vector to differentiate
///
const Vector2D Curl(const Vector2D &v);
const Vector3D Curl(const Vector3D &v);

// Upwinding routines

/// Advection of a scalar field \p f by a velocity vector \p v
///
/// The vector and the field must be at the same location, which
/// cannot be CELL_VSHIFT
const Coordinates::metric_field_type V_dot_Grad(const Vector2D &v, const Field2D &f);
const Field3D V_dot_Grad(const Vector2D &v, const Field3D &f);
const Field3D V_dot_Grad(const Vector3D &v, const Field2D &f);
const Field3D V_dot_Grad(const Vector3D &v, const Field3D &f);

/// Advection of a vector field \p a by a velocity vector \p v
///
/// Both vectors must be at the same location, which cannot be CELL_VSHIFT
const Vector2D V_dot_Grad(const Vector2D &v, const Vector2D &a);
const Vector3D V_dot_Grad(const Vector2D &v, const Vector3D &a);
const Vector3D V_dot_Grad(const Vector3D &v, const Vector2D &a);
const Vector3D V_dot_Grad(const Vector3D &v, const Vector3D &a);

#endif // __VECOPS_H__
