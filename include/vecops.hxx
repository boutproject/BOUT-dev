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

#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

/// Gradient of scalar field \p f, returning a covariant vector
///
/// All locations supported
///
/// @param[in] f  The field to differentiate
/// @param[in] outloc The location where the result is desired (if staggered meshes are enabled)
///                   By default this is the same location as the input \p f
const Vector2D Grad(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT);
const Vector3D Grad(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT);

/// Gradient of scalar field \p f, returning a covariant vector
///
/// All locations supported
///
/// @param[in] f  The field to differentiate
/// @param[in] outloc_x  The cell location where the X component should be defined
/// @param[in] outloc_y  The cell location where the Y component should be defined
/// @param[in] outloc_z  The cell location where the Z component should be defined
const Vector3D DEPRECATED(Grad(const Field3D &f, 
			       CELL_LOC outloc_x, CELL_LOC outloc_y, CELL_LOC outloc_z));

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
///
const Vector3D Grad_perp(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT);

/// Perpendicular gradient of scalar field \p f
///
///
/// outloc must all be the same and must be either CELL_DEFAULT or f.getLocation() --> arguments can be removed
///
/// result.x = df/dx - g_12/(JB)^2 df/dy
/// result.y = 0
/// result.z = df/dz - g_23/(JB)^2 df/dy
/// 
/// @param[in] f  The field to differentiate
/// @param[in] outloc_x  The cell location where the X component should be defined
/// @param[in] outloc_y  The cell location where the Y component should be defined
/// @param[in] outloc_z  The cell location where the Z component should be defined
///
const Vector3D DEPRECATED(Grad_perp(const Field3D &f, 
				    CELL_LOC outloc_x, 
				    CELL_LOC outloc_y,
				    CELL_LOC outloc_z));

/// Divergence of a vector \p v, returning a scalar
///
/// All locations except CELL_VSHIFT supported
///
/// @param[in] v  The vector to differentiate
/// @param[in] outloc  The cell location where the result is desired
///
const Field2D Div(const Vector2D &v, CELL_LOC outloc = CELL_DEFAULT);
const Field3D Div(const Vector3D &v, CELL_LOC outloc = CELL_DEFAULT);

const Field2D Div(const Vector2D &v, const Field2D &f);
const Field3D Div(const Vector3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc = CELL_DEFAULT);
const Field3D Div(const Vector3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method = DIFF_DEFAULT);
const Field3D Div(const Vector3D &v, const Field3D &f);

/// Curl of a vector
///
/// All locations except CELL_VSHIFT supported
///
/// @param[in] v  The vector to differentiate
/// @param[in] outloc  The cell location where the result is desired
///
const Vector2D Curl(const Vector2D &v, CELL_LOC outloc = CELL_DEFAULT);
const Vector3D Curl(const Vector3D &v, CELL_LOC outloc = CELL_DEFAULT);
const Vector3D DEPRECATED(Curl(const Vector3D &v, 
			       CELL_LOC outloc_x, CELL_LOC outloc_y, CELL_LOC outloc_z));

// Upwinding routines

/// Advection of a scalar field \p f by a velocity vector \p v
///
/// The vector and the field must be at the same location, which
/// cannot be CELL_VSHIFT
const Field2D V_dot_Grad(const Vector2D &v, const Field2D &f);
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
