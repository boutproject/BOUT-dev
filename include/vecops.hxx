/**************************************************************************
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

const Vector2D Grad(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT);
const Vector3D Grad(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT);
const Vector3D Grad(const Field3D &f, 
                    CELL_LOC outloc_x, CELL_LOC outloc_y, CELL_LOC outloc_z = CELL_DEFAULT);

const Vector3D Grad_perp(const Field3D &f, 
			 CELL_LOC outloc_x = CELL_DEFAULT, 
			 CELL_LOC outloc_y = CELL_DEFAULT,
			 CELL_LOC outloc_z = CELL_DEFAULT);

const Field2D Div(const Vector2D &v, CELL_LOC outloc = CELL_DEFAULT);
const Field3D Div(const Vector3D &v, CELL_LOC outloc = CELL_DEFAULT);

const Field2D Div(const Vector2D &v, const Field2D &f);
const Field3D Div(const Vector3D &v, const Field3D &f, DIFF_METHOD method, CELL_LOC outloc = CELL_DEFAULT);
const Field3D Div(const Vector3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method = DIFF_DEFAULT);
const Field3D Div(const Vector3D &v, const Field3D &f);

const Vector2D Curl(const Vector2D &v, CELL_LOC outloc = CELL_DEFAULT);
const Vector3D Curl(const Vector3D &v, CELL_LOC outloc = CELL_DEFAULT);
const Vector3D Curl(const Vector3D &v, 
                    CELL_LOC outloc_x, CELL_LOC outloc_y, CELL_LOC outloc_z);

// Upwinding routines

const Field2D V_dot_Grad(const Vector2D &v, const Field2D &f);
const Field3D V_dot_Grad(const Vector2D &v, const Field3D &f);
const Field3D V_dot_Grad(const Vector3D &v, const Field2D &f);
const Field3D V_dot_Grad(const Vector3D &v, const Field3D &f);

const Vector2D V_dot_Grad(const Vector2D &v, const Vector2D &a);
const Vector3D V_dot_Grad(const Vector2D &v, const Vector3D &a);
const Vector3D V_dot_Grad(const Vector3D &v, const Vector2D &a);
const Vector3D V_dot_Grad(const Vector3D &v, const Vector3D &a);

#endif // __VECOPS_H__
