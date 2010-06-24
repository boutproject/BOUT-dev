/**************************************************************************
 *
 * Set variable boundaries
 * Currently just very simple options
 * Need to consider more complicated boundary conditions
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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

#ifndef __BOUNDARY_H__
#define __BOUNDARY_H__

#include "field3d.h"
#include "vector3d.h"

// boundary conditions
const int BNDRY_NONE         = 0;
const int BNDRY_ZERO         = 1; // Zero value
const int BNDRY_GRADIENT     = 2; // Zero gradient
const int BNDRY_LAPLACE      = 3; // This is laplace2 - laplacian zero everywhere
const int BNDRY_LAPLACE_GRAD = 4; // laplace - zero-gradient for last two points
const int BNDRY_ROTATE       = 5; // Rotates 180 degrees
const int BNDRY_ZAVERAGE     = 6; // Z average for last point
const int BNDRY_ROTATE_NEG   = 7; // Rotate 180 degrees, reverse sign
const int BNDRY_DIVCURL      = 8; // For vectors: Div B = Curl B = 0
const int BNDRY_RELAX_VAL    = 9; // Relax to given value
const int BNDRY_LAPLACE_ZERO = 10; // laplace - zero value on last point
const int BNDRY_LAPLACE_DECAY = 11; // Zero Laplacian, decaying solution
const int BNDRY_C_LAPLACE_DECAY = 12; // Const Laplacian, decaying


/// Prints boundary condition
void print_boundary(const char* fullname, const char* shortname);
void print_boundary(const char* name);
void print_boundary(const char* name, bool covariant);

/// Applies a boundary condition, depending on setting in BOUT.inp
void apply_boundary(Field2D &var, const char* fullname, const char* shortname);
void apply_boundary(Field3D &var, const char* fullname, const char* shortname);
void apply_boundary(Field2D &var, const char* name);
void apply_boundary(Field3D &var, const char* name);
void apply_boundary(Vector3D &var, const char* name);

// Inner x core boundary
void bndry_core_zero(Field2D &var);
void bndry_core_flat(Field2D &var);

void bndry_core_zero(Field3D &var);
void bndry_core_flat(Field3D &var);

void bndry_core_zero(Vector3D &var);
void bndry_core_flat(Vector3D &var);

void bndry_core_laplace(Field3D &var);
void bndry_core_laplace2(Field3D &var);

// Inner x PF boundar
void bndry_pf_zero(Field2D &var);
void bndry_pf_flat(Field2D &var);

void bndry_pf_zero(Field3D &var);
void bndry_pf_flat(Field3D &var);

void bndry_pf_zero(Vector3D &var);
void bndry_pf_flat(Vector3D &var);

void bndry_pf_laplace(Field3D &var);

// Inner x (core and PF) boundary
void bndry_inner_zero(Field2D &var);
void bndry_inner_flat(Field2D &var);

void bndry_inner_zero(Field3D &var);
void bndry_inner_flat(Field3D &var);

void bndry_inner_zero(Vector3D &var);
void bndry_inner_flat(Vector3D &var);

void bndry_inner_laplace(Field3D &var);

void bndry_inner_divcurl(Vector3D &var);

void bndry_inner_periodic(Field3D &var);

// Outer x boundary
void bndry_sol_zero(Field2D &var);
void bndry_sol_flat(Field2D &var);

void bndry_sol_zero(Field3D &var);
void bndry_sol_flat(Field3D &var);

void bndry_sol_zero(Vector3D &var);
void bndry_sol_flat(Vector3D &var);

void bndry_sol_laplace(Field3D &var);

void bndry_sol_divcurl(Vector3D &var);

void bndry_sol_periodic(Field3D &var);

// y boundary

void bndry_ydown_flat(Field2D &var);
void bndry_yup_flat(Field2D &var);

void bndry_ydown_zero(Field2D &var);
void bndry_yup_zero(Field2D &var);

void bndry_ydown_flat(Field3D &var);
void bndry_yup_flat(Field3D &var);

void bndry_ydown_zero(Field3D &var);
void bndry_yup_zero(Field3D &var);

void bndry_ydown_zero(Vector3D &var);
void bndry_yup_zero(Vector3D &var);

void bndry_ydown_rotate(Field3D &var, bool reverse = false);

void bndry_ydown_zaverage(Field3D &var);

// Relaxing boundaries
void bndry_core_relax_val(Field3D &F_var, const Field3D &var, BoutReal value, BoutReal rate=10.);
void bndry_pf_relax_val(Field3D &F_var, const Field3D &var, BoutReal value, BoutReal rate=10.);
void bndry_inner_relax_val(Field3D &F_var, const Field3D &var, BoutReal value, BoutReal rate=10.);
void bndry_inner_relax_val2(Field3D &F_var, const Field3D &var, BoutReal value, BoutReal rate=10.);
void bndry_sol_relax_val(Field3D &F_var, const Field3D &var, BoutReal value, BoutReal rate=10.);
void bndry_sol_relax_val2(Field3D &F_var, const Field3D &var, BoutReal value, BoutReal rate=10.);
void bndry_ydown_relax_val(Field3D &F_var, const Field3D &var, BoutReal value, BoutReal rate=10.);
void bndry_yup_relax_val(Field3D &F_var, const Field3D &var, BoutReal value, BoutReal rate=10.);

// Relax to zero gradient
void bndry_core_relax_flat(Field3D &F_var, const Field3D &var, BoutReal rate=10.);
void bndry_pf_relax_flat(Field3D &F_var, const Field3D &var, BoutReal rate=10.);
void bndry_inner_relax_flat(Field3D &F_var, const Field3D &var, BoutReal rate=10.);
void bndry_sol_relax_flat(Field3D &F_var, const Field3D &var, BoutReal rate=10.);
void bndry_ydown_relax_flat(Field3D &F_var, const Field3D &var, BoutReal rate=10.);
void bndry_yup_relax_flat(Field3D &F_var, const Field3D &var, BoutReal rate=10.);

// Symmetric boundary condition
void bndry_core_sym(Field3D &var);
void bndry_pf_sym(Field3D &var);
void bndry_inner_sym(Field3D &var);
void bndry_sol_sym(Field3D &var);
void bndry_ydown_sym(Field3D &var);
void bndry_yup_sym(Field3D &var);

// Relax to symmetric conditions
void bndry_core_relax_sym(Field3D &F_var, const Field3D &var, BoutReal rate=10.);
void bndry_pf_relax_sym(Field3D &F_var, const Field3D &var, BoutReal rate=10.);
void bndry_inner_relax_sym(Field3D &F_var, const Field3D &var, BoutReal rate=10.);
void bndry_sol_relax_sym(Field3D &F_var, const Field3D &var1, BoutReal rate=10.);
void bndry_ydown_relax_sym(Field3D &F_var, const Field3D &var, BoutReal rate=10.);
void bndry_yup_relax_sym(Field3D &F_var, const Field3D &var, BoutReal rate=10.);

// z boundary

void bndry_toroidal(Field3D &var);
void bndry_toroidal(Vector3D &var);

// Zero laplace 
void bndry_inner_zero_laplace(Field3D &var);
void bndry_core_zero_laplace(Field3D &var);
void bndry_pf_zero_laplace(Field3D &var);

// Zero Laplacian, decaying solution
void bndry_inner_laplace_decay(Field3D &var);
void bndry_outer_laplace_decay(Field3D &var);

// Constant Laplacian, decaying solution
void bndry_inner_const_laplace_decay(Field3D &var);
void bndry_outer_const_laplace_decay(Field3D &var);

#endif // __BOUNDARY_H__
