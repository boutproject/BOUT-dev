/**************************************************************************
 * Perpendicular Laplacian inversion using FFT and Tridiagonal solver
 *
 * Equation solved is: d*\nabla^2_\perp x + (1/c)\nabla_perp c\cdot\nabla_\perp x + a x = b
 * 
 * Where a, c and d are functions of x and y only (not z)
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

#ifndef __LAPLACE_H__
#define __LAPLACE_H__

#ifdef BOUT_HAS_PETSC
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE
#endif

#include "fieldperp.h"
#include "field3d.h"
#include "field2d.h"

#include "dcomplex.h"

// Inversion flags

const int INVERT_DC_IN_GRAD  = 1;
const int INVERT_AC_IN_GRAD  = 2;
const int INVERT_DC_OUT_GRAD = 4;
const int INVERT_AC_OUT_GRAD = 8;
const int INVERT_ZERO_DC     = 16;
const int INVERT_START_NEW   = 32;
const int INVERT_BNDRY_ONE   = 64; // Sets the width of the boundary to 1
const int INVERT_4TH_ORDER   = 128; // Use band solver for 4th order in x

const int INVERT_AC_IN_LAP   = 256;
const int INVERT_AC_OUT_LAP  = 512;

const int INVERT_IN_SYM  =  1024; // Use symmetry to enforce either zero-value or zero-gradient
const int INVERT_OUT_SYM =  2048; // Same for outer boundary
const int INVERT_IN_SET  =  4096; // Set inner boundary
const int INVERT_OUT_SET =  8192; // Set outer boundary
const int INVERT_IN_RHS  = 16384; // Use input value in RHS at inner boundary
const int INVERT_OUT_RHS = 32768; // Use input value in RHS at outer boundary
const int INVERT_KX_ZERO = 65536; // Zero the kx=0, n = 0 component

int invert_init();

void laplace_tridag_coefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b, dcomplex &c, const Field2D *ccoef = NULL, const Field2D *d=NULL);

int invert_laplace(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *c=NULL, const Field2D *d=NULL);
int invert_laplace(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c=NULL, const Field2D *d=NULL);

/// More readable API for calling Laplacian inversion. Returns x
const Field3D invert_laplace(const Field3D &b, int flags, 
                             const Field2D *a = NULL, const Field2D *c=NULL, const Field2D *d=NULL);

#endif // __LAPLACE_H__

